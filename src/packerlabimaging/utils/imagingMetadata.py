# collection of functions for parsing PrairieView metadata from imaging and alloptical experiments
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass

@dataclass
class PrairieViewMetadata:
    tiff_path_dir: str  # s2pResultsPath to the directory containing the .tiff imaging output from PrairieView. the .xml PrairieView files should also be in the same directory.
    def __post_init__(self):
        self.n_frames: int = 0  # number of imaging frames in the current trial
        self.fps = None  # rate of imaging acquisition (frames per second)
        self.frame_x = None  # num of pixels in the x direction of a single frame
        self.frame_y = None  # num of pixels in the y direction of a single frame
        self.n_planes = None  # num of FOV planes in imaging acquisition
        self.pix_sz_x = None  # size of a single imaging pixel in x direction (microns)
        self.pix_sz_y = None  # size of a single imaging pixel in y direction (microns)
        self.scan_x = None  # TODO ROB - not sure what the comment for this is
        self.scan_y = None  # TODO ROB - not sure what the comment for this is
        self.zoom: float = 0.0 # zoom level on Bruker microscope
        self.last_good_frame = None  # indicates when the last good frame was during the t-series recording, if nothing was wrong the value is 0, otherwise it is >0 and that indicates that PV is not sure what happened after the frame listed, but it could be corrupt data

        self._parsePVMetadata()

    def __repr__(self):
        return(f'PrairieViewMetadata from: {self.tiff_path_dir}\n\tn planes: {self.n_planes} \n\tn frames: {self.n_frames} '
               f'\n\tfps: {self.fps} \n\tframe size (px): {self.frame_x} x {self.frame_y}  \n\tzoom: {self.zoom} \n\t pixel size (um): {self.pix_sz_x}, {self.pix_sz_y} '
               f'\n\tscan centre (V):', self.scan_x, self.scan_y)

    @staticmethod
    def _getPVStateShard(root, key):
        '''
        Find the value, description and indices of a particular parameter from an xml file

        Inputs:
            s2pResultsPath        - s2pResultsPath to xml file
            key         - string corresponding to key in xml tree
        Outputs:
            value       - value of the key
            description - unused
            index       - index that the key was found at
        '''
        value = []
        description = []
        index = []

        pv_state_shard = root.find('PVStateShard')  # find pv state shard element in root

        for elem in pv_state_shard:  # for each element in pv state shard, find the value for the specified key
            if elem.get('key') == key:
                if len(elem) == 0:  # if the element has only one subelement
                    value = elem.get('value')
                    break

                else:  # if the element has many subelements (i.e. lots of entries for that key)
                    for subelem in elem:
                        value.append(subelem.get('value'))
                        description.append(subelem.get('description'))
                        index.append(subelem.get('index'))
            else:
                for subelem in elem:  # if key not in element, try subelements
                    if subelem.get('key') == key:
                        value = elem.get('value')
                        break

            if value:  # if found key in subelement, break the loop
                break

        if not value:  # if no value found at all, raise exception
            raise Exception('ERROR: no element or subelement with that key')

        return value, description, index

    def _parsePVMetadata(self):
        '''
        Parse all of the relevant acquisition metadata from the PrairieView xml file for this recording

        '''

        print('\n\----- Parsing PV Metadata for Bruker microscope...')


        tiff_path = self.tiff_path_dir  # starting s2pResultsPath to search for the .xml PrairieView files
        xml_path = []  # searching for xml s2pResultsPath

        try:  # look for xml file in s2pResultsPath, or two paths up (achieved by decreasing count in while loop)
            count = 2
            while count != 0 and not xml_path:
                count -= 1
                for file in os.listdir(tiff_path):
                    if file.endswith('.xml'):
                        xml_path = os.path.join(tiff_path, file)
                tiff_path = os.path.dirname(tiff_path)  # re-assign tiff_path as next folder up

        except:
            raise Exception('ERROR: Could not find xml for this acquisition, check it exists')

        xml_tree = ET.parse(xml_path)  # parse xml from a s2pResultsPath
        root = xml_tree.getroot()  # make xml tree structure

        sequence = root.find('Sequence')
        acq_type = sequence.get('type')

        if 'ZSeries' in acq_type:
            n_planes = len(sequence.findall('Frame'))
        else:
            n_planes = 1

        frame_branch = root.findall('Sequence/Frame')[-1]
        #         frame_period = float(self._getPVStateShard(root,'framePeriod')[0])
        frame_period = float(self._getPVStateShard(frame_branch, 'framePeriod')[0])
        fps = 1 / frame_period

        frame_x = int(self._getPVStateShard(root, 'pixelsPerLine')[0])
        frame_y = int(self._getPVStateShard(root, 'linesPerFrame')[0])
        zoom = float(self._getPVStateShard(root, 'opticalZoom')[0])

        scan_volts, _, index = self._getPVStateShard(root, 'currentScanCenter')
        for scan_volts, index in zip(scan_volts, index):
            if index == 'XAxis':
                scan_x = float(scan_volts)
            if index == 'YAxis':
                scan_y = float(scan_volts)

        pixel_size, _, index = self._getPVStateShard(root, 'micronsPerPixel')
        for pixel_size, index in zip(pixel_size, index):
            if index == 'XAxis':
                pix_sz_x = float(pixel_size)
            if index == 'YAxis':
                pix_sz_y = float(pixel_size)

        if n_planes == 1:
            n_frames = root.findall('Sequence/Frame')[-1].get('index')  # use suite2p output instead later
        else:
            n_frames = root.findall('Sequence')[-1].get('cycle')

        extra_params = root.find('Sequence/Frame/ExtraParameters')
        last_good_frame = extra_params.get('lastGoodFrame')


        self.fps = fps / n_planes
        self.frame_x = frame_x
        self.frame_y = frame_y
        self.n_planes = n_planes
        self.pix_sz_x = pix_sz_x
        self.pix_sz_y = pix_sz_y
        self.scan_x = scan_x
        self.scan_y = scan_y
        self.zoom = zoom
        self.n_frames = int(n_frames)
        self.last_good_frame = last_good_frame


        print('\tn planes:', self.n_planes,
              '\n\tn frames:', self.n_frames,
              '\n\tfps:', self.fps,
              '\n\tframe size (px):', self.frame_x, 'x', self.frame_y,
              '\n\tzoom:', self.zoom,
              '\n\tpixel size (um):', self.pix_sz_x, self.pix_sz_y,
              '\n\tscan centre (V):', self.scan_x, self.scan_y
              )


        # return {
        # "fps": fps / n_planes,
        # "frame_x": frame_x,
        # "frame_y": frame_y,
        # "n_planes": n_planes,
        # "pix_sz_x": pix_sz_x,
        # "pix_sz_y": pix_sz_y,
        # "scan_x": scan_x,
        # "scan_y": scan_y,
        # "zoom": zoom,
        # "n_frames": int(n_frames),
        # "last_good_frame": last_good_frame
        # }

@dataclass
class ImagingMetadata:
    n_frames: int  # number of imaging frames in the current trial
    fps: str    # rate of imaging acquisition (frames per second)
    frame_x: int  # num of pixels in the x direction of a single frame
    frame_y: int  # num of pixels in the y direction of a single frame
    n_planes: int  # num of FOV planes in imaging acquisition
    pix_sz_x: float  # size of a single imaging pixel in x direction (microns)
    pix_sz_y: float  # size of a single imaging pixel in y direction (microns)