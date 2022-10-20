# collection of functions for parsing PrairieView metadata from imaging and alloptical experiments
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass

# need to implement alternative class constructor using cls.method
from typing import List, Union, Any
from imagingplus.main.subcore import ImagingMetadata


# class ImagingMetadata:
#     """Class containing metadata about imaging microscope"""
#
#     def __init__(self, microscope, n_frames, fps, frame_x, frame_y, n_planes, pix_sz_x, pix_sz_y, **kwargs):
#
#         self.microscope = microscope  #: given name of microscope (use 'Bruker' to process Metadata from PrairieView)
#         self.n_frames = n_frames  # number of imaging frames in the current trial
#         self.fps = fps  # rate of imaging acquisition (frames per second)
#         self.frame_x = frame_x  # num of pixels in the x direction of a single frame
#         self.frame_y = frame_y  # num of pixels in the y direction of a single frame
#         self.n_planes = n_planes  # num of FOV planes in imaging acquisition
#         self.pix_sz_x = pix_sz_x  # size of a single imaging pixel in x direction (microns)
#         self.pix_sz_y = pix_sz_y  # size of a single imaging pixel in y direction (microns)
#         for key, value in kwargs.items():
#             setattr(self, key, value)
#
#     def __repr__(self):
#         return f'ImagingMetadata for imaging cellsdata collected with {self.microscope}.'


class PrairieViewMetadata(ImagingMetadata):
    """class for parsing metadata of imaging microscope system."""

    def __init__(self, pv_xml_dir: str, microscope: str = 'Bruker'):

        print(f'\n\- Adding Imaging Acquisition Metadata from {microscope} ...')
        self.pv_xml_dir = pv_xml_dir  #: path to the directory containing the .xml imaging output from PrairieView for a given trial
        pv_metadata = self._parsePVMetadata()
        super().__init__(microscope=microscope, **pv_metadata)  # call to ImagingMetadata parent class

        self.scan_x = pv_metadata['scan_x']  #: resonant scan center in x axis
        self.scan_y = pv_metadata['scan_y']  #: resonant scan center in y axis
        self.zoom: float = pv_metadata['zoom']  #: zoom level on Bruker microscope
        self.last_good_frame = pv_metadata[
            'last_good_frame']  #: indicates when the last good frame was during the t-series recording, if nothing was wrong the value is 0, otherwise it is >0 and that indicates that PV is not sure what happened after the frame listed, but it could be corrupt cellsdata
        # for key, value in kwargs.items():  # todo removing kwargs
        #     setattr(self, key, value)

    def __repr__(self):
        return f'PrairieViewMetadata from: {self.pv_xml_dir}\n\tn planes: {self.n_planes} \n\tn frames: {self.n_frames}' \
               f'\n\tfps: {self.fps} \n\tframe size (px): {self.frame_x} x {self.frame_y}  \n\tzoom: {self.zoom} \n\t' \
               f'pixel size (um): {self.pix_sz_x}, {self.pix_sz_y} ' \
               f'\n\tscan centre (V): {self.scan_x}, {self.scan_y}'

    @staticmethod
    def _getPVStateShard(root, key):
        """
        Find the value, description and indices of a particular parameter from an xml file
        """
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
            raise Exception(f'ERROR: no element or subelement with the key {key} under root: {root}')

        return value, description, index

    def _parsePVMetadata(self) -> dict:
        """
        Parse all of the relevant acquisition metadata from the PrairieView xml file for this recording

        :return: dict; parsed values related to metadata from PrairieView xml file.
        """

        print('\n\t\- Parsing PV Metadata for Bruker microscope...')

        xml_dir_path = self.pv_xml_dir  # starting path to search for the .xml PrairieView files
        xml_path = []  # searching for xml path
        print(f'\t searching for xml path in tiff path directory at: {xml_dir_path} ... ')
        try:  # look for xml file in path, or two paths up (achieved by decreasing count in while loop)
            count = 2
            while count != 0 and not xml_path:
                count -= 1
                for file in os.listdir(xml_dir_path):
                    if file.endswith('.xml'):
                        xml_path = os.path.join(xml_dir_path, file)
                xml_dir_path = os.path.dirname(xml_path)  # re-assign xml_dir_path as next folder up
            assert os.path.exists(xml_path), f'xml_path could not be found: {xml_path}'
        except Exception:
            raise Exception(f'ERROR: Could not find xml for this acquisition at: {xml_dir_path}')

        xml_tree = ET.parse(xml_path)  # parse xml from a path
        root = xml_tree.getroot()  # make xml tree structure

        sequence = root.find('Sequence')
        acq_type = sequence.get('type')

        if 'ZSeries' in acq_type:
            n_planes = len(sequence.findall('Frame'))
        else:
            n_planes = 1

        frame_branch = root.findall('Sequence/Frame')[-1]
        try:
            frame_period = float(self._getPVStateShard(frame_branch, 'framePeriod')[0])
        except Exception:
            frame_period = float(self._getPVStateShard(root, 'framePeriod')[0])
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

        self.fps: float = fps / n_planes
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
              '\n\tframe size (px):', self.frame_x, 'y', self.frame_y,
              '\n\tzoom:', self.zoom,
              '\n\tpixel size (um):', self.pix_sz_x, self.pix_sz_y,
              '\n\tscan centre (V):', self.scan_x, self.scan_y
              )

        return_dict = {'fps': round(self.fps, 5),
                       'n_planes': self.n_planes,
                       'n_frames': self.n_frames,
                       'frame_x': self.frame_x,
                       'frame_y': self.frame_y,
                       'zoom': self.zoom,
                       'pix_sz_x': round(pix_sz_x, 5),
                       'pix_sz_y': round(pix_sz_y, 5),
                       'scan_x': scan_x,
                       'scan_y': scan_y,
                       'last_good_frame': last_good_frame}

        return return_dict

