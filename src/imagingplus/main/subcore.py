from __future__ import absolute_import
from dataclasses import dataclass, field
from typing import Optional, MutableMapping, Union, TypedDict, List, Dict, Any

import numpy as np
import pandas as pd


@dataclass
class TemporalData:
    """1-D time series datasets corresponding with an imaging trial. Accomodates all 1-D series that are temporally synchronized to each other."""
    file_path: str  #: path to cellsdata file
    sampling_rate: float  #: rate of cellsdata collection (Hz)
    channels: List[str]  #: list of cellsdata channel names.
    data: pd.DataFrame  #: N columns x Time array of N x 1D cellsdata channels collected over Time at the same sampling rate
    units: List[str] = None  #: list of units corresponding to channels.
    frame_times: Union[
        list, np.ndarray] = None  #: timestamps representing imaging frame times. must be of same time duration as imaging dataset.
    sparse_data: pd.DataFrame = None  #: dataframe of timeseries channels containing cellsdata keyed at imaging frames for each timeseries channel
    crop_offset_time: float = 0.0  #: length of time (secs) of the offset after cropping temporal data

    def __post_init__(self):
        if self.units:
            if len(self.units) != len(self.channels): raise AttributeError(f'Units must be provided for all channels.')
        print(f"Created new TemporalData of {self.n_channels} x {self.n_timepoints} (sampled at {self.sampling_rate}")
        pass

    def __repr__(self):
        return f"TemporalData module, containing {self.n_timepoints} timepoints and {self.n_channels} channels."

    def __str__(self):
        return f"TemporalData module, containing {self.n_timepoints} timepoints and {self.n_channels} channels: \n{self.channels}"

    @property
    def n_frames(self):
        """number of timed imaging frames
        """
        return len(self.frame_times)

    @property
    def n_timepoints(self):
        """number of cellsdata collection timepoints
        """
        return self.data.shape[0]

    @property
    def n_channels(self):
        """number of cellsdata collection channels
        """
        return self.data.shape[1]


    def get_sparse_data(self, frame_times: Union[list, np.ndarray] = None) -> pd.DataFrame:
        """
        Returns dictionary of numpy array keyed on channels from paq_data timed to 2photon imaging frame_times.
        In effect this works as a downsampling algorithm.

        :param frame_times: imaging frame timestamps to use for collecting sparse temporal data
        :return: Dataframe of sparse temporal data across all temporal data channels collected
        """

        assert hasattr(self,
                       'frame_times') or frame_times, 'no frame_times given to retrieve temporal sync. data from those timestamps.'

        frame_times = self.frame_times if frame_times is None else frame_times

        print(f"\n\t\- Getting imaging key frames timed cellsdata from {len(frame_times)} frames ... ")

        # read in and save sparse version of all tmdata channels (only save tmdata from timepoints at frame clock times)
        sparse_data = {}
        for idx, chan in enumerate(self.channels):
            print(f'\t\t\- Adding sparse tmdata for channel: {chan} ')
            try:
                data = self.data.loc[frame_times, chan]
            except ValueError:
                # assert self.crop_offset_time == 0, 'not able to collect sparse data from cropped data.'
                frame_idxs = np.searchsorted(self.data.index.to_numpy(), self.frame_times)
                data = self.data.iloc[frame_idxs, idx]
            sparse_data[chan] = data

        sparse_data = pd.DataFrame(sparse_data)
        print(f"\t|- Collected sparse data:  {sparse_data.shape} ... ")

        return sparse_data

    def cropData(self, begin: int, end: int, channels: Union[str, List[str]] = 'all', replace: bool = False):
        """
        Crops saved temporal cellsdata channels to the timestamps of begin and end.

        :param channels: channels to crop cellsdata, default is 'all' (all channels in dataset).
        :param begin: timestamp to begin cropping at
        :param end: timestamp to end cropping at
        :param replace: if True, then replace .cellsdata attribute with newly cropped cellsdata. if False, return the cropped dataframe.
        """

        channels = self.channels if channels == 'all' else channels
        print(f"\- cropping {channels} to {begin} and {end} paq clock times.")

        assert self.data.shape[0] >= (
                end - begin), f'Not enough time series samples in cellsdata to crop between the provided times.'
        data_ = self.data.loc[begin: end, channels]

        if replace:
            self.data = data_
        else:
            return data_

        self.crop_offset_time = begin / self.sampling_rate

        # older approach...
        # for channel in channels:
        #     print(f"\- cropping {channel} to {begin} and {end} paq clock times.")
        #     cellsdata = self.cellsdata[channel].to_numpy()
        #     assert len(cellsdata) >= (end - begin), f'{channel} paq cellsdata is not long enough to crop between the provided clock times.'
        #     cropdata = cellsdata[begin: end]
        #     setattr(self, channel, cropdata)


@dataclass
class CellAnnotations:
    """Annotations of cells in an imaging trial."""

    def __init__(self, cells_array: Union[List[int], pd.Index, pd.RangeIndex, np.ndarray], annotations: Union[List[str], pd.Index],
                 cellsdata: Union[pd.DataFrame, pd.Series], multidim_data: Dict[str, List[Any]] = None):
        """

        :param cells_array: ID of all cells in imaging dataset. must be of same cell length as imaging dataset.
        :param annotations: list of names of annotations
        :param cellsdata: M x Cells array of an arbritrary number (M) 1D annotations channels collected for all Cells. must contain same number of cells as cells_array.
        :param multidim_data: annotations with cellsdata of unconstrained dimensions for all cells. Structured as dictionary with keys corresponding to annotation name and a list of the length of cells containing cellsdata in any format.
        """

        self.cells_array = cells_array
        self.annotations = annotations
        self.cellsdata = cellsdata
        self.multidim_data = multidim_data

        assert len(self.cellsdata) == self.n_cells, 'mismatch of # of cells and cellsdata field.'
        print(f'\- added CellAnnotations module. consisting of {self.n_annotations} annotations x {self.n_cells} cells.')
        if self.multidim_data:
            for label, data in self.multidim_data.items():
                if not len(data) == self.n_cells:
                    raise ValueError(f"length of {label} of multidimdata does not match number of cells.")

    def __repr__(self):
        return f"CellAnnotations module, containing {self.n_cells} cells and {self.n_annotations} annotations."

    def __str__(self):
        return f"CellAnnotations module, containing {self.n_cells} cells and {self.n_annotations} annotations: \n{self.annotations}"

    # constructors:
    # @classmethod
    # def s2p_to_CellAnnotations(cls, s2pTrial):
    #     """alternative constructor for CellAnnotations from suite2p results."""
    #     cellsdata = s2pTrial.setCellsAnnotations()
    #     obj = cls(cells_array=cellsdata['original_index'], annotations=cellsdata.columns, cellsdata=cellsdata)
    #     return obj

    # properties:
    @property
    def n_cells(self):
        """number of cells
        """
        return len(self.cells_array)

    @property
    def n_annotations(self):
        """number of annotations
        """
        return len(self.annotations)

    @property
    def cell_id(self):
        """ID of cells. redundancy of .cells_array.
        """
        if 'cell_id' in self.cellsdata:
            return list(self.cellsdata['cell_id'])
        else:
            return list(self.cells_array)

    # functions:


@dataclass
class ImagingData:
    """Imaging dataset."""

    def __init__(self, imdata: Union[np.ndarray, pd.DataFrame], **kwargs):
        """
        Imaging dataset.

        :param imdata: N rois x num_frames, table of imaging data of cells (N) collected over time (Frames)
        """
        self.imdata = imdata  #: N rois x num_frames, table of imaging data of cells (N) collected over time (Frames)
        for i, values in kwargs.items():
            setattr(self, i, values)
        print(f'\- added ImagingData module. consisting of {self.imdata.shape} ROIs x Frames.')


    @property
    def n_frames(self):
        """number of imaging frames in imaging data
        """
        return self.imdata.shape[1]

    @property
    def n_rois(self):
        """number of ROIs in imaging data
        :return:
        """
        return self.imdata.shape[0]


class ImagingMetadata:
    """Metadata about imaging system parameters."""

    PIXEL_SIZE_UNITS: str = 'microns per pixel'  #: units for the values of the size of imaging pixels
    IMAGING_SIGNAL_UNITS: str = 'A.U.'  #: set as arbitrary units

    def __init__(self, microscope, n_frames, fps, frame_x, frame_y, n_planes, pix_sz_x, pix_sz_y, emission_lambda: float = 0., excitation_lambda: float = 0., **kwargs):
        self.microscope = microscope  #: given name of microscope
        self.n_frames = n_frames  #: number of imaging frames in the current trial
        self.fps = fps  #: rate of imaging acquisition (frames per second)
        self.frame_x = frame_x  #: num of pixels in the x direction of a single frame
        self.frame_y = frame_y  #: num of pixels in the y direction of a single frame
        self.n_planes = n_planes  #: num of FOV planes in imaging acquisition
        self.pix_sz_x = pix_sz_x  #: size of a single imaging pixel in x direction (microns per pixel)
        self.pix_sz_y = pix_sz_y  #: size of a single imaging pixel in y direction (microns per pixel)
        self.emission_lambda = emission_lambda  #: wavelength of optical channel imaging collection
        self.excitation_lambda = excitation_lambda  #: wavelength of excitation during imaging
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return f'ImagingMetadata for imaging cellsdata collected with {self.microscope}.'

    @property
    def fov_size(self):
        x_length_um = self.frame_x * self.pix_sz_x
        y_length_um = self.frame_y * self.pix_sz_y
        print(f'FOV size: (X): {round(x_length_um, 2)}microns, (Y): {round(y_length_um, 2)}microns')
        return x_length_um, y_length_um

