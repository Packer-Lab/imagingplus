## - archived on may 21 2022

import os.path
from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from scipy import io

from packerlabimaging.utils.utils import threshold_detect


def paq2py(file_path=None, plot=False):
    """
    Read PAQ file (from PackIO) into python
    Lloyd Russell 2015. Minor update for numpy >1.18 by Prajay Shah 2021.
    Parameters
    ==========
    file_path : str, optional
        full path to file to read in. if none is supplied a load file dialog
        is opened, buggy on mac osx - Tk/matplotlib. Default: None.
    plot : bool, optional
        plot the cellsdata after reading? Default: False.
    Returns
    =======
    cellsdata : ndarray
        the cellsdata as a m-by-n array where m is the number of channels and n is
        the number of datapoints
    chan_names : list of str
        the names of the channels provided in PackIO
    hw_chans : list of str
        the hardware lines corresponding to each channel
    units : list of str
        the units of measurement for each channel
    rate : int
        the acquisition sample rate, in Hz
    """

    # file load gui
    if file_path is None:
        import Tkinter
        import tkFileDialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename()
        root.destroy()

    # open file
    fid = open(file_path, 'rb')

    # get sample rate
    rate = int(np.fromfile(fid, dtype='>f', count=1))

    # get number of channels
    num_chans = int(np.fromfile(fid, dtype='>f', count=1))

    # get channel names
    chan_names = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        chan_name = ''
        for j in range(num_chars):
            chan_name = chan_name + chr(np.fromfile(fid, dtype='>f', count=1)[0])
        chan_names.append(chan_name)

    # get channel hardware lines
    hw_chans = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        hw_chan = ''
        for j in range(num_chars):
            hw_chan = hw_chan + chr(np.fromfile(fid, dtype='>f', count=1)[0])
        hw_chans.append(hw_chan)

    # get acquisition units
    units = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        unit = ''
        for j in range(num_chars):
            unit = unit + chr(np.fromfile(fid, dtype='>f', count=1)[0])
        units.append(unit)

    # get cellsdata
    temp_data = np.fromfile(fid, dtype='>f', count=-1)
    num_datapoints = int(len(temp_data) / num_chans)
    data = np.reshape(temp_data, [num_datapoints, num_chans]).transpose()

    # close file
    fid.close()

    # plot
    if plot:
        # import matplotlib
        # matplotlib.use('QT4Agg')
        import matplotlib.pylab as plt
        f, axes = plt.subplots(num_chans, 1, sharex=True, figsize=(15, num_chans * 5), frameon=False)
        for idx, ax in enumerate(axes):
            ax.plot(data[idx])
            ax.set_xlim([0, num_datapoints - 1])
            ax.set_ylim([data[idx].min() - 1, data[idx].max() + 1])
            # ax.set_ylabel(units[idx])
            ax.set_title(chan_names[idx])

            # -- Prajay edit
            # change x axis ticks to seconds
            label_format = '{:,.0f}'
            labels = [item for item in ax.get_xticks()]
            for item in labels:
                labels[labels.index(item)] = int(round(item / rate))
            ticks_loc = ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax.set_xticklabels([label_format.format(x) for x in labels])
            ax.set_xlabel('Time (secs)')
            # --

        plt.suptitle(file_path)
        plt.tight_layout()
        plt.show()

    # make pandas cellsdata frame using cellsdata in channels
    df = pd.DataFrame(data.T, columns=chan_names)

    return {"cellsdata": data,
            "chan_names": chan_names,
            "hw_chans": hw_chans,
            "units": units,
            "rate": rate,
            "num_datapoints": num_datapoints}, df



@dataclass
class PaqData:
    """access and storage of cellsdata from .paq files."""

    paq_path: str
    paq_channels: List[str] = None
    paq_rate: float = None
    def __post_init__(self):
        print(f"\n\- ADDING PAQ MODULE from {self.paq_path}... ")

    #     _, paq_rate, paq_channels = self.paq_read(paq_path=self.paq_path, plot=False)
    #     self.paq_channels = paq_channels
    #     self.paq_rate = paq_rate

    @classmethod
    def import_paqdata(cls, paq_path, plot=False):
        """
        Alternative constructor for PaqData. This is the preferred method for loading in paqdata.

        :param paq_path: path to .paq file
        :param plot: whether to plot output of reading .paq file
        :return: PaqData object, and raw cellsdata from .paq file in numpy array

        todo add example in docstring
        """
        paqData_obj = cls(paq_path=paq_path)

        paq_data, paq_rate, paq_channels = paqData_obj.paq_read(paq_path=paq_path, plot=plot)
        paqData_obj.paq_channels = paq_channels
        paqData_obj.paq_rate = paq_rate

        for chan_name in paqData_obj.paq_channels:
            chan_name_idx = paq_channels.index(chan_name)
            print(f"\t- adding '{chan_name}' channel cellsdata as attribute")
            setattr(paqData_obj, chan_name, paq_data['cellsdata'][chan_name_idx])

        return paqData_obj, paq_data

    def __repr__(self):
        information = ""
        for i in self.__dict__:
            if type(self.__dict__[i]) != dict:
                information += f"\n\t{i}: {self.__dict__[i]}"
            else:
                information += f"\n\t{i}: {[*self.__dict__[i]]}"
        return f"imagingplus.processing.Paq.PaqData: {information}"

    @staticmethod
    def paq_read(paq_path: str, plot: bool = False):
        """
        Loads .Paq file and saves cellsdata from individual channels.

        :param paq_path: path to the .Paq file for this cellsdata object
        :param plot: (optional) whether to plot
        """
        assert os.path.exists(paq_path), f'File path not found {paq_path}'

        print(f'\tloading Paq cellsdata from: {paq_path}')
        paq, _ = paq2py(paq_path, plot=plot)
        paq_rate = paq['rate']
        paq_channels = paq['chan_names']
        print(f"\t - loaded {len(paq['chan_names'])} channels from .Paq file: {paq['chan_names']}")

        return paq, paq_rate, paq_channels

    def storePaqChannel(self, chan_name):
        """add a specific channel's (`chan_name`) cellsdata from the .paq file as attribute of the same name for
        PaqData object.

        :param chan_name: name of paq channel to add.
        """

        paq_data, _, paq_channels = self.paq_read(paq_path=self.paq_path)
        chan_name_idx = paq_channels.index(chan_name)
        print(f"\t|- adding '{chan_name}' channel cellsdata as attribute")
        setattr(self, chan_name, paq_data['cellsdata'][chan_name_idx])

    def cropPaqData(self, begin: int, end: int, channels: List[str] = 'all'):
        """
        Crops saved paq cellsdata channels to the .paq clock timestamps of begin and end.

        :param begin: paq clock time to begin cropping at
        :param end: paq clock time to end cropping at
        :param channels: channels to crop paq cellsdata.
        """

        channels = self.paq_channels if channels == 'all' else channels
        for channel in channels:
            print(f"\- cropping {channel} to {begin} and {end} paq clock times.")
            data = getattr(self, channel)
            assert len(data) >= (end - begin), f'{channel} paq cellsdata is not long enough to crop between the provided clock times.'
            cropdata = data[begin: end]
            setattr(self, channel, cropdata)


    ## refactor these methods to their respective Trial code locations
    def paq_frame_times(self, frame_channel: str):
        """
        Retrieve two-photon imaging frame times from .paq signal found in frame_channel.

        :param paq_data: cellsdata loaded from .paq file (use .paq_read method)
        :param frame_channel: channel to retrieve frame clock times from
        :return: numpy array of frame clock times
        """
        print(f"\n\t\- Retrieving two-photon imaging frame times from .paq channel: {frame_channel} ... ")

        # if frame_channel not in paq_data['chan_names']:
        if frame_channel not in self.paq_channels:
            raise KeyError(f'{frame_channel} not found in .Paq channels. Specify channel containing frame signals.')

        # find frame times
        # clock_idx = paq_data['chan_names'].index(frame_channel)
        # clock_voltage = paq_data['cellsdata'][clock_idx, :]
        clock_voltage = getattr(self, frame_channel)

        __frame_clock = threshold_detect(clock_voltage, 1)
        __frame_clock = __frame_clock

        # find start and stop __frame_clock times -- there might be multiple 2p imaging starts/stops in the Paq trial (hence multiple frame start and end times)
        frame_start_times = [__frame_clock[0]]  # initialize list
        frame_end_times = []
        i = len(frame_start_times)
        for idx in range(1, len(__frame_clock) - 1):
            if (__frame_clock[idx + 1] - __frame_clock[idx]) > 2e3:
                i += 1
                frame_end_times.append(__frame_clock[idx])
                frame_start_times.append(__frame_clock[idx + 1])
        frame_end_times.append(__frame_clock[-1])

        # handling cases where 2p imaging clock has been started/stopped >1 in the Paq trial
        if len(frame_start_times) > 1:
            diff = [frame_end_times[idx] - frame_start_times[idx] for idx in
                    range(len(frame_start_times))]
            idx = diff.index(max(diff))
            frame_start_time_actual = frame_start_times[idx]
            frame_end_time_actual = frame_end_times[idx]
            frame_clock_actual = [frame for frame in __frame_clock if
                                  frame_start_time_actual <= frame <= frame_end_time_actual]
        else:
            frame_start_time_actual = frame_start_times[0]
            frame_end_time_actual = frame_end_times[0]
            frame_clock_actual = __frame_clock

        return frame_clock_actual

    def get_sparse_paq(self, frame_clock):
        """
        Returns dictionary of numpy array keyed on channels from paq_data timed to 2photon imaging frame_times.

        :param paq_data:
        :param frame_clock:
        :return:
        """
        print(f"\n\t\- Getting imaging frames timed .paq cellsdata from {len(frame_clock)} frames ... ")

        # read in and save sparse version of all Paq channels (only save cellsdata from timepoints at frame clock times)
        sparse_paq_data = {}
        for idx, chan in enumerate(self.paq_channels):
            data = getattr(self, chan)
            sparse_paq_data[chan] = data[frame_clock]
        return sparse_paq_data

    @staticmethod
    def paq_alloptical_stims(paq_data, frame_clock: List[int], stim_channel: str, plot: bool = False):
        if stim_channel not in paq_data['chan_names']:
            raise KeyError(f'{stim_channel} not found in .Paq channels. Specify channel containing frame signals.')

        # find stim times
        stim_idx = paq_data['chan_names'].index(stim_channel)
        stim_volts = paq_data['cellsdata'][stim_idx, :]
        stim_times = threshold_detect(stim_volts, 1)
        stim_start_times = stim_times
        print(f'# of stims found on {stim_channel}: {len(stim_start_times)}')

        if plot:
            plt.figure(figsize=(10, 5))
            plt.plot(stim_volts)
            plt.plot(stim_times, np.ones(len(stim_times)), '.')
            plt.suptitle('stim times')
            sns.despine()
            plt.show()

        # TODO need to figure out how to handle multi-plane imaging and retrieving stim frame times
        # find stim frames
        stim_start_frames = []
        for stim in stim_times:
            # the index of the frame immediately preceeding stim
            stim_start_frame = next(
                i - 1 for i, sample in enumerate(frame_clock) if sample - stim >= 0)
            stim_start_frames.append(stim_start_frame)

        return stim_start_frames, stim_start_times

    # TODO review code
    def _1p_stims(self, paq_data, plot: bool = False, optoloopback_channel: str = 'opto_loopback'):
        "find 1p stim times"
        if optoloopback_channel not in paq_data['chan_names']:
            raise KeyError(
                f'{optoloopback_channel} not found in .Paq channels. Specify channel containing 1p stim TTL loopback signals.')

        opto_loopback_chan = paq_data['chan_names'].index('opto_loopback')
        stim_volts = paq_data['cellsdata'][opto_loopback_chan, :]
        stim_times = threshold_detect(stim_volts, 1)

        self.stim_times = stim_times
        self.stim_start_times = [self.stim_times[0]]  # initialize ls
        self.stim_end_times = []
        i = len(self.stim_start_times)
        for stim in self.stim_times[1:]:
            if (stim - self.stim_start_times[i - 1]) > 1e5:
                i += 1
                self.stim_start_times.append(stim)
                self.stim_end_times.append(self.stim_times[np.where(self.stim_times == stim)[0] - 1][0])
        self.stim_end_times.append(self.stim_times[-1])

        print("\nNumber of 1photon stims found: ", len(self.stim_start_times))

        if plot:
            plt.figure(figsize=(50, 2))
            plt.plot(stim_volts)
            plt.plot(stim_times, np.ones(len(stim_times)), '.')
            plt.suptitle('1p stims from Paq, with detected 1p stim instances as scatter')
            plt.xlim([stim_times[0] - 2e3, stim_times[-1] + 2e3])
            plt.show()

        # find all stim frames
        self.stim_frames = []
        for stim in range(len(self.stim_start_times)):
            stim_frames_ = [frame for frame, t in enumerate(self._frame_clock_actual) if
                            self.stim_start_times[stim] - 100 / self.paq_rate <= t <= self.stim_end_times[
                                stim] + 100 / self.paq_rate]

            self.stim_frames.append(stim_frames_)

        # if >1 1p stims per trial, find the start of all 1p trials
        self.stim_start_frames = [stim_frames[0] for stim_frames in self.stim_frames if len(stim_frames) > 0]
        self.stim_end_frames = [stim_frames[-1] for stim_frames in self.stim_frames if len(stim_frames) > 0]
        self.stim_duration_frames = int(np.mean(
            [self.stim_end_frames[idx] - self.stim_start_frames[idx] for idx in range(len(self.stim_start_frames))]))

        print(f"\nStim duration of 1photon stim: {self.stim_duration_frames} frames")

    def _shutter_times(self, paq_data, shutter_channel: str = 'shutter_loopback'):
        """find shutter loopback frames from .Paq cellsdata
        :param paq_data:
        :param shutter_channel:
        """

        if shutter_channel not in paq_data['chan_names']:
            raise KeyError(f'{shutter_channel} not found in .Paq channels. Specify channel containing shutter signals.')

        shutter_idx = paq_data['chan_names'].index('shutter_loopback')
        shutter_voltage = paq_data['cellsdata'][shutter_idx, :]

        shutter_times = np.where(shutter_voltage > 4)
        self.shutter_times = shutter_times[0]
        self.shutter_frames = []
        self.shutter_start_frames = []
        self.shutter_end_frames = []

        shutter_frames_ = [frame for frame, t in enumerate(self.paq_frame_times()) if
                           t in self.shutter_times]
        self.shutter_frames.append(shutter_frames_)

        shutter_start_frames = [shutter_frames_[0]]
        shutter_end_frames = []
        i = len(shutter_start_frames)
        for frame in shutter_frames_[1:]:
            if (frame - shutter_start_frames[i - 1]) > 5:
                i += 1
                shutter_start_frames.append(frame)
                shutter_end_frames.append(shutter_frames_[shutter_frames_.index(frame) - 1])
        shutter_end_frames.append(shutter_frames_[-1])
        self.shutter_start_frames.append(shutter_start_frames)
        self.shutter_end_frames.append(shutter_end_frames)

    def frames_discard(self, input_array, total_frames, discard_all=False):
        '''
        calculate which 2P imaging frames to discard (or use as bad frames input into suite2p) based on the bad frames
        identified by manually inspecting the Paq files in EphysViewer.m

        :param input_array: .m file path to read that contains the timevalues for signal to remove
        :param total_frames: the number of frames in the TIFF file of the actual 2p imaging recording
        :param discard_all: bool; if True, then add all 2p imaging frames from this Paq file as bad frames to discard
        :return: array that contains the indices of bad frames (in format ready to input into suite2p processing)
        '''

        frame_times = self.paq_frame_times()
        frame_times = frame_times[
                      0:total_frames]  # this is necessary as there are more TTL triggers in the Paq file than actual frames (which are all at the end)

        all_btwn_paired_frames = []
        paired_frames_first = []
        paired_frames_last = []
        if input_array is not None:
            print('\nadding seizure frames loaded up from: ', input_array)
            measurements = io.loadmat(input_array)
            for set_ in range(len(measurements['PairedMeasures'])):
                # calculate the sample value for begin and end of the set
                begin = int(measurements['PairedMeasures'][set_][3][0][0] * self.paq_rate)
                end = int(measurements['PairedMeasures'][set_][5][0][0] * self.paq_rate)
                frames_ = list(np.where(np.logical_and(frame_times >= begin, frame_times <= end))[0])
                if len(frames_) > 0:
                    all_btwn_paired_frames.append(frames_)
                    paired_frames_first.append(frames_[0])
                    paired_frames_last.append(frames_[-1])

            all_btwn_paired_frames = [item for x in all_btwn_paired_frames for item in x]

        if discard_all and input_array is None:
            frames_to_discard = list(range(len(frame_times)))
            return frames_to_discard
        elif not discard_all and input_array is not None:
            frames_to_discard = all_btwn_paired_frames
            return frames_to_discard, all_btwn_paired_frames, paired_frames_first, paired_frames_last
        elif discard_all and input_array is not None:
            frames_to_discard = list(range(len(frame_times)))
            return frames_to_discard, all_btwn_paired_frames, paired_frames_first, paired_frames_last
        else:
            raise ReferenceError('something wrong....No frames selected for discarding')

    ## refactor these methods to their respective Trial code locations // END
