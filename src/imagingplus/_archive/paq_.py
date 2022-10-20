# TEMP copy of processing/paq.py to use for experimenting with and creating paq as child of temporal cellsdata.

import os.path
from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt

from packerlabimaging.main.subcore import TemporalData
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
            chan_name = chan_name + chr(int(np.fromfile(fid, dtype='>f', count=1)[0]))
        chan_names.append(chan_name)

    # get channel hardware lines
    hw_chans = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        hw_chan = ''
        for j in range(num_chars):
            hw_chan = hw_chan + chr(int(np.fromfile(fid, dtype='>f', count=1)[0]))
        hw_chans.append(hw_chan)

    # get acquisition units
    units = []
    for i in range(num_chans):
        num_chars = int(np.fromfile(fid, dtype='>f', count=1))
        unit = ''
        for j in range(num_chars):
            unit = unit + chr(int(np.fromfile(fid, dtype='>f', count=1)[0]))
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


# noinspection DuplicatedCode
class PaqData(TemporalData):
    """access and storage of cellsdata from .paq files."""

    def __init__(self, file_path, **kwargs):
        print(f"\n\- ADDING PAQ DATA from {file_path}... ")
        self.sparse_paq_data = None  #: array of paq cellsdata that corresponds to

        if 'channels' not in kwargs or 'sampling_rate' not in kwargs or 'paq_data' not in kwargs:
            paq_data, paq_rate, channels = self.paq_read(file_path=file_path, plot=True)
        else:
            paq_data = kwargs['paq_data']
            channels = kwargs['channels']
            paq_rate = kwargs['sampling_rate']

        # # todo switch this out in favour of the pandas dataframes below:
        # for chan_name in channels:
        #     chan_name_idx = channels.index(chan_name)
        #     print(f"\t- adding '{chan_name}' channel cellsdata as attribute")
        #     setattr(self, chan_name, paq_data['cellsdata'][chan_name_idx])

        data = pd.DataFrame(data=paq_data['cellsdata'].T, columns=channels, index=range(paq_data['cellsdata'].shape[1]))
        super(PaqData, self).__init__(file_path=file_path, channels=channels, sampling_rate=paq_rate, data=data)

        # init attr's
        self.stim_start_times = []  #: paq clock time of photostimulation trial start



    @classmethod
    def import_paqdata(cls, file_path, frame_times_channel=None, plot=False):
        """
        Alternative constructor for PaqData. This is the preferred method for loading in paqdata.

        :param file_path: path to .paq file
        :param frame_times_channel: channel to retrieve frame clock times from
        :param plot: whether to plot output of reading .paq file
        :return: PaqData object, and raw cellsdata from .paq file in numpy array

        todo add example in docstring
        """
        paq_data, paq_rate, channels = cls.paq_read(file_path=file_path, plot=plot)

        paqData_obj = cls(file_path=file_path, channels=channels, sampling_rate=paq_rate)

        if frame_times_channel:
            paqData_obj.getPaqFrameTimes(frame_times_channel=frame_times_channel)
            paqData_obj.sparse_paq_data = paqData_obj.get_sparse_data(frame_times=paqData_obj.frame_times)

        return paqData_obj

    def __str__(self):
        information = ""
        for i in self.__dict__:
            if type(self.__dict__[i]) != dict:
                information += f"\n\t{i}: {self.__dict__[i]}"
            else:
                information += f"\n\t{i}: {[*self.__dict__[i]]}"
        return f"imagingplus.processing.Paq.PaqData: {information}"

    @staticmethod
    def paq_read(file_path: str, plot: bool = False):
        """
        Loads .Paq file and saves cellsdata from individual channels.

        :param file_path: path to the .Paq file for this cellsdata object
        :param plot: (optional) whether to plot
        """
        assert os.path.exists(file_path), f'File path not found {file_path}'

        print(f'\tloading Paq cellsdata from: {file_path}')
        paq, _ = paq2py(file_path, plot=plot)
        paq_rate = paq['rate']
        channels = paq['chan_names']
        print(f"\t - loaded {len(paq['chan_names'])} channels from .Paq file: {paq['chan_names']}")

        return paq, paq_rate, channels

    def storePaqChannel(self, chan_name):
        """add a specific channel's (`chan_name`) cellsdata from the .paq file as attribute of the same name for
        PaqData object.

        :param chan_name: name of paq channel to add.
        """

        paq_data, _, channels = self.paq_read(file_path=self.file_path)
        chan_name_idx = channels.index(chan_name)
        print(f"\t|- adding '{chan_name}' channel cellsdata as attribute")
        setattr(self, chan_name, paq_data['cellsdata'][chan_name_idx])

    ## refactor these methods to their respective Trial code locations
    def getPaqFrameTimes(self, frame_times_channel: str):
        """
        Retrieve two-photon imaging frame times from .paq signal found in frame_times_channel.

        :param paq_data: cellsdata loaded from .paq file (use .paq_read method)
        :param frame_times_channel: channel to retrieve frame clock times from
        :return: numpy array of frame clock times
        """
        print(f"\n\t\- Retrieving two-photon imaging frame times from .paq channel: {frame_times_channel} ... ")

        # if frame_times_channel not in paq_data['chan_names']:
        if frame_times_channel not in self.channels:
            raise KeyError(f'{frame_times_channel} not found in .Paq channels. Specify channel containing frame signals.')

        # find frame times
        clock_voltage = self.data[frame_times_channel].to_numpy()

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

    def plot__paq_channel(self):
        """temp placeholder incase you need specific plotting code compared to plotting with the general temporal cellsdata function"""
        pass

    @classmethod
    def paqProcessingTwoPhotonImaging(cls, paq_path, frame_channel, plot: bool = False):
        """
        Alternative constructor for paq module for working with two photon imaging trials.

        :param plot:
        :param paq_path: path to .paq file
        :param frame_channel: channel to use for measuring frame times from .paq cellsdata

        :return: PAQ cellsdata object
        """

        paq_data_obj = cls.import_paqdata(file_path=paq_path, plot=plot)
        assert frame_channel in paq_data_obj.channels, f"frame_channel argument: '{frame_channel}', not found in channels in .paq cellsdata."
        paq_data_obj.frame_times = paq_data_obj.getPaqFrameTimes(frame_times_channel=frame_channel)
        paq_data_obj.sparse_paq_data = paq_data_obj.get_sparse_data(frame_times=paq_data_obj.frame_times)

        return paq_data_obj

    @classmethod
    def paqProcessingAllOptical(cls, paq_path: str, stim_channel: str, frame_channel: str, plot: bool =False):
        """
        Alternative constructor for paq module for working with all optical trials.

        :param plot:
        :param stim_channel:
        :param paq_path: path to .paq file
        :param frame_channel: channel to use for measuring frame times from .paq cellsdata

        :return: PAQ cellsdata object
        """

        paq_data_obj = cls.paqProcessingTwoPhotonImaging(paq_path=paq_path, frame_channel=frame_channel, plot=plot)

        assert stim_channel in paq_data_obj.channels, f"stim_channel argument: '{stim_channel}', not found in channels in .paq cellsdata."

        # find stim times
        stim_volts = paq_data_obj.data[stim_channel].to_numpy()
        stim_start_times = threshold_detect(stim_volts, 1)
        print('# of stims found on %s: %s' % (stim_channel, len(stim_start_times)))

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(stim_volts)
        ax.plot(stim_start_times, np.ones(len(stim_start_times)), '.')
        fig.suptitle(f'detected photostim start times from {stim_channel}')
        ax.set_xlabel(f'paq clock (sampled at {paq_data_obj.sampling_rate} Hz)')
        ax.set_ylabel(f"{stim_channel} volt")
        sns.despine()
        plt.show()

        stim_starts_ = []
        for time, v in enumerate(paq_data_obj.data[stim_channel]):
            if time in stim_start_times:
                stim_starts_.append(True)
            else:
                stim_starts_.append(False)

        paq_data_obj.data['stim_start_times'] = stim_starts_

        return paq_data_obj


    # # TODO review code
    # def _1p_stims(self, paq_data, plot: bool = False, optoloopback_channel: str = 'opto_loopback'):
    #     "find 1p stim times"
    #     if optoloopback_channel not in paq_data['chan_names']:
    #         raise KeyError(
    #             f'{optoloopback_channel} not found in .Paq channels. Specify channel containing 1p stim TTL loopback signals.')
    #
    #     opto_loopback_chan = paq_data['chan_names'].index('opto_loopback')
    #     stim_volts = paq_data['cellsdata'][opto_loopback_chan, :]
    #     stim_times = threshold_detect(stim_volts, 1)
    #
    #     self.stim_times = stim_times
    #     self.stim_start_times = [self.stim_times[0]]  # initialize ls
    #     self.stim_end_times = []
    #     i = len(self.stim_start_times)
    #     for stim in self.stim_times[1:]:
    #         if (stim - self.stim_start_times[i - 1]) > 1e5:
    #             i += 1
    #             self.stim_start_times.append(stim)
    #             self.stim_end_times.append(self.stim_times[np.where(self.stim_times == stim)[0] - 1][0])
    #     self.stim_end_times.append(self.stim_times[-1])
    #
    #     print("\nNumber of 1photon stims found: ", len(self.stim_start_times))
    #
    #     if plot:
    #         plt.figure(figsize=(50, 2))
    #         plt.plot(stim_volts)
    #         plt.plot(stim_times, np.ones(len(stim_times)), '.')
    #         plt.suptitle('1p stims from Paq, with detected 1p stim instances as scatter')
    #         plt.xlim([stim_times[0] - 2e3, stim_times[-1] + 2e3])
    #         plt.show()
    #
    #     # find all stim frames
    #     self.stim_frames = []
    #     for stim in range(len(self.stim_start_times)):
    #         stim_frames_ = [frame for frame, t in enumerate(self._frame_clock_actual) if
    #                         self.stim_start_times[stim] - 100 / self.paq_rate <= t <= self.stim_end_times[
    #                             stim] + 100 / self.paq_rate]
    #
    #         self.stim_frames.append(stim_frames_)
    #
    #     # if >1 1p stims per trial, find the start of all 1p trials
    #     self.stim_start_frames = [stim_frames[0] for stim_frames in self.stim_frames if len(stim_frames) > 0]
    #     self.stim_end_frames = [stim_frames[-1] for stim_frames in self.stim_frames if len(stim_frames) > 0]
    #     self.stim_duration_frames = int(np.mean(
    #         [self.stim_end_frames[idx] - self.stim_start_frames[idx] for idx in range(len(self.stim_start_frames))]))
    #
    #     print(f"\nStim duration of 1photon stim: {self.stim_duration_frames} frames")
    #
    # def _shutter_times(self, paq_data, shutter_channel: str = 'shutter_loopback'):
    #     """find shutter loopback frames from .Paq cellsdata
    #     :param paq_data:
    #     :param shutter_channel:
    #     """
    #
    #     if shutter_channel not in paq_data['chan_names']:
    #         raise KeyError(f'{shutter_channel} not found in .Paq channels. Specify channel containing shutter signals.')
    #
    #     shutter_idx = paq_data['chan_names'].index('shutter_loopback')
    #     shutter_voltage = paq_data['cellsdata'][shutter_idx, :]
    #
    #     shutter_times = np.where(shutter_voltage > 4)
    #     self.shutter_times = shutter_times[0]
    #     self.shutter_frames = []
    #     self.shutter_start_frames = []
    #     self.shutter_end_frames = []
    #
    #     shutter_frames_ = [frame for frame, t in enumerate(self.getPaqFrameTimes()) if
    #                        t in self.shutter_times]
    #     self.shutter_frames.append(shutter_frames_)
    #
    #     shutter_start_frames = [shutter_frames_[0]]
    #     shutter_end_frames = []
    #     i = len(shutter_start_frames)
    #     for frame in shutter_frames_[1:]:
    #         if (frame - shutter_start_frames[i - 1]) > 5:
    #             i += 1
    #             shutter_start_frames.append(frame)
    #             shutter_end_frames.append(shutter_frames_[shutter_frames_.index(frame) - 1])
    #     shutter_end_frames.append(shutter_frames_[-1])
    #     self.shutter_start_frames.append(shutter_start_frames)
    #     self.shutter_end_frames.append(shutter_end_frames)
    #
    # def frames_discard(self, input_array, total_frames, discard_all=False):
    #     '''
    #     calculate which 2P imaging frames to discard (or use as bad frames input into suite2p) based on the bad frames
    #     identified by manually inspecting the Paq files in EphysViewer.m
    #
    #     :param input_array: .m file path to read that contains the timevalues for signal to remove
    #     :param total_frames: the number of frames in the TIFF file of the actual 2p imaging recording
    #     :param discard_all: bool; if True, then add all 2p imaging frames from this Paq file as bad frames to discard
    #     :return: array that contains the indices of bad frames (in format ready to input into suite2p processing)
    #     '''
    #
    #     frame_times = self.getPaqFrameTimes()
    #     frame_times = frame_times[
    #                   0:total_frames]  # this is necessary as there are more TTL triggers in the Paq file than actual frames (which are all at the end)
    #
    #     all_btwn_paired_frames = []
    #     paired_frames_first = []
    #     paired_frames_last = []
    #     if input_array is not None:
    #         print('\nadding seizure frames loaded up from: ', input_array)
    #         measurements = io.loadmat(input_array)
    #         for set_ in range(len(measurements['PairedMeasures'])):
    #             # calculate the sample value for begin and end of the set
    #             begin = int(measurements['PairedMeasures'][set_][3][0][0] * self.paq_rate)
    #             end = int(measurements['PairedMeasures'][set_][5][0][0] * self.paq_rate)
    #             frames_ = list(np.where(np.logical_and(frame_times >= begin, frame_times <= end))[0])
    #             if len(frames_) > 0:
    #                 all_btwn_paired_frames.append(frames_)
    #                 paired_frames_first.append(frames_[0])
    #                 paired_frames_last.append(frames_[-1])
    #
    #         all_btwn_paired_frames = [item for x in all_btwn_paired_frames for item in x]
    #
    #     if discard_all and input_array is None:
    #         frames_to_discard = list(range(len(frame_times)))
    #         return frames_to_discard
    #     elif not discard_all and input_array is not None:
    #         frames_to_discard = all_btwn_paired_frames
    #         return frames_to_discard, all_btwn_paired_frames, paired_frames_first, paired_frames_last
    #     elif discard_all and input_array is not None:
    #         frames_to_discard = list(range(len(frame_times)))
    #         return frames_to_discard, all_btwn_paired_frames, paired_frames_first, paired_frames_last
    #     else:
    #         raise ReferenceError('something wrong....No frames selected for discarding')
    #
    # ## refactor these methods to their respective Trial code locations // END
