# use this in a tutorial as example of building new workflow from existing TwoPhotonImaging

import os
import time

import numpy as np
import matplotlib.pyplot as plt

# grabbing functions from .utils_funcs that are used in this script - Prajay's edits (review based on need)
from imagingplus._archive.paq import PaqData, paq2py

from imagingplus._archive.TwoPhotonImagingMain import TwoPhotonImagingTrial


class OnePhotonStim(TwoPhotonImagingTrial):
    def __init__(self, trial, metainfo, microscope: str, analysis_save_path: str = None):
        print('----------------------------------------')
        print('-----Processing trial # %s-----' % trial)
        print('----------------------------------------\n')

        super().__init__(self.tiff_path, self._paq_path, metainfo=metainfo, analysis_save_path=analysis_save_path, microscope=microscope, **kwargs)

        self.paq = PaqData(paq_path=self._paq_path)
        self.paq._1p_stim()

        self.save_pkl(pkl_path=self.pkl_path)

        print('\n-----DONE OnePhotonStim init of trial # %s-----' % trial)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if self.pkl_path:
            lastmod = time.ctime(os.path.getmtime(self.pkl_path))
            information = self.t_series_name
            return f"({information}) TwoPhotonImagingTrial.OnePhotonStim experimental trial object, last saved: {lastmod}"
        else:
            return f" -- unsaved TwoPhotonImagingTrial.OnePhotonStim experimental trial object -- "

    def _paqProcessingTwoPhotonImaging(self, **kwargs):

        print('\n-----processing Paq file for 1p photostim...')

        print('loading', self._paq_path)

        paq, _ = paq2py(self._paq_path, plot=True)
        self.paq_rate = paq['rate']

        frame_rate = self.fps / self.n_planes

        # if 'shutter_loopback' in Paq['chan_names']:
        #     ans = input('shutter_loopback in this Paq found, should we continue')
        #     if ans is True or 'Yes':
        #         pass
        #     else:
        #         raise Exception('need to write code for using the shutter loopback')

        # find frame_times times
        clock_idx = paq['chan_names'].index('frame_times')
        clock_voltage = paq['cellsdata'][clock_idx, :]

        frame_clock = pj.threshold_detect(clock_voltage, 1)
        self.frame_times = frame_clock

        # find start and stop frame_times times -- there might be multiple 2p imaging starts/stops in the Paq trial (hence multiple frame start and end times)
        self.frame_start_times = [self.frame_times[0]]  # initialize ls
        self.frame_end_times = []
        i = len(self.frame_start_times)
        for idx in range(1, len(self.frame_times) - 1):
            if (self.frame_times[idx + 1] - self.frame_times[idx]) > 2e3:
                i += 1
                self.frame_end_times.append(self.frame_times[idx])
                self.frame_start_times.append(self.frame_times[idx + 1])
        self.frame_end_times.append(self.frame_times[-1])

        # for frame in self.frame_times[1:]:
        #     if (frame - self.frame_start_times[i - 1]) > 2e3:
        #         i += 1
        #         self.frame_start_times.append(frame)
        #         self.frame_end_times.append(self.frame_times[np.where(self.frame_times == frame)[0] - 1][0])
        # self.frame_end_times.append(self.frame_times[-1])

        # handling cases where 2p imaging clock has been started/stopped >1 in the Paq trial
        if len(self.frame_start_times) > 1:
            diff = [self.frame_end_times[idx] - self.frame_start_times[idx] for idx in
                    range(len(self.frame_start_times))]
            idx = diff.index(max(diff))
            self.frame_start_time_actual = self.frame_start_times[idx]
            self.frame_end_time_actual = self.frame_end_times[idx]
            self.frame_times_actual = [frame for frame in self.frame_times if
                                       self.frame_start_time_actual <= frame <= self.frame_end_time_actual]
        else:
            self.frame_start_time_actual = self.frame_start_times[0]
            self.frame_end_time_actual = self.frame_end_times[0]
            self.frame_times_actual = self.frame_times

        f, ax = plt.subplots(figsize=(20, 2))
        # plt.figure(figsize=(50, 2))
        ax.plot(clock_voltage)
        ax.plot(frame_clock, np.ones(len(frame_clock)), '.', color='orange')
        ax.plot(self.frame_times_actual, np.ones(len(self.frame_times_actual)), '.', color='red')
        ax.set_title('frame clock from Paq, with detected frame clock instances as scatter')
        ax.set_xlim([1e6, 1.2e6])
        f.tight_layout(pad=2)
        f.show()



        # i = len(self.stim_start_frames)
        # for stim in self.stim_frames[1:]:
        #     if (stim - self.stim_start_frames[i-1]) > 100:
        #         i += 1
        #         self.stim_start_frames.append(stim)

        # # sanity check
        # assert max(self.stim_start_frames[0]) < self.raw[plane].shape[1] * self.n_planes

        # find voltage channel and save as lfp_signal attribute
        voltage_idx = paq['chan_names'].index('voltage')
        self.lfp_signal = paq['cellsdata'][voltage_idx]

    def collect_seizures_info(self, seizures_lfp_timing_matarray=None, discard_all=True):
        print('\ncollecting information about seizures...')
        self.seizures_lfp_timing_matarray = seizures_lfp_timing_matarray  # path to the matlab array containing paired measurements of seizures onset and offsets

        # retrieve seizure onset and offset times from the seizures info array input
        paq = paq2py(file_path=self._paq_path, plot=False)

        # print(Paq[0]['cellsdata'][0])  # print the frame clock signal from the .Paq file to make sure its being read properly
        # NOTE: the output of all of the following function is in dimensions of the FRAME CLOCK (not official Paq clock time)
        if seizures_lfp_timing_matarray is not None:
            print('-- using matlab array to collect seizures %s: ' % seizures_lfp_timing_matarray)
            bad_frames, self.seizure_frames, self.seizure_lfp_onsets, self.seizure_lfp_offsets = frames_discard(
                paq=paq[0], input_array=seizures_lfp_timing_matarray, total_frames=self.n_frames,
                discard_all=discard_all)
        else:
            print('-- no matlab array given to use for finding seizures.')
            self.seizure_frames = []
            bad_frames = frames_discard(paq=paq[0], input_array=seizures_lfp_timing_matarray,
                                        total_frames=self.n_frames,
                                        discard_all=discard_all)

        print('\nTotal extra seizure/CSD or other frames to discard: ', len(bad_frames))
        print('|- first and last 10 indexes of these frames', bad_frames[:10], bad_frames[-10:])

        if seizures_lfp_timing_matarray is not None:
            # print('|-now creating raw movies for each sz as well (saved to the /Analysis folder) ... ')
            # self.subselect_tiffs_sz(onsets=self.seizure_lfp_onsets, offsets=self.seizure_lfp_offsets,
            #                         on_off_type='lfp_onsets_offsets')

            print('|-now classifying photostims at phases of seizures ... ')
            self.stims_in_sz = [stim for stim in self.stim_start_frames if stim in self.seizure_frames]
            self.stims_out_sz = [stim for stim in self.stim_start_frames if stim not in self.seizure_frames]
            self.stims_bf_sz = [stim for stim in self.stim_start_frames
                                for sz_start in self.seizure_lfp_onsets
                                if -2 * self.fps < (
                                        sz_start - stim) < 2 * self.fps]  # select stims that occur within 2 seconds before of the sz onset
            self.stims_af_sz = [stim for stim in self.stim_start_frames
                                for sz_start in self.seizure_lfp_offsets
                                if -2 * self.fps < -1 * (
                                        sz_start - stim) < 2 * self.fps]  # select stims that occur within 2 seconds afterof the sz offset
            print(' \n|- stims_in_sz:', self.stims_in_sz, ' \n|- stims_out_sz:', self.stims_out_sz,
                  ' \n|- stims_bf_sz:', self.stims_bf_sz, ' \n|- stims_af_sz:', self.stims_af_sz)

        else:
            print('|- No matlab measurement array given so setting all stims as outside of sz ... ')
            self.stims_in_sz = []
            self.stims_out_sz = [stim for stim in self.stim_start_frames if stim not in self.seizure_frames]
            self.stims_bf_sz = []
            self.stims_af_sz = []

        aoplot.plot_lfp_stims(self, x_axis='time')
        self.save_pkl()
