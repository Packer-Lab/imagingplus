import os
import numpy as np
import pandas as pd
import matplotlib.ticker as mticker

from .utils import threshold_detect

# paq2py by Llyod Russel
def paq_read(file_path=None, plot=False):
    """
    Read PAQ file (from PackIO) into python
    Lloyd Russell 2015
    Parameters
    ==========
    file_path : str, optional
        full path to file to read in. if none is supplied a load file dialog
        is opened, buggy on mac osx - Tk/matplotlib. Default: None.
    plot : bool, optional
        plot the data after reading? Default: False.
    Returns
    =======
    data : ndarray
        the data as a m-by-n array where m is the number of channels and n is
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

    # get data
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
        f, axes = plt.subplots(num_chans, 1, sharex=True, figsize=(15, num_chans*5), frameon=False)
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

    # make pandas data frame using data in channels
    df = pd.DataFrame(data.T, columns=chan_names)

    return {"data": data,
            "chan_names": chan_names,
            "hw_chans": hw_chans,
            "units": units,
            "rate": rate,
            "num_datapoints": num_datapoints}, df


from dataclasses import dataclass

@dataclass
class paqData:
    paq_path: str
    def __post_init__(self):
        self.paqProcessing(paq_path=self.paq_path, plot=False)

    def paqProcessing(self, paq_path: str = None, plot: bool = False):
        """
        Loads .paq file and saves data from individual channels.

        :param paq_path: (optional) path to the .paq file for this data object
        """

        print('\n\----- Processing paq file ...')

        if not hasattr(self, 'paq_path'):
            if paq_path is not None:
                self.paq_path = paq_path
            else:
                ValueError(
                    'ERROR: no paq_path defined for data object, please provide paq_path to load .paq file from.')
        elif paq_path is not None and paq_path != self.paq_path:
            assert os.path.exists(paq_path), print('ERROR: paq_path provided was not found')
            print(f"|- Updating paq_path to newly provided path: {paq_path}")
            self.paq_path = paq_path  # update paq_path if provided different path

        print(f'\tloading paq data from: {self.paq_path}')

        paq, _ = paq_read(self.paq_path, plot=plot)
        self.paq_rate = paq['rate']
        self.paq_channels = paq['chan_names']
        ## TODO print the paq channels that were loaded. and some useful metadata about the paq channels.
        print(f"\t|- loaded {len(paq['chan_names'])} channels from .paq file: {paq['chan_names']}")

        # find frame times
        clock_idx = paq['chan_names'].index('frame_clock')
        clock_voltage = paq['data'][clock_idx, :]

        __frame_clock = threshold_detect(clock_voltage, 1)
        self.__frame_clock = __frame_clock

        # find start and stop __frame_clock times -- there might be multiple 2p imaging starts/stops in the paq trial (hence multiple frame start and end times)
        self.frame_start_times = [self.__frame_clock[0]]  # initialize list
        self.frame_end_times = []
        i = len(self.frame_start_times)
        for idx in range(1, len(self.__frame_clock) - 1):
            if (self.__frame_clock[idx + 1] - self.__frame_clock[idx]) > 2e3:
                i += 1
                self.frame_end_times.append(self.__frame_clock[idx])
                self.frame_start_times.append(self.__frame_clock[idx + 1])
        self.frame_end_times.append(self.__frame_clock[-1])

        # handling cases where 2p imaging clock has been started/stopped >1 in the paq trial
        if len(self.frame_start_times) > 1:
            diff = [self.frame_end_times[idx] - self.frame_start_times[idx] for idx in
                    range(len(self.frame_start_times))]
            idx = diff.index(max(diff))
            self.frame_start_time_actual = self.frame_start_times[idx]
            self.frame_end_time_actual = self.frame_end_times[idx]
            self.__frame_clock_actual = [frame for frame in self.__frame_clock if
                                         self.frame_start_time_actual <= frame <= self.frame_end_time_actual]
        else:
            self.frame_start_time_actual = self.frame_start_times[0]
            self.frame_end_time_actual = self.frame_end_times[0]
            self.__frame_clock_actual = self.__frame_clock

        # read in and save sparse version of all paq channels (only save data from timepoints at frame clock times)
        self.sparse_paq_data = {}
        for idx, chan in enumerate(self.paq_channels):
            self.sparse_paq_data[chan] = paq['data'][idx, self.__frame_clock_actual]

