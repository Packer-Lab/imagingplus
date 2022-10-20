# add electrophysiology cellsdata to the trial

# retrieving and processing on LFP recordings from the .Paq file
from typing import Union

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

from imagingplus.main.core import ImagingTrial
from imagingplus.plotting._utils import plotting_decorator


class LFP:
    def __init__(self, data, rate, **kwargs):
        self.data = data  #: 1-D LFP cellsdata
        self.rate = rate  #: rate of sampling

    @classmethod
    def lfp_from_trial(cls, trialobj: ImagingTrial, channel: str = 'voltage'):
        """Alternative constructor for LFP. """
        print('\n----- Creating new LFP analysis submodule from TemporalData ...')

        if channel not in trialobj.tmdata.channels:
            raise ValueError('channel not found in provided tmdata.sparse_data')

        return cls(data=trialobj.tmdata.sparse_data[channel], rate=trialobj.imparams.fps)

    def detrend(self, data: Union[list, np.ndarray] = None):
        data = self.data if data is None else data
        return signal.detrend(data)

    @plotting_decorator(figsize=(4, 4))
    def run_spectrogram(self, V, t, fs, start=None, stop=None, **kwargs):
        ax = kwargs['ax']
        fig = kwargs['fig']
        # create a power spectrum of selected cellsdata
        if start or stop:
            t_sub = np.where((t > start) & (t < stop))
            V = V[0][t_sub]

        powerSpectrum, freqenciesFound, time, imageAxis = ax.specgram(V, Fs=fs,
                                                                      vmin=-60)
        ax.set_xlabel('Time (secs)')
        ax.set_ylabel('Frequency')
        ax.set_ylim([0, 200])
        if 'title' in kwargs.keys():
            ax.set_title(kwargs['title'])
        if 'colorbar' in kwargs.keys():
            if kwargs['colorbar']:
                fig.colorbar(imageAxis, ax=ax)

    def run_correlogram(self, signalA: Union[list, np.ndarray] = None, signalB: Union[list, np.ndarray] = None):
        signalA = self.detrend() if signalA is None else signalA
        # use scipy.signal.correlate to correlate the two timeseries
        correlation_strength = signal.correlate(signalA, signalB)
        lags = signal.correlation_lags(len(signalA), len(signalB))
        correlation_strength /= np.max(correlation_strength)

    def power_spectral_density(self, data: Union[list, np.ndarray] = None):

        data = self.detrend() if data is None else data
        fs = self.rate

        fourier_transform = np.fft.rfft(data)

        abs_fourier_transform = np.abs(fourier_transform)

        power_spectrum = np.square(abs_fourier_transform)


        frequency = np.linspace(0, fs / 2, len(power_spectrum))
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot(frequency, power_spectrum)
        ax.set_xlim([0.1, fs / 2])
        ax.set_yscale('log')
        fig.show()

