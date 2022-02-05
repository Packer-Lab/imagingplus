# STA Movie Maker
# Lloyd Russell 2017
# No GUI edit Rob Lees 2020

import sys; sys.path.append('/home/pshah/Documents/code/Vape/')
sys.path.append('/home/pshah/Documents/code/PackerLab_pycharm/')
from _utils_ import paq_utils
import os
import ctypes
import json
import numpy as np
import difflib
import glob
import colorsys
import scipy
from scipy import io as sio
# from skimage.external import tifffile
import tifffile
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter
import time
import tensorflow as tf

from Vape.utils import sta, PrairieLink, paq2py, ThorLink
from skimage import exposure


class STAMovieMaker():

    def __init__(self, arg_dict):
        self.p = arg_dict
        self.work()

    def work(self):
        error = False

        # make save directory (/STA)
        base_path, file_name = os.path.split(self.p['moviePath'])
        file_name_no_ext = os.path.splitext(file_name)[0]
        save_dir = self.p['savePath']
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # get params
        sync_frame_channel = self.p['syncFrameChannel']
        sync_stim_channel = self.p['syncStimChannel']
        sync_start = self.p['syncStartSec']
        sync_stop = self.p['syncStopSec']
        num_diff_stims = self.p['numDiffStims']
        start_on_stim = self.p['startOnStim'] - 1
        only_every_x_stims = self.p['everyXStims']
        pre_sec = self.p['preSeconds']
        post_sec = self.p['postSeconds']
        frame_rate = self.p['frameRate']
        avg_image_start_sec = pre_sec + self.p['averageImageStart']
        avg_image_stop_sec = pre_sec + self.p['averageImageStop']
        avg_image_frames = range(int(avg_image_start_sec*frame_rate), int(avg_image_stop_sec*frame_rate))

        if sync_stop == 0:
            sync_stop = np.Inf

        # decide which normalisation methods
        methods = []
        if self.p['methodDF']:
            methods.append('dF')
        if self.p['methodDFF']:
            methods.append('dFF')
        if self.p['methodZscore']:
            methods.append('Zscore')

        # load sync file
        print('Loading sync file')
        sync_ext = os.path.splitext(self.p['syncPath'])[1]
        if sync_ext == '.paq':
            try:
                paq, _ = paq_utils.paq_read(self.p['syncPath'], plot=False)
                frame_trace = paq['data'][paq['chan_names'].index(sync_frame_channel)]
                stim_trace = paq['data'][paq['chan_names'].index(sync_stim_channel)]
                rate = paq['rate']

            except:
                print('Error. Channel names: ' + str(paq['chan_names']))
                error = True

        elif sync_ext == '.h5':
            sync_data = ThorLink.ReadSyncFile(self.p['syncPath'],
                    datasets=[sync_frame_channel, sync_stim_channel], listdatasets=True)
            frame_trace = sync_data[sync_frame_channel]
            stim_trace = sync_data[sync_stim_channel]
            rate = 250000
        
        elif sync_ext == '.txt':
            # comma separated ls of stim frames (rows can be diff stim types)
            pass

        if not error:
            # load movie
            print('Loading movie')
            movie_ext = os.path.splitext(self.p['moviePath'])[1]
            if movie_ext == '.bin':
                movie = PrairieLink.ReadRawFile(self.p['moviePath'])
            elif movie_ext == '.raw':
                movie = ThorLink.ReadRawFile(self.p['moviePath'])
            elif movie_ext == '.tif':
                movie = tifffile.TiffFile(self.p['moviePath'], multifile=True).asarray()

            # get movie dimensions
            num_frames = movie.shape[0]
            frame_dims = movie.shape[1:]

            # get frame times
            frame_times = sta.threshold_detect(
            frame_trace, 1).astype(np.float) / rate
            frame_times = frame_times[(frame_times > sync_start) & (
            frame_times < sync_stop)]
            frame_times = frame_times[0:num_frames]
            frame_times_all = frame_times
            
            # get stim times
            all_stim_times = sta.threshold_detect(stim_trace, 1).astype(np.float) / rate
            all_stim_times = all_stim_times[(all_stim_times > sync_start) & (all_stim_times < sync_stop)]
            all_stim_times = all_stim_times[(all_stim_times-pre_sec > min(frame_times)) & 
                                            (all_stim_times+post_sec < max(frame_times))]

            # NUMBER OF PLANES HERE
            numZPlanes = self.p['zPlanes']

            for z in range(numZPlanes):
                frame_times = frame_times_all[z::numZPlanes]
                frame_indices = np.arange(z, num_frames, numZPlanes)

                est_frame_rate = np.round(1/np.diff(frame_times).mean())
                if not est_frame_rate == frame_rate:
                    print('Error. Estimated frame rate: ' + str(
                                            est_frame_rate) + ', Specified frame rate: ' + str(frame_rate))
                    time.sleep(2)

                # make STA template
                pre_samples = round(pre_sec * frame_rate)
                post_samples = round(post_sec * frame_rate)
                sta_template = np.arange(-pre_samples*numZPlanes, post_samples*numZPlanes, numZPlanes)

                # search for stim order file
                if self.p['useStimOrder']:
                    stim_list = sio.loadmat(self.p['stimOrder'])
                    # this should change!
                    stim_order = list(stim_list['oris'])
                    stim_order = stim_order[:len(all_stim_times)]
                    unique_stims = np.unique(stim_order)
                    print('Using stim order file')

                # store all average images in dict to combine later
                avg_img = {}
                if methods:
                    for method in methods:
                        avg_img[method] = np.ndarray(
                                            [num_diff_stims, frame_dims[0], frame_dims[1]], dtype=np.float32)

                for i in range(num_diff_stims):

                    # get stim times
                    if self.p['useStimOrder']:
                        use_stims = unique_stims[i::num_diff_stims]
                        indices = np.in1d(stim_order, use_stims)
                        indices = indices[:len(all_stim_times)]
                        stim_times = all_stim_times[indices[start_on_stim::only_every_x_stims]]
                    else:
                        # start:stop:step
                        stim_times = all_stim_times[start_on_stim + i::only_every_x_stims]
                    num_trials = len(stim_times)

                    # make and display status message
                    msg = 'Plane ' + str(z+1) + ' of ' + str(numZPlanes) + '. Stim ' + str(
                            i+1) + ' of ' + str(num_diff_stims) + ' (' + str(num_trials) + ' trials)'
                    print(msg)

                    # make frame indices
                    frames_with_stim = np.searchsorted(frame_times, stim_times)
                    frames_with_stim = frame_indices[frames_with_stim]

                    all_trials_sta_frames = []
                    for stim_frame_idx in frames_with_stim:
                        all_trials_sta_frames.append(sta_template + stim_frame_idx)

                    # get data
                    trials = np.zeros([len(sta_template), frame_dims[0], frame_dims[1], num_trials], dtype=np.float32)
                    for j, trial_sta_frames in enumerate(all_trials_sta_frames):
                        # print(trial_sta_frames)
                        for k, frame_idx in enumerate(trial_sta_frames):
                            trials[k, :, :, j] = movie[frame_idx]

                    print(msg + ' - Raw')

                    avg_movie = trials.mean(axis=3)
                    save_name = save_dir + os.sep + file_name_no_ext + \
                                '_Stim' + str(i+1) + '_STA_Plane' + \
                                str(z+1) + '_raw.tif'
                    tifffile.imsave(save_name, avg_movie.astype(np.uint16))

                    if methods:
                        for method in methods:
                            print(msg + ' - ' + method)

                            if self.p['useSingleTrials']:
                                norm_trials = sta.normalise_movie(trials, range(pre_samples), method=method)

                                if self.p['doThreshold']:
                                    norm_trials = sta.threshold(norm_trials, threshold=self.p['threshold'])

                                trial_imgs = sta.make_avg_image(norm_trials, avg_image_frames)
                                trial_imgs = trial_imgs.transpose([2, 0, 1])
                                save_name = save_dir + os.sep + file_name_no_ext + '_Stim' + \
                                        str(i+1) + '_STA_Plane' + str(z+1) + '_' + method + \
                                        '_All' + str(num_trials) + 'Trials' + '.tif'
                                tifffile.imsave(save_name, trial_imgs)

                            norm_movie = sta.normalise_movie(avg_movie, range(pre_samples), method=method)

                            if self.p['doThreshold']:
                                norm_movie = sta.threshold(norm_movie, threshold=self.p['threshold'])

                            save_name = save_dir + os.sep + file_name_no_ext + '_Stim' + \
                                    str(i+1) + '_STA_Plane' + \
                                    str(z+1) + '_' + method + '.tif'
                            tifffile.imsave(save_name, norm_movie)

                            avg_img[method][i] = sta.make_avg_image(norm_movie, avg_image_frames)
                            save_name = save_dir + os.sep + file_name_no_ext + '_Stim' + \
                                    str(i+1) + '_STA_Plane' + str(z+1) + \
                                    '_' + method + '_AvgImage' + '.tif'
                            tifffile.imsave(save_name, avg_img[method][i])

                            if self.p['colourByTime']:
                                results = sta.colour_by_time(norm_movie, avg_image_frames, smooth=int(
                                            frame_rate/2), useCorrelationImage=self.p['useCorrelationImage'],
                                            blurHandS=self.p['blurHandS'])
                                rgb = results['RGB']
                                hsv = results['HSV']
                                corr = results['Corr']

                                save_name = save_dir + os.sep + file_name_no_ext + '_Stim' + \
                                        str(i+1) + '_STA_Plane' + str(z+1) + \
                                        '_' + method + '_TimePeak' + '.tif'
                                tifffile.imsave(save_name, rgb)

                                save_name = save_dir + os.sep + file_name_no_ext + '_Stim' + \
                                        str(i+1) + '_STA_Plane' + str(z+1) + \
                                        '_' + method + '_Corr' + '.tif'
                                tifffile.imsave(save_name, (corr*65535).astype(np.uint16))

                                # tifffile.imsave(save_name.replace('.tif', '_H.tif'), H)
                                # tifffile.imsave(save_name.replace('.tif', '_S.tif'), S)
                                # tifffile.imsave(save_name.replace('.tif', '_V.tif'), V)

                if num_diff_stims > 1:
                    for method in methods:
                        if self.p['makeMaxImage']:
                            save_name = save_dir + os.sep + file_name_no_ext + '_' + \
                                    str(num_diff_stims) + 'Stims' + '_STA_Plane' + \
                                    str(z+1) + '_' + method + \
                                        '_MaxResponseImage' + '.tif'
                            tifffile.imsave(save_name, avg_img[method].max(axis=0))

                        if self.p['makeColourImage']:
                            msg = 'Making colour image'
                            print(msg)

                            # calculate indices
                            pref_stim = np.argmax(avg_img[method], axis=0).astype(np.int)

                            # this will only work with 4 visual orientations!
                            orth_stim = pref_stim + 2
                            orth_stim[orth_stim]
                            orth_stim[orth_stim >3] = orth_stim[orth_stim > 3]-3

                            x, y = np.meshgrid(np.arange(frame_dims[0]), np.arange(frame_dims[1]))

                            pref_img = avg_img[method][pref_stim, y, x]
                            pref_img[pref_img < 0] = 0
                            if num_diff_stims == 4:
                                orth_img = avg_img[method][orth_stim, y, x]
                                orth_img[orth_img < 0] = 0

                            other_imgs = avg_img[method].copy()
                            other_imgs[pref_stim, y, x] = np.nan
                            mean_other_imgs = np.nanmedian(other_imgs, axis=0).astype(np.float32)
                            if num_diff_stims == 4:
                                mean_other_imgs = orth_img

                            # hue
                            # from scipy.stats import vonmises
                            # data = avg_img[method]
                            # blank = np.zeros(frame_dims)
                            # for y_px in range(frame_dims[0]):
                            # 	for x_px in range(frame_dims[1]):
                            # 		this_data = data[:,y_px,x_px]
                            # 		print(this_data)
                            # 		kappa, loc, scale = vonmises.fit(data[:,y_px,x_px], fscale=1)
                            # 		blank[y_px, x_px] = kappa

                            H = pref_stim.astype(np.float32) / num_diff_stims
                            H[H == 0.0] = 0.025
                            H[H == 0.75] = 0.8
                            # H = H + 0.05
                            # H = H/2
                            # H = gaussian_filter(H,1)

                            # H = blank

                            # brightness
                            V = pref_img
                            # V = V / np.percentile(V, 99)
                            if method.lower() == 'df':
                                V = V/1000
                            elif method.lower() == 'dff':
                                V = V/100
                            elif method.lower() == 'zscore':
                                V = V/5

                            V[V < 0] = 0
                            V[V > 1] = 1

                            if self.p['useCorrelationImage']:
                                V = sta.makeCorrImg(sta.downsample(
                                        movie[frame_indices], 10), 4)
                                tifffile.imsave(save_name.replace(
                                        '.tif', '_Corr.tif'), (V*65535).astype(np.uint16))
                            v_min, v_max = np.percentile(V, (1, 99))
                            V = exposure.rescale_intensity(
                            V, in_range=(v_min, v_max))
                            # V = exposure.equalize_adapthist(V)
                            V = V/V.max()

                            # saturation
                            S = (pref_img - mean_other_imgs) / (pref_img + mean_other_imgs)
                            # S[V<np.nanpercentile(V,90)] = S[V<np.nanpercentile(V,90)]/10
                            # # S[S > np.percentile(S, 95)] = 0
                            # # S = S - np.percentile(S, 1)
                            # S = S / np.nanpercentile(S, 90)
                            # S = S*V
                            S[np.isnan(S)] = 0
                            # print(S.max())
                            # print(S.min())
                            # S = S*2
                            S[S < 0] = 0
                            S[S > 1] = 1
                            # S = gaussian_filter(S,1)

                            # convert HSV to RGB
                            # rgb_img = sta.hsv2rgb(H, S, V)

                            hsv = np.stack((H, S, V), axis=-1)
                            hsv_tf = tf.convert_to_tensor(hsv, np.float32)

                            rgb_tf = tf.image.hsv_to_rgb(hsv_tf)
                            sess = tf.Session()
                            
                            with sess.as_default():
                                rgb = rgb_tf.eval()

                            # blur the rgb image then convert back to hsv and keep the blurred H and S channels.
                            if self.p['blurHandS']:
                                rgb = gaussian_filter(rgb, (1, 1, 0))
                                rgb_tf = tf.convert_to_tensor(
                                        rgb, np.float32)
                                hsv_tf = tf.image.rgb_to_hsv(rgb_tf)

                                with sess.as_default():
                                    hsv2 = hsv_tf.eval()

                                hsv2[:, :, 2] = hsv[:, :, 2]
                                hsv_tf = tf.convert_to_tensor(
                                        hsv2, np.float32)
                                rgb_tf = tf.image.hsv_to_rgb(hsv_tf)
                                with sess.as_default():
                                    rgb = rgb_tf.eval()

                            rgb_img = (rgb * 65535).astype(np.uint16)

                            # tifffile.imsave(save_name.replace('.tif', 'pref_stim.tif'), pref_stim)
                            # tifffile.imsave(save_name.replace('.tif', 'mean_other_imgs.tif'), mean_other_imgs)
                            # tifffile.imsave(save_name.replace('.tif', 'H.tif'), H)
                            # tifffile.imsave(save_name.replace('.tif', 'S.tif'), S)
                            # tifffile.imsave(save_name.replace('.tif', 'V.tif'), (V*65535).astype(np.uint16))

                            # for v in range(num_diff_stims):
                            # 	idx = np.where(H==(np.float(v)/num_diff_stims))
                            # 	# V[idx] = 0.2 * V[idx] / np.mean(V[idx])
                            # 	single_img = np.zeros(frame_dims[0]*frame_dims[1], dtype=np.float)
                            # 	single_img[idx] = V[idx]
                            # 	single_img = np.reshape(single_img, [frame_dims[0],frame_dims[1]])
                            # 	save_name = save_dir + os.sep + file_name.replace('.bin', '_' + str(num_diff_stims) + 'stims_STA_PrefImage_' + str(v+1) + '.tif')
                            # 	tifffile.imsave(save_name, single_img.astype(np.int16))
                            # 	# V[idx] = V[idx] / np.max(V[idx])

                            save_name = save_dir + os.sep + file_name_no_ext + '_' + \
                                str(num_diff_stims) + 'Stims' + '_STA_Plane' + \
                                str(z+1) + '_' + method + '_PrefImage' + '.tif'
                            tifffile.imsave(save_name, rgb_img)

                            # hsv_img = np.stack([H,S,V], axis=2).astype(np.float16)
                            # print(hsv_img.shape)
                            # save_name = save_dir + os.sep + file_name_no_ext + '_' + str(num_diff_stims) + 'Stims' + '_STA_Plane' + str(z+1) + '_' + method + '_PrefImage_HSV' + '.tif'
                            # tifffile.imsave(save_name, hsv_img)