# image processing
import os
from typing import Union

import cv2 as cv
import numpy as np
import tifffile as tf
from matplotlib import pyplot as plt, patches
from skimage import io as skio


# simple ZProfile function for any sized square in the frame (equivalent to ZProfile function in Fiji)
def ZProfile(movie, area_center_coords: tuple = None, area_size: int = -1, plot_trace: bool = True,
             plot_image: bool = True, plot_frame: int = 1, vasc_image: np.array = None, **kwargs):
    """
    Plot a z-profile of a movie, averaged over space inside a square area. by S. Armstrong.

    :param movie: can be np.array of the TIFF stack or a tiff path from which it is read in
    :param area_center_coords: coordinates of pixel at center of box (x,y)
    :param area_size: int, length and width of the square in pixels
    :param plot_frame: which movie frame to take as a reference to plot the area boundaries on
    :param vasc_image: optionally include a vasculature image tif of the correct dimensions to plot the coordinates on.

TODO add paramters
    :param plot_trace:
    :param plot_image:
    :param kwargs:
    :return:
    """

    if type(movie) is str:
        movie = tf.imread(movie)
    print('plotting zprofile for TIFF of shape: ', movie.shape)

    # assume 15fps for 1024x1024 movies and 30fps imaging for 512x512 movies
    if movie.shape[1] == 1024:
        img_fps = 15
    elif movie.shape[1] == 512:
        img_fps = 30
    else:
        img_fps = None

    assert area_size <= movie.shape[1] and area_size <= movie.shape[2], "area_size must be smaller than the image"
    if area_size == -1:  # this parameter used to plot whole FOV area
        area_size = movie.shape[1]
        area_center_coords = (movie.shape[1] / 2, movie.shape[2] / 2)
    assert area_size % 2 == 0, "pls give an even area size"

    x = area_center_coords[0]
    y = area_center_coords[1]
    x1 = int(x - 1 / 2 * area_size)
    x2 = int(x + 1 / 2 * area_size)
    y1 = int(y - 1 / 2 * area_size)
    y2 = int(y + 1 / 2 * area_size)
    smol_movie = movie[:, y1:y2, x1:x2]
    smol_mean = np.nanmean(smol_movie, axis=(1, 2))
    print('|- Output shape for z profile: ', smol_mean.shape)

    if plot_image:
        f, ax1 = plt.subplots()
        ref_frame = movie[plot_frame, :, :]
        if vasc_image is not None:
            assert vasc_image.shape == movie.shape[1:], 'vasculature image has incompatible dimensions'
            ax1.imshow(vasc_image, cmap="binary_r")
        else:
            ax1.imshow(ref_frame, cmap="binary_r")

        rect1 = patches.Rectangle(
            (x1, y1), area_size, area_size, linewidth=1.5, edgecolor='r', facecolor="none")

        ax1.add_patch(rect1)
        ax1.set_title("Z-profile area")
        plt.show()

    if plot_trace:
        if 'figsize' in kwargs:
            figsize = kwargs['figsize']
        else:
            figsize = [10, 4]
        fig, ax2 = plt.subplots(figsize=figsize)
        if img_fps is not None:
            ax2.plot(np.arange(smol_mean.shape[0]) / img_fps, smol_mean, linewidth=0.5, color='black')
            ax2.set_xlabel('Time (sec)')
        else:
            ax2.plot(smol_mean, linewidth=0.5, color='black')
            ax2.set_xlabel('frames')
        ax2.set_ylabel('Flu (a.u.)')
        if 'title' in kwargs:
            ax2.set_title(kwargs['title'])
        plt.show()

    return smol_mean

def subselect_tiff(tiff_path: str = None, tiff_stack: np.array = None, select_frames: tuple = (0, 0),
                   save_as: str = None):
    """
TODO fill explanation and add parameters
    :param tiff_path:
    :param tiff_stack:
    :param select_frames:
    :param save_as:
    :return:
    """
    if tiff_stack is None:
        # open tiff file
        print('running subselecting tiffs')
        print('|- working on... %s' % tiff_path)
        tiff_stack = ImportTiff(tiff_path)

    stack_cropped = tiff_stack[select_frames[0]:select_frames[1]]

    if save_as is not None:
        tf.imwrite(save_as, stack_cropped, photometric='minisblack')

    return stack_cropped


def convert_to_8bit(img, target_type_min=0, target_type_max=255):
    """
    TODO fill documentation and add parameters
    :param img:
    :param target_type:
    :param target_type_min:
    :param target_type_max:
    :return:
    """
    imin = img.min()
    imax = img.max()

    a = (target_type_max - target_type_min) / (imax - imin)
    b = target_type_max - a * imax
    new_img = (a * img + b).astype(np.uint8)
    return new_img


def fourierImage(img: np.ndarray):
    """
TODO fill documentation and add parameters
    :param img:
    :return:
    """
    dft = cv.dft(np.float32(img), flags=cv.DFT_COMPLEX_OUTPUT)
    dft_shift = np.fft.fftshift(dft)
    return dft_shift


def highPass(img: np.ndarray, fshift, filter=30):
    """high pass filtering of image"""
    rows, cols = img.shape
    crow, ccol = rows // 2, cols // 2

    fshift[crow - filter:crow + (filter + 1), ccol - filter:ccol + (filter + 1)] = 0
    f_ishift = np.fft.ifftshift(fshift)
    img_back = np.fft.ifft2(f_ishift)
    img_back = np.real(img_back)
    return img_back


def thresholdImage(img: np.ndarray, low_threshold: int = 0, high_threshold: int = 250):
    filter_out = np.where((low_threshold > img) | (img > high_threshold))
    img[filter_out] = 0
    return img


def lowPass(img: np.ndarray, fshift, filter=30):
    rows, cols = img.shape
    crow, ccol = rows // 2, cols // 2

    # create a mask first, center square is 1, remaining all zeros
    mask = np.zeros((rows, cols, 2), np.uint8)
    mask[crow - filter:crow + filter, ccol - filter:ccol + filter] = 1

    # apply mask and inverse DFT
    fshift = fshift * mask
    f_ishift = np.fft.ifftshift(fshift)
    img_back = cv.idft(f_ishift)
    img_back = cv.magnitude(img_back[:, :, 0], img_back[:, :, 1])

    return img_back


def bandPass(img: np.ndarray, fshift, lowfilter=10, highfilter=3):
    rows, cols = img.shape
    crow, ccol = rows // 2, cols // 2

    # create a mask first, center square is 1, remaining all zeros
    mask = np.zeros((rows, cols, 2), np.uint8)
    mask[crow - (lowfilter - highfilter):crow + (lowfilter - highfilter),
    ccol - (lowfilter - highfilter):ccol + (lowfilter - highfilter)] = 1

    # apply mask and inverse DFT
    fshift = fshift * mask
    f_ishift = np.fft.ifftshift(fshift)
    img_back = cv.idft(f_ishift)
    img_back = cv.magnitude(img_back[:, :, 0], img_back[:, :, 1])

    return img_back

def z_score_img(key_image: np.ndarray, mean_img: np.ndarray = None, std_img: np.ndarray = None, normal_img_stack: np.ndarray = None, plot: bool = False):
    """
    z-score normalization of the input key_image to the input mean img and std img, or to image stack provided in normal_img_stack.

    :param key_image: input image
    :param mean_img: input mean image
    :param std_img: input std image
    :param normal_img_stack: input stack of images to normalize the key image to.
    :return:
    """
    if not mean_img or not std_img:
        assert normal_img_stack is not None, 'must provide normal_img_stack to use as baseline to z-normalize to.'
        assert normal_img_stack.shape[0] > 1, 'normal_img_stack must be a stack of 2-D arrays.'
        mean_img = np.mean(normal_img_stack, axis=0)
        std_img = np.std(normal_img_stack, axis=0)
        assert key_image.shape == mean_img.shape == std_img.shape, 'Shape of img must match mean_img and std_img'
    else:
        assert key_image.shape == mean_img.shape == std_img.shape, 'Shape of img must match mean_img and std_img'

    z_scored = (key_image - mean_img) / std_img
    if plot: plt.imshow(z_scored, cmap='bwr'), plt.colorbar(), plt.show()
    return z_scored


def makeFrameAverageTiff(frames: Union[int, list, tuple], tiff_path: str = None, stack: np.ndarray = None,
                         peri_frames: int = 100, save_dir: str = None, to_plot=False, **kwargs):
    """Creates, plots and/or saves an average image of the specified number of peri-key_frames around the given frame from either the provided tiff_path or the stack array.
    TODO add parameters
    :param frames:
    :param tiff_path:
    :param stack:
    :param peri_frames:
    :param save_dir:
    :param to_plot:
    :param kwargs:
    :return:
    """

    if type(frames) == int:
        frames = [frames]

    stack = ImportTiff(tiff_path) if not stack else stack

    imgs = []
    for idx, frame in enumerate(frames):
        # im_batch_reg = tf.imread(tif_path, key=range(0, self.output_ops['batch_size']))

        if 0 > frame - peri_frames // 2:
            peri_frames_low = frame
        else:
            peri_frames_low = peri_frames // 2
        if stack.shape[0] < frame + peri_frames // 2:
            peri_frames_high = stack.shape[0] - frame
        else:
            peri_frames_high = peri_frames // 2
        im_sub_reg = stack[frame - peri_frames_low: frame + peri_frames_high]

        avg_sub = np.mean(im_sub_reg, axis=0)

        # convert to 8-bit
        avg_sub = convert_to_8bit(avg_sub, 0, 255)

        if save_dir:
            if '.tif' in save_dir: save_dir = os.path.dirname(save_dir) + '/'
            save_path = save_dir + f'/{frames[idx]}_s2preg_frame_avg.tif'
            os.makedirs(save_dir, exist_ok=True)

            print(f"\t\- Saving averaged s2p registered tiff for frame: {frames[idx]}, to: {save_path}")
            tf.imwrite(save_path, avg_sub, photometric='minisblack')

        if to_plot:  # todo replace with proper plotting average tiff frame code
            kwargs['title'] = f'{peri_frames} key_frames avg from s2p reg tif, frame: {frames[idx]}' if not kwargs[
                'title'] else kwargs['title']
            from packerlabimaging.plotting.plotting import plotImg
            plotImg(img=avg_sub, **kwargs)

        imgs.append(avg_sub)

    return np.asarray(imgs)


def ImportTiff(tiff_path, frames: Union[tuple, list, int] = None):
    """
TODO fill explanation and add parameters
    :param tiff_path:
    :param frames:
    :return:
    """
    if frames and type(frames) == tuple:
        im_stack = tf.imread(tiff_path, key=range(frames[0], frames[1]))
    elif frames and type(frames) == int:
        im_stack = tf.imread(tiff_path, key=frames)
    else:
        # import cv2
        # ret, images = cv2.imreadmulti(tiff_path, [], cv2.IMREAD_ANYCOLOR)
        # if len(images) > 0:
        #     im_stack = np.asarray(images)
        # im_stack = tf.imread(tiff_path)
        try:
            stack = []
            with tf.TiffFile(tiff_path) as tif:
                for page in tif.pages:
                    image = page.asarray()
                    stack.append(image)
            im_stack = np.array(stack)
            if len(im_stack) == 1: im_stack = im_stack[0]
        except Exception as ex:
            try:
                im_stack = skio.imread(tiff_path, plugin='pil')
            except Exception as ex:
                raise ImportError('unknown error in loading tiff stack.')

    return im_stack


def SaveDownsampledTiff(tiff_path: str = None, stack: np.array = None, group_by: int = 4, save_as: str = None):
    """
    Create and save a downsampled version of the original tiff file. Original tiff file can be given as a numpy array stack
    or a str path to the tiff.

    :param tiff_path: path to the tiff to downsample
    :param stack: numpy array stack of the tiff file already read in
    :param group_by: specified interval for grouped averaging of the TIFF
    :param save_as: .tif path to save the downsampled tiff to, if none provided it will save to the same parent directory as the provided tiff_path
    :param plot_zprofile: if True, plot the zaxis profile using the full TIFF stack provided.
    :return: numpy array containing the downsampled TIFF stack
    """
    print(f'\- downsampling tiff stack {stack.shape}...') if stack is not None else None

    if save_as is None:
        assert tiff_path is not None, "please provide a save path to save_as"
        save_as = tiff_path[:-4] + f'{group_by}x_downsampled.tif'

    if stack is None and tiff_path is not None:
        # open tiff file
        print('|- working on... %s' % tiff_path)
        stack = ImportTiff(tiff_path)
        print(f'\- downsampling tiff stack [{stack.shape}]...') if stack else None

    resolution = stack.shape[1]

    # downsample to 8-bit
    stack8 = np.full_like(stack, fill_value=0)
    for frame in np.arange(stack.shape[0]):
        stack8[frame] = convert_to_8bit(stack[frame], 0, 255)

    # stack8 = stack

    # grouped average by specified interval
    num_frames = stack8.shape[0] // group_by
    # avgd_stack = np.empty((num_frames, resolution, resolution), dtype='uint16')
    avgd_stack = np.empty((num_frames, resolution, resolution), dtype='uint8')
    frame_count = np.arange(0, stack8.shape[0], group_by)
    for i in np.arange(num_frames):
        frame = frame_count[i]
        avgd_stack[i] = np.mean(stack8[frame:frame + group_by], axis=0)

    avgd_stack = avgd_stack.astype(np.uint8)

    # bin down to 512 x 512 resolution if higher resolution - not functional so far
    shape = np.shape(avgd_stack)
    if shape[1] != 512:
        # input_size = avgd_stack.shape[1]
        # output_size = 512
        # bin_size = input_size // output_size
        # final_stack = avgd_stack.reshape((shape[0], output_size, bin_size,
        #                                   output_size, bin_size)).mean(4).mean(2)
        final_stack = avgd_stack
    else:
        final_stack = avgd_stack

    # write output
    # print(f"\n\- saving [{final_stack.shape}] tiff to... {save_as}")
    WriteTiff(save_path=save_as, stack=final_stack) if save_as else None

    return final_stack


def WriteTiff(save_path, stack: np.array):
    """use Tifffile imwrite function to save a numpy array to tiff file
    TODO add parameters
    :param save_path:
    :param stack:
    """
    print(f"\n\- saving array [{stack.shape}] to: {save_path}", end="\r")
    if not os.path.exists(os.path.dirname(save_path)):
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
    tf.imwrite(file=save_path, data=stack, photometric='minisblack')
    print(f"\n|- saved array [{stack.shape}] to: {save_path}")


def make_tiff_stack(tiff_paths: list, save_as: str = None) -> np.ndarray:
    """
    Read in a bunch of tiffs and stack them together, and save the output as the save_as

    :param tiff_paths:
    :return:

    TODO old:
    :param sorted_paths: ls of string paths for tiffs to stack
    :param save_as: .tif file path to where the tif should be saved
    """

    num_tiffs = len(tiff_paths)
    # TEST_TIFFS = tiff_paths[:5]
    # tiff_paths = TEST_TIFFS
    print('working on tifs to stack: ', num_tiffs)

    data = []
    for i, tif_ in enumerate(tiff_paths):
        img = ImportTiff(tiff_path=tif_)
        # print(f'imported tif stack: {img.shape}')
        data.extend(img)
        # print('length of data: ', len(data))
    data = np.array(data)
    print(f'\- made tiff stack: {data.shape}')
    WriteTiff(save_path=save_as, stack=data) if save_as else None
    return data

