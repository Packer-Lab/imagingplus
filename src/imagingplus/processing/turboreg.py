# TODO TURN INTO PROPER CODE THAT CAN BE INTEGRATED INTO THE PACKAGE WORKFLOW

# run python implementation of turboreg
import numpy as np
from matplotlib import pyplot as plt

from imagingplus.main.core import ImagingTrial
from pystackreg import StackReg
from skimage import io

from imagingplus.utils.images import ImportTiff


def run__turboreg(stack_path, type: str, **kwargs):

    # # plot reference image
    # fig, ax = plt.subplots(figsize=(6,6), dpi=300)
    # ax.imshow(ref_image)
    # ax.title('reference image')
    # fig.show()

    data_tiff = ImportTiff(stack_path, frames=(0, 1000))  # 3 dimensions : frames x width x height
    print(data_tiff.shape)

    sr = StackReg(StackReg.TRANSLATION)

    # todo need to figure out how to register to a custom reference image!!!
    # register each frame to an input reference image
    out_tra = sr.register_transform_stack(data_tiff, ref_image)

    out_first30 = sr.register_transform_stack(data_tiff, reference='first', n_frames=30, verbose=True)
    print(out_first30.shape)


# %% 1) register and get transformation matrices for red channel stack in trials where both red and green channel imaging was used

sr = StackReg(StackReg.AFFINE)

multi_tiff_path_reg_transform = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-004/2021-01-28_PS14_t-004_Cycle00001_Ch2.tif'
multi_tiff_path_transform = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-004/2021-01-28_PS14_t-004_Cycle00001_Ch3.tif'

reg_transform = ImportTiff(multi_tiff_path_reg_transform, frames=(0, 1000))
transform = ImportTiff(multi_tiff_path_transform, frames=(0, 1000))


# %% register
tmats = sr.register_stack(reg_transform, reference='first', n_frames = 60, verbose=True)

# transform both
reg_transform_out = sr.transform_stack(reg_transform)
transform_out = sr.transform_stack(transform)

# %% show registered images
plt.imshow(np.mean(reg_transform[:60], axis=0), cmap='Greys_r')
plt.title('first 60fr - reg + transform')
plt.show()


plt.imshow(np.mean(reg_transform_out[:60], axis=0), cmap='Greys_r')
plt.title('first 60fr - reg + transform OUT')
plt.show()

plt.imshow(np.mean(reg_transform_out[-60:], axis=0), cmap='Greys_r')
plt.title('last 60fr - reg + transform OUT')
plt.show()


plt.imshow(np.mean(transform[:60], axis=0), cmap='Greys_r')
plt.title('first 60fr - transform')
plt.show()

plt.imshow(np.mean(transform_out[:60], axis=0), cmap='Greys_r')
plt.title('first 60fr - transform OUT')
plt.show()

plt.imshow(np.mean(transform_out[-60:], axis=0), cmap='Greys_r')
plt.title('last 60fr - transform OUT')
plt.show()


# %% back register + transform the single scan red ref image
red_scan = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_s-013/2021-01-28_PS14_s-013_Cycle00001_Ch2_000001.ome.tif'

red_img = ImportTiff(red_scan)[0]

sr = StackReg(StackReg.AFFINE)
red_img_out = sr.register_transform(np.mean(reg_transform_out[:60], axis=0), red_img)


plt.imshow(red_img, cmap='Greys_r')
plt.title('red img un reg')
plt.show()

plt.imshow(red_img_out, cmap='Greys_r')
plt.title('red img back registered')
plt.show()






# %%
ref_image_path = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_s-013/2021-01-28_PS14_s-013_Cycle00001_Ch2_000001.ome.tif'
multi_tiff_path = '/home/pshah/mnt/qnap/Data/2021-01-28/2021-01-28_PS14_t-004/2021-01-28_PS14_t-004_Cycle00001_Ch2.tif'

ref_img = ImportTiff(ref_image_path)[0]
data_tiff = ImportTiff(multi_tiff_path, frames=(0, 1000))  # 3 dimensions : frames x width x height
print(data_tiff.shape)
unreg_data = np.mean(data_tiff[:30], axis=0)

plt.imshow(unreg_data, cmap='Greys_r')
plt.title('first 30fr')
plt.show()

plt.imshow(ref_img, cmap='Greys_r')
plt.title('single scan')
plt.show()

sr = StackReg(StackReg.RIGID_BODY)
out_first30 = sr.register_transform_stack(data_tiff, reference='first', n_frames=30, verbose=True)
print(out_first30.shape)

reg_30 = np.mean(out_first30[:30], axis=0)
plt.imshow(reg_30, cmap='Greys_r')
plt.title('reg first 30fr')
plt.show()



# %% register first 30 data frames to ref_image

unreg_data = np.mean(data_tiff[:30], axis=0)
sr = StackReg(StackReg.RIGID_BODY)
tmat = sr.register(ref_img, unreg_data)
out = sr.transform(unreg_data, tmat)
# look at out - and see if it has been registered to ref_img
plt.imshow(unreg_data, cmap='Greys_r')
plt.title('first 30fr - unreg')
plt.show()

plt.imshow(out, cmap='Greys_r')
plt.title('first 30fr - reg to ref_img')
plt.show()

plt.imshow(ref_img, cmap='Greys_r')
plt.title('single scan')
plt.show()

# transform full stack to this registration mat
full_out = sr.transform_stack(data_tiff, tmats=tmat) # - doesn't work... (need same length of tmat as data to be transformed)
















