# ~ based off code originally written by Adam R. (in the Josselyn Lab ~ 2023)

"""
!!!! READ FIRST !!!! ~ HOW TO USE THIS SCRIPT

As of now, this script can be run with two options. Note that this script is also not tested yet, so there maybe be errors that you should fix and git commit/push please :)

Option 1) use for converting a set of raw data (i.e. pairs of .mat + .sbx files) to .tiff files
    This is run as:
    1) at the command line, cd to the location of your raw data
    2) from the command line, run: `python <enter-path-to-this-script>`
    3) this should then run this script and create .tiff files (as long as there are no other errors)

Option 2) use for creating and saving an ImagingPlus experiment object to use for further analysis with this package.
- NOTE: this options doesn't work yet - the code isn't built out for it - please contribute to it below to make it functional!
    - to enable this, we will use the createExperiment flag and set it to true as shown:
    1) at the command line, cd to the location of your raw data
    2) from the command line, run: `python  <enter-path-to-this-script> createExperiment=True`
    3) this should then run this script and create .tiff files, and also create an Experiment object (as long as there are no other errors)
    - note the Experiment object will be saved as a .pkl file in the same directory

"""


import pathlib
from dataclasses import dataclass, field
from datetime import datetime
import os

import numpy as np
# import skimage.io as io
from scipy import io

from imagingplus import Experiment
from imagingplus.main.core import ImagingTrial
from imagingplus.utils.images import WriteTiff
from imagingplus.utils.utils import get_args

start_time = datetime.now()


@dataclass
class Neurolabware:
    sbx_list: list[pathlib.Path] = field(default_factory=list)  #: list of sbx paths in current batch
    tiff_list: list[str] = field(default_factory=list)  #: list of tiff paths (after conversion) in current batch

    def __post_init__(self):
        nbFiles = 0
        for i in self.sbx_list:
            fpath = i.as_posix().split('.')[0]
            # save_name_base = '_'.join([i.parts[-4], i.parts[-3], i.stem])
            save_name_base = i.parts[-1][:-4]
            save_name = i.parent.joinpath(save_name_base).as_posix()

            print('Converting ' + fpath)
            self.sbx_to_tiff(fpath, save_name)
            nbFiles += 1

        print(f'Conversion of {nbFiles} .sbx files to TIFF completed, {(datetime.now() - start_time)} seconds elapsed --- ')

    @staticmethod
    def _loadmatfile(matfile):
        info = io.loadmat(matfile)['info']

        # The following variables are taken from the corresponding scanbox code
        info_channels = info['channels'][0][0][0]
        if info_channels == 1:
            info_nchan = 2
            factor = 1
        if info_channels == 2 or info_channels == 3:
            info_nchan = 1
            factor = 2

        info_sz = info['sz'][0][0][0].astype(np.uint16)
        info_scanmode = info['scanmode'][0][0]

        info_bytesperbuffer = info['bytesPerBuffer'][0][0][0][0]
        info_recordsperbuffer = info['recordsPerBuffer'][0][0][0][0]
        info_scanbox_version = info['scanbox_version'][0][0][0]

        # Computed values
        # info_nsamples = info_sz[1].astype(int) * info_recordsperbuffer * 2 * info_nchan
        # info_max_idx =  fid.size/info_recordsperbuffer/info_sz[1]*factor/4 - 1

        return info_sz

    def sbx_to_tiff(self, fpath, s_name):
        start = datetime.now()
        print('--- Start time = ', start)

        print('--- Loading matfile')
        try:
            info_sz = self._loadmatfile(fpath + '.mat')
        except:
            pass  # TODO : armaan : remove this `pass` and raise the appropriate Error message when the .mat file is not found in the same folder location as the .sbx file; hint: use FileNotFoundError() and add a helpful message for user to identify what the issue is
        height, width = info_sz
        maxint = np.iinfo(
            np.uint16).max  # This is the largest value a uint16 can store. We have to remove it from each value in the binary, for some reason.
        print("--- Height = {}, Width = {}".format(height, width))

        print('--- Reading sbx file to memory...')
        fid = np.fromfile(fpath + '.sbx', dtype=np.uint16)

        print('--- Done reading. Reshaping array...')
        vid = fid.reshape(-1, height, width)
        vid = maxint - vid
        fid = None

        WriteTiff(save_path=s_name + '.tiff', stack=vid)
        self.tiff_list.append(s_name + '.tiff')
        # print('--- Saving as TIFF...')
        # io.imsave(s_name + '.tiff', vid)
        # print('--- Saving complete')
        end = datetime.now()
        print('--- End time = ', end)
        print('--- Runtime = ', end - start)

    @classmethod
    def newExperimentFromNeurolabware(cls, sbx: list[pathlib.Path], date: str, dataPath: str, expID: str, saveDir: str,
                                      expobj_pkl: str = None):
        """
        Alternative Constructor:
        Create a new experiment object (including imaging trials) for use in analysis from Neurolabware data.

        :param sbx: list of pathnames for .sbx data files to import into imagingplus (note that matfiles are also required)
        :param expobj_pkl: if using a preexisting imagingplus Experiment object, provide .pkl file path
        :param dataPath:
        :param expID:
        :param saveDir:
        :param kwargs: see core.Experiment for additional optional arguments.
        """

        nb_exp = cls(sbx_list=sbx)
        if expobj_pkl is None:
            expobj = Experiment(dataPath= dataPath,
                                expID= expID,
                                saveDir= saveDir)
        else:
            from imagingplus.utils.io import import_obj
            expobj = import_obj(expobj_pkl)
        for tiff in nb_exp.tiff_list:
            trialID = tiff.split('/')[-1].split('.')[0]
            trialobj = ImagingTrial.newImagingTrialfromExperiment(experiment=expobj, trialID=trialID, dataPath=tiff,
                                                                  date=date)

        expobj.save()
        return nb_exp, expobj


if __name__ == "__main__":
    args = get_args()
    # folderpath = args.pathname

    folderpath = os.getcwd()

    exten = '*.sbx'  # CHANGE THIS IF YOU WANT TO CONVERT DIFFERENT VIDEOS
    print(folderpath)

    p = pathlib.Path(folderpath)

    sbx_list = p.rglob(exten)

    # process neurolabware experiment with the sbx list of files to process in current batch
    if 'createExperiment' in args and args['createExperiment'] is True:
        assert ('date' in args and 'dataPath' in args and 'expID' in args and
                'saveDir' in args), 'missing data or dataPath or expId or saveDir - all arguments are required to create an iamgingplus Experiment'
        _, _ = Neurolabware.newExperimentFromNeurolabware(sbx=list(sbx_list), date=args['date'],
                                                          dataPath=args['dataPath'], expID=args['expID'],
                                                          saveDir=args['saveDir'])
    else:
        _ = Neurolabware(sbx_list=list(sbx_list))





