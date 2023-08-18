# ~ based off code originally written by Adam R. (in the Josselyn Lab ~ 2023)

import pathlib
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
from scipy import io

from imagingplus import Experiment
from imagingplus.main.core import ImagingTrial
from imagingplus.utils.images import WriteTiff
from imagingplus.utils.utils import get_args


@dataclass
class Neurolabware:
    sbx_list: list[pathlib.Path] = field(default_factory=list)  #: list of sbx paths in current batch
    tiff_list: list[str] = field(default_factory=list)  #: list of tiff paths (after conversion) in current batch

    def __post_init__(self):
        nbFiles = 0
        for i in sbx_list:
            fpath = i.as_posix().split('.')[0]
            save_name_base = '_'.join([i.parts[-4], i.parts[-3], i.stem])
            save_name = i.parent.joinpath(save_name_base).as_posix()

            print('Converting ' + fpath)
            print('Saving as ' + save_name_base)
            self.sbx_to_tiff(fpath, save_name)
            nbFiles += 1

        print(f'Conversion of {nbFiles} files over!')
        print((datetime.now() - start_time), "seconds elapsed ---")
        print('Done!')

    @staticmethod
    def _loadmatfile(matfile):
        info = io.loadmatfile(matfile)['info']

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
        info_sz = self._loadmatfile(fpath + '.mat')
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
    def newExperimentFromNeurolabware(cls, sbx: list[pathlib.Path], **kwargs):
        """
        Alternative Constructor:
        Create a new experiment object (including imaging trials) for use in analysis from Neurolabware data.

        :param kwargs: see core.Experiment for required and optional arguments.
        """

        nb_exp = cls(sbx_list=sbx)

        expobj = Experiment(dataPath=kwargs['dataPath'],
                            expID=kwargs['expID'],
                            saveDir=kwargs['saveDir'],
                            **kwargs)
        for tiff in nb_exp.tiff_list:
            trialID = tiff.split('.')[0]
            trialobj = ImagingTrial.newImagingTrialfromExperiment(experiment=expobj, trialID=trialID, **kwargs)

        expobj.save()
        return nb_exp, expobj


if __name__ == "__main__":
    args = get_args()
    folderpath = args.pathname

    exten = '*.sbx'  # CHANGE THIS IF YOU WANT TO CONVERT DIFFERENT VIDEOS
    print(folderpath)

    start_time = datetime.now()

    p = pathlib.Path(folderpath)

    sbx_list = p.rglob(exten)

    # process neurolabware experiment with the sbx list of files to process in current batch
    if 'createExperiment' in args and args['createExperiment'] is True:
        _, _ = Neurolabware.newExperimentFromNeurolabware(sbx=list(sbx_list))
    else:
        _ = Neurolabware(sbx_list=list(sbx_list))





