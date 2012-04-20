import os
from nipype.interfaces.base import (BaseInterface,
                                    BaseInterfaceInputSpec, File,
                                    Directory, traits, TraitedSpec)

from nipype.utils.misc import package_check
import warnings

from nipype import logging
iflogger = logging.getLogger('interface')

have_cmp = True
try:
    package_check('cmp')
except Exception, e:
    have_cmp = False
    warnings.warn('cmp not installed')
else:
    from cmp.util import runCmd

def generate_WM_and_GM_mask(subjects_dir, subject_id, wm_mask_filename, roi_filename):
    """Taken from cmtk and adapted for running as isolated function in
    this interface.

    http://github.com/tnez/cmp/blob/master/cmp/stages/parcellation/maskcreation.py

    """
    import os.path as op
    import nibabel as ni
    import numpy as np
    
    fs_dir = op.join(subjects_dir,subject_id)
    output_dir = op.abspath(op.curdir)
    
    # need to convert
    mri_cmd = 'mri_convert -i "%s/mri/aparc+aseg.mgz" -o "%s/mri/aparc+aseg.nii.gz"' % (fs_dir, fs_dir)
    runCmd( mri_cmd, None )

    fout = op.join(fs_dir, 'mri', 'aparc+aseg.nii.gz') ##OUTPUT
    niiAPARCimg = ni.load(fout)
    niiAPARCdata = niiAPARCimg.get_data()

    # mri_convert aparc+aseg.mgz aparc+aseg.nii.gz
    WMout = op.join(output_dir, wm_mask_filename)

    #%% label mapping
    # Using FreesurferColorLUT.txt
    # mappings are stored in mappings.ods

#    CORTICAL = {1 : [ 1, 2, 3, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],
#                2 : [31,13, 9,21,27,25,19,29,15,23, 1,24, 4,30,26,11, 6, 2, 5,22,16,14,10,20,12, 7, 8,18,30,17, 3,28,33]}
#
#
#    SUBCORTICAL = {1:[48,49,50,51,52,53,54,58,59,60, 9,10,11,12,13,17,18,26,27,28],
#                   2:[34,34,35,36,37,40,41,38,39,39,75,75,76,77,78,81,82,79,80,80]}
#
#    OTHER = {1:[16],
#             2:[83]}

    MAPPING = [[1,2012],[2,2019],[3,2032],[4,2014],[5,2020],[6,2018],[7,2027],[8,2028],[9,2003],[10,2024],
               [11,2017],[12,2026],[13,2002],[14,2023],[15,2010],[16,2022],[17,2031],[18,2029],[19,2008],
               [20,2025],[21,2005],[22,2021],[23,2011],[24,2013],[25,2007],[26,2016],[27,2006],[28,2033],
               [29,2009],[30,2015],[31,2001],[32,2030],[33,2034],[34,2035],[35,49],[36,50],[37,51],[38,52],
               [39,58],[40,53],[41,54],[42,1012],[43,1019],[44,1032],[45,1014],[46,1020],[47,1018],[48,1027],
               [49,1028],[50,1003],[51,1024],[52,1017],[53,1026],[54,1002],[55,1023],[56,1010],[57,1022],
               [58,1031],[59,1029],[60,1008],[61,1025],[62,1005],[63,1021],[64,1011],[65,1013],[66,1007],
               [67,1016],[68,1006],[69,1033],[70,1009],[71,1015],[72,1001],[73,1030],[74,1034],[75,1035],
               [76,10],[77,11],[78,12],[79,13],[80,26],[81,17],[82,18],[83,16]]

    WM = [2, 29, 32, 41, 61, 64, 59, 60, 27, 28] + range(77,86+1) + range(100, 117+1) + \
         range(155,158+1) + range(195,196+1) + range(199,200+1) + range(203,204+1) + \
         [212, 219, 223] + range(250,255+1)
    # add
    # 59  Right-Substancia-Nigra
    # 60  Right-VentralDC
    # 27  Left-Substancia-Nigra
    # 28  Left-VentralDC

    #%% create WM mask    
    niiWM = np.zeros( niiAPARCdata.shape, dtype = np.uint8 )

    for i in WM:
         niiWM[niiAPARCdata == i] = 1

    # we do not add subcortical regions
#    for i in SUBCORTICAL[1]:
#         niiWM[niiAPARCdata == i] = 1

    img = ni.Nifti1Image(niiWM, niiAPARCimg.get_affine(), niiAPARCimg.get_header())
    ni.save(img, WMout)

    #%% create GM mask (CORTICAL+SUBCORTICAL)
    #%  -------------------------------------
    GMout = op.join(output_dir, roi_filename)
    niiGM = np.zeros( niiAPARCdata.shape, dtype = np.uint8 )
    for ma in MAPPING:
        niiGM[ niiAPARCdata == ma[1]] = ma[0]
#        # % 33 cortical regions (stored in the order of "parcel33")
#        for idx,i in enumerate(CORTICAL[1]):
#            niiGM[ niiAPARCdata == (2000+i)] = CORTICAL[2][idx] # RIGHT
#            niiGM[ niiAPARCdata == (1000+i)] = CORTICAL[2][idx] + 41 # LEFT
#        #% subcortical nuclei
#        for idx,i in enumerate(SUBCORTICAL[1]):
#            niiGM[ niiAPARCdata == i ] = SUBCORTICAL[2][idx]
#        # % other region to account for in the GM
#        for idx, i in enumerate(OTHER[1]):
#            niiGM[ niiAPARCdata == i ] = OTHER[2][idx]
        img = ni.Nifti1Image(niiGM, niiAPARCimg.get_affine(), niiAPARCimg.get_header())
        ni.save(img, GMout)

def crop_and_move_WM_and_GM(subjects_dir, subject_id, wm_mask_filename, roi_filename):
    """Taken from cmtk and adapted for running as isolated function in
    this interface.

    http://github.com/tnez/cmp/blob/master/cmp/stages/parcellation/maskcreation.py

    """
    import os.path as op
    import nibabel as ni
    import numpy as np

    fs_dir = op.join(subjects_dir,subject_id)
    output_dir = op.abspath(op.curdir)

    # datasets to crop and move: (from, to)
    ds = [
          (op.join(fs_dir, 'mri', 'fsmask_1mm.nii.gz'), op.join(output_dir, 'fsmask_1mm.nii.gz') ),
          ]

    ds.append( (op.join(fs_dir, 'mri', 'ROIv_%s.nii.gz' % p), op.join(ouput_dir, 'ROIv_HR_th_%s.nii.gz')) )

    orig = op.join(fs_dir, 'mri', 'orig', '001.mgz')

    for d in ds:
        log.info("Processing %s:" % d[0])

        # does it exist at all?
        if not op.exists(d[0]):
            raise Exception('File %s does not exist.' % d[0])
        # reslice to original volume because the roi creation with freesurfer
        # changed to 256x256x256 resolution
        mri_cmd = 'mri_convert -rl "%s" -rt nearest "%s" -nc "%s"' % (orig, d[0], d[1])
        runCmd( mri_cmd,log )
    
class ROIGenInputSpec(BaseInterfaceInputSpec):
    subjects_dir = Directory(exists=True, mandatory=True,
                             desc='freesurfer subjects directory for given run')
    subject_id = traits.Str(mandatory=True, desc='freesurfer subject id (this will be FREESURFER when using acmtk directory structure)')
    roi_filename = File("freesurferaparc.nii.gz", desc='Filename of ROI file to be generated',
                        genfile=True, usedefault=True)
    wm_mask_filename = File("fs_mask_1mm.nii.gz", desc='Filename of white-matter mask to be generated',
                            genfile=True, usesdefault=True)

class ROIGenOutputSpec(TraitedSpec):
    roi_file = File(desc="Generated ROI volume in structural space", exists=True)
    wm_mask_file = File(desc="Generated WM mask in structural space", exists=True)

class ROIGen(BaseInterface):
    """Generates an ROI volume given a FREESURFER dir that has been
    run through recon-all

    Example
    -------
    >>> import muscip.interfaces.freesurfer as muscfree
    >>> rg = muscfree.ROIGen()
    >>> rg.inputs.subjects_dir = 'path/to/cmtk/style/subject'
    >>> rg.inputs.subject_id = 'FREESURFER' # given cmtk style dir structure
    >>> rg.inputs.roi_filename = 'freesurferaparc.nii.gz'
    >>> rg.run() # doctest: +SKIP

    """

    input_spec = ROIGenInputSpec
    output_spec = ROIGenOutputSpec

    def _run_interface(self, runtime):
        iflogger.info('Generating WM and GM masks...')
        generate_WM_and_GM_mask()
        iflogger.info('Cropping and moving WM and GM masks...')
        crop_and_move_WM_and_GM_mask()
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['roi_file'] = self.inputs.roi_filename
        outputs['wm_mask_file'] = self.inputs.wm_mask_filename
        return outputs
