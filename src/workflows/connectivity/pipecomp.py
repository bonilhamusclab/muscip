"""Temporary workflow built to accomplish pipe comparisons -- either
replacement or heavy refactoring intented for future.

"""

import os

def create_pipecomp_workflow(base_dir, subject_ids, name='pipecomp_workflow'):
    """Process structural diffusion with the intention of comparing
    results - this is a temporary workflow.

    Example
    -------

    >>> import tn_image_processing.workflows.connectivity.pipecomp as pipes
    >>> W = pipes.create_pipecomp_workflow('~/path/to/base/dir',
    >>>                                    ['s01','s02'])
    >>> W.run()

    Inputs::

      [Mandatory]
      base_dir: directory containing all subjects
      subject_ids: list of strings representing subject dirs to be processed

    Outputs::

      None, but files written to disk
    
    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio

    this_workflow = pe.Workflow(name=name)

    # Get subjects
    subjects = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                       name='get_subjects')
    subjects.iterables = ('subject_id', subject_ids)

    # Grab needed files
    inputfiles = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                   outfields=['dwi', 'fs_mask', 'b0_mask']),
                         name='get_input_files')
    inputfiles.inputs.base_directory = os.path.abspath(base_dir)
    inputfiles.inputs.template = ''
    inputfiles.inputs.field_template = dict(dwi='%s/NIFTI/DTI.nii.gz',
                                            fs_mask='%s/CMP/fs_output/HR__registered-TO-b0/fsmask_1mm.nii.gz',
                                            b0_mask='%s/NIFTI/DTI_first.nii.gz')
    inputfiles.inputs.template_args = dict(dwi=[['subject_id']],
                                           fs_mask=[['subject_id']],
                                           b0_mask=[['subject_id']])
    this_workflow.connect(subjects, 'subject_id', inputfiles, 'subject_id')

    # Prep masks - we need to convert voxel size for fs_mask, and
    # gunzip both images
    import nipype.interfaces.freesurfer.preprocess as fspre
    resample = pe.Node(interface=fspre.MRIConvert(), name='resample')
    resample.inputs.vox_size = tuple((3.,3.,3.))
    this_workflow.connect(inputfiles, 'fs_mask', resample, 'in_file')
    
    # DTK w/ fs mask
    dtk_fsmask = dtk_workflow(name='dtk_fsmask')
    this_workflow.connect(inputfiles, 'dwi', dtk_fsmask, 'inputnode.dwi')
    this_workflow.connect(resample, 'out_file', dtk_fsmask, 'inputnode.mask')
    # DTK w/ dwi mask
    # dtk_dwimask = dtk_workflow(name='dtk_dwimask')
    # this_workflow.connect(inputfiles, 'dwi', dtk_dwimask, 'inputnode.dwi')
    # this_workflow.connect(inputfiles, 'b0_mask', dtk_dwimask, 'inputnode.mask1')

    # Output
    rename1 = pe.Node(interface=util.Rename(
        format_string = "%(base_dir)s/%(subject_id)s/CMP/fibers/track_fs_%(seeds)s.trk"),
        name = 'rename1')
    rename1.inputs.base_dir = base_dir
    this_workflow.connect(subjects, 'subject_id', rename1, 'subject_id')
    this_workflow.connect(dtk_fsmask, 'outputnode.seeds', rename1, 'seeds')
    this_workflow.connect(dtk_fsmask, 'outputnode.track_file', rename1, 'in_file')
    
    # rename2 = pe.Node(interface=util.Rename(
    #     format_string = "%(base_dir)s/%(subject_id)s/CMP/fibers/track_b0_%(seeds)s.trk"),
    #     name = 'rename2')
    # rename2.inputs.base_dir = base_dir
    # this_workflow.connect(subjects, 'subject_id', rename2, 'subject_id')
    # this_workflow.connect(dtk_dwimask, 'outputnode.seeds', rename2, 'seeds')
    # this_workflow.connect(dtk_dwimask, 'outputnode.track_file', rename2, 'in_file')

    return this_workflow

    
def dtk_workflow(name='dtk_workflow'):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.diffusion_toolkit as dtk

    this_workflow = pe.Workflow(name=name)
    
    # Get input
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['dwi', 'mask', 'rseed_value'])
                        , name='inputnode')
    inputnode.inputs.rseed_value = 1

    # Diffusion recon
    recon = pe.Node(interface=dtk.dti.DTIRecon(), name='recon')
    recon.inputs.bvals = '/home/tnesland/Data/PipeComp/bvals.txt'
    recon.inputs.bvecs = '/home/tnesland/Data/PipeComp/bvecs.txt'
    recon.inputs.oblique_correction = True
    recon.inputs.output_type = 'nii.gz'
    this_workflow.connect(inputnode, 'dwi', recon, 'DWI')

    # Fiber tracking
    tracker = pe.Node(interface=dtk.dti.DTITracker(), name='tracker')
    tracker.inputs.input_type = 'nii.gz'
    this_workflow.connect(inputnode, 'mask', tracker, 'mask1_file')
    this_workflow.connect(inputnode, 'rseed_value', tracker, 'random_seed')
    this_workflow.connect(recon, 'tensor', tracker, 'tensor_file')

    # Spline filtering
    spline_filter = pe.Node(interface=dtk.postproc.SplineFilter(), name='spline_filter')
    spline_filter.inputs.step_length = 1.0
    this_workflow.connect(tracker, 'track_file', spline_filter, 'track_file')

    # Output
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['track_file', 'seeds'])
                         , name='outputnode')
    this_workflow.connect(inputnode, 'rseed_value', outputnode, 'seeds')
    this_workflow.connect(spline_filter, 'smoothed_track_file', outputnode, 'track_file')

    return this_workflow
