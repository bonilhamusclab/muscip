"""Collection of workflows relating to structural connectivity."""

import os

def create_default_dti_workflow(base_dir,
                                subject_ids,
                                name='default_dti_workflow'):
    """Process DTI using our default workflow.

    Example
    -------

    >>> from tn_image_processing.workflows.connectivity import structural
    >>> dti_connectivity = structural.create_default_dti_workflow()
    >>> dti_connectivity.inputnode.subject_dir = ['path/to/subject/dir']
    >>> dti_connectivity.inputnode.subject_list = ['s001','s002','s003']

    Inputs::

        [Mandatory]
        base_dir: directory containing all subjects
        subject_ids: list of subject strings

    Outputs::

        None

    """

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    # Create our workflow so we can connect as we go ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this_workflow = pe.Workflow(name=name)
    
    # Grab data for all subject ids from base directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['subject_id']),
                         name='inputnode')
    inputnode.iterables = ('subject_id', subject_ids)
    
    import nipype.interfaces.io as nio
    import os.path
    datagrab = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                 outfields=['struct_image',
                                                            'dwi_image']),
                                                            name='datagrab')
    datagrab.inputs.base_directory = os.path.abspath(base_dir)    
    datagrab.inputs.template = ''
    datagrab.inputs.field_template = dict(struct_image='%s/NIFTI/T1.nii.gz',
                                            dwi_image='%s/NIFTI/DTI.nii.gz')
    datagrab.inputs.template_args = dict(struct_image=[['subject_id']],
                                           dwi_image=[['subject_id']])
    
    this_workflow.connect(inputnode, 'subject_id', datagrab, 'subject_id')
    
    # Registration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    registration = create_registration()
    this_workflow.connect(datagrab, 'dwi_image', registration, 'inputnode.dwi_image')
    this_workflow.connect(datagrab, 'struct_image', registration, 'inputnode.structural_image')

    # Segmentation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    segmentation = create_segmentation(base_dir)
    # this_workflow.connect(
    #     [datagrab, segmentation,
    #      [('struct_image', 'inputnode.structural_image'),
    #       ('subject_id', 'inputnode.subject_id')],
    # ])
    this_workflow.connect(datagrab, 'struct_image', segmentation, 'inputnode.structural_image')
    this_workflow.connect(inputnode, 'subject_id', segmentation, 'inputnode.subject_id')
    
    # Write out data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datasink = pe.Node(nio.DataSink(parameterization=False), name='nifti-sinker')
    datasink.inputs.base_directory = os.path.abspath(base_dir)
    datasink.inputs.substitutions = [('DTI_roi.nii.gz', 'DTI_first.nii.gz'),
                                     ('DTI_roi_out.nii.gz', 'Diffusion_b0_resampled.nii.gz'),
                                     ('T1_flirt.nii.gz', 'T1-TO-b0.nii.gz'),
                                     ('T1_flirt.mat', 'transformations/T1-TO-b0.mat')]
    niftiMerge = pe.Node(interface=util.Merge(4), name='nifti-file-merge')
    this_workflow.connect(registration, 'outputnode.dwi_first', niftiMerge, 'in1')
    this_workflow.connect(registration, 'outputnode.dwi_resampled', niftiMerge, 'in2')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_image', niftiMerge, 'in3')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_mat', niftiMerge, 'in4')
    this_workflow.connect(inputnode, 'subject_id', datasink, 'container')
    this_workflow.connect(niftiMerge, 'out', datasink, 'NIFTI')
    
    # Return this workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return this_workflow

def create_registration(name='registration'):
    """Create a workflow designed to register structural image to dwi.

    Example
    -------
    >>> import tn_image_processing.workflows as pipes
    >>> registration = pipes.create_registration()
    >>> registration.inputnode.dwi_image = ['path/to/dwi/image']
    >>> registration.inputnode.struct_image = ['path/to/struct/image']
    >>> registration.run()

    Inputs::

      [mandatory]
      inputnode.dwi_image
      inputnode.structural_image
      inputnode.working_dir

    Outputs::

      outputnode.dwi_first
      outputnode.dwi_resampled
      outputnode.struct_to_dwi_image
      outputnode.struct_to_dwi_mat

    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.utils as fslutil
    import nipype.interfaces.fsl.preprocess as fslpre    
    
    # Setup input node
    inputnode = pe.Node(interface=util.IdentityInterface(
        fields=['dwi_image','structural_image']),
                        name='inputnode')
    
    # Get first volume from 4d diffusion image (we assume this is b0)
    extract_b0 = pe.Node(interface=fslutil.ExtractROI(t_min=0,
                                                      t_size=1), name='extract_b0')

    # Resample b0 to voxel size of 1mm
    import nipype.interfaces.freesurfer.preprocess as fspre
    resample_b0 = pe.Node(interface=fspre.MRIConvert(),
                          name='resample_b0')
    resample_b0.inputs.vox_size = tuple((1.,1.,1.))

    # Align structural image to b0
    flirt = pe.Node(interface=fslpre.FLIRT( dof=6, cost='mutualinfo',
                                            uses_qform=True),
                    name="flirt")

    # Setup output node
    outputnode = pe.Node(interface=util.IdentityInterface(
        fields=['dwi_first', 'dwi_resampled', 'struct_to_dwi_image',
                'struct_to_dwi_mat']), name='outputnode')

    # Create workflow
    this_workflow = pe.Workflow(name=name)
    # Connect workflow
    this_workflow.connect(inputnode, "dwi_image", extract_b0, "in_file")
    this_workflow.connect(inputnode, "structural_image", flirt, "in_file")
    this_workflow.connect(extract_b0, "roi_file", resample_b0, "in_file")
    this_workflow.connect(resample_b0, "out_file", flirt, "reference")
    this_workflow.connect(extract_b0, "roi_file", outputnode, "dwi_first")
    this_workflow.connect(resample_b0, "out_file", outputnode, "dwi_resampled")
    this_workflow.connect(flirt, "out_file", outputnode, "struct_to_dwi_image")
    this_workflow.connect(flirt, "out_matrix_file", outputnode, "struct_to_dwi_mat")    
    # Return workflow
    return this_workflow

def create_segmentation(base_dir, name='segmentation'):
    """Create a workflow designed to segment structural image.

    Example
    -------
    >>> import tn_image_processing.workflows as pipes
    >>> segmentation = pipes.create_segmentation()
    >>> segmentation.inputnode.dwi_image = ['path/to/dwi/image']
    >>> segmentation.inputnode.struct_image = ['path/to/struct/image']
    >>> segmentation.run()

    Inputs::

      [mandatory]
      inputnode.structural_image
      inputnode.subject_id

    Outputs::

    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.freesurfer.preprocess as fs

    # Get base dir
    import os
    _base_dir = os.path.abspath(base_dir)
    
    # Create this workflow
    this_workflow = pe.Workflow(name=name)
    
    # Setup input node
    inputnode = pe.Node(interface=util.IdentityInterface(
        fields=['structural_image', 'subject_id']),
                        name='inputnode')
    
    # Freesurfer recon_all
    recon_all = pe.Node(interface=fs.ReconAll(directive="all",
                                              subject_id='FREESURFER')
                        , name='recon_all')

    # Connect workflow
    this_workflow.connect(inputnode, ('subject_id', freesurfer_dir, _base_dir), recon_all, 'subjects_dir')
    this_workflow.connect(inputnode, 'structural_image', recon_all, 'T1_files')
    
    # Return workflow
    return this_workflow

def freesurfer_dir(subject_id, base_dir):
    """Return a freesurfer dir for given params"""
    import os.path
    return os.path.join(os.path.abspath(base_dir),subject_id)
