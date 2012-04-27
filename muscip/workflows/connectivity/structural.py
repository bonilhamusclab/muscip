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

    # Create our workflow so we can connect as we go ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this_workflow = pe.Workflow(name=name)
    # Grab data for all subject ids from base directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datasource = create_datasource(base_dir, subject_ids)
    # Registration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    registration = create_registration()
    this_workflow.connect(datasource, 'outputnode.dwi_image', registration, 'inputnode.dwi_image')
    this_workflow.connect(datasource, 'outputnode.struct_image', registration, 'inputnode.structural_image')
    # Segmentation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    segmentation = create_segmentation(base_dir)
    this_workflow.connect(datasource, 'outputnode.struct_image', segmentation, 'inputnode.structural_image')
    this_workflow.connect(datasource, 'outputnode.subject_id', segmentation, 'inputnode.subject_id')
    # Atlas Creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fs_atlas = create_fs_atlas()
    this_workflow.connect(datasource, ('outputnode.subject_id', real_freesurfer_dir, base_dir)
                          , fs_atlas, 'inputnode.freesurfer_dir')
    # Apply Registration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    apply_registration = create_apply_registration()
    this_workflow.connect(datasource, 'outputnode.dwi_image', apply_registration, 'inputnode.reference')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_mat', apply_registration, 'inputnode.xform')
    this_workflow.connect(fs_atlas, 'outputnode.roi_img', apply_registration, 'inputnode.roi_file')
    this_workflow.connect(fs_atlas, 'outputnode.wm_img', apply_registration, 'inputnode.wm_mask')
    # Tractography ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tractography = create_tractography()
    this_workflow.connect(datasource, 'outputnode.dwi_image', tractography, 'inputnode.dwi')
    this_workflow.connect(datasource, 'outputnode.bvals', tractography, 'inputnode.bvals')
    this_workflow.connect(datasource, 'outputnode.bvecs', tractography, 'inputnode.bvecs')
    this_workflow.connect(datasource, ('outputnode.bvals', number_of_b0_for_bvals),
                          tractography, 'inputnode.number_b0')
    # Connectome Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    connectome_mapper = create_connectome()
    this_workflow.connect(tractography, 'outputnode.spline_filtered', connectome_mapper, 'inputnode.track_file')
    this_workflow.connect(apply_registration, 'outputnode.roi_file', connectome_mapper, 'inputnode.roi_file')
    this_workflow.connect(tractography, 'outputnode.mask_file', connectome_mapper, 'inputnode.wm_file')
    # Write out data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datawriter = create_datawriter(base_dir)
    this_workflow.connect(datasource, 'outputnode.subject_id', datawriter, 'inputnode.subject_id')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_image', datawriter, 'inputnode.t1_to_b0_img')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_mat', datawriter, 'inputnode.t1_to_b0_mat')
    this_workflow.connect(fs_atlas, 'outputnode.roi_img', datawriter, 'inputnode.roi')
    this_workflow.connect(fs_atlas, 'outputnode.wm_img', datawriter, 'inputnode.fs_mask')
    this_workflow.connect(apply_registration, 'outputnode.roi_file', datawriter, 'inputnode.roi_to_b0')
    this_workflow.connect(apply_registration, 'outputnode.wm_mask', datawriter, 'inputnode.fs_mask_to_b0')
    this_workflow.connect(tractography, 'outputnode.mask_file', datawriter, 'inputnode.fa_mask')
    this_workflow.connect(tractography, 'outputnode.streamline', datawriter, 'inputnode.streamline')
    this_workflow.connect(tractography, 'outputnode.spline_filtered', datawriter, 'inputnode.spline_filtered')
    this_workflow.connect(connectome_mapper, 'outputnode.network_file', datawriter, 'inputnode.network_file')
    # Return this workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return this_workflow

def create_post_fs_dti_workflow(base_dir, subject_ids, name='post_fs_workflow'):
    """Generate connectomes given that freesurfer recon-all has
    already been run and we have access to folders

    Example
    -------

    >>> from muscip.workflows.connectivity import structural
    >>> workflow = structural.create_post_fs_dti_workflow(base_dir, subject_ids)
    >>> workflow.run()

    Inputs::

      base_dir: directory containing all subjects
      subject_ids: list of strings corresponding to subject ids
    """
    import nipype.pipeline.engine as pe
    # Create our workflow so we can connect as we go ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this_workflow = pe.Workflow(name=name)
    # Grab data for all subject ids from base directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datasource = create_datasource(base_dir, subject_ids)
    # Registration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    registration = create_registration()
    this_workflow.connect(datasource, 'outputnode.dwi_image', registration, 'inputnode.dwi_image')
    this_workflow.connect(datasource, 'outputnode.struct_image', registration, 'inputnode.structural_image')
    # Atlas Creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    fs_atlas = create_fs_atlas()
    this_workflow.connect(datasource, ('outputnode.subject_id', real_freesurfer_dir, base_dir)
                          , fs_atlas, 'inputnode.freesurfer_dir')
    # Apply Registration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    apply_registration = create_apply_registration()
    this_workflow.connect(datasource, 'outputnode.dwi_image', apply_registration, 'inputnode.reference')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_mat', apply_registration, 'inputnode.xform')
    this_workflow.connect(fs_atlas, 'outputnode.roi_img', apply_registration, 'inputnode.roi_file')
    this_workflow.connect(fs_atlas, 'outputnode.wm_img', apply_registration, 'inputnode.wm_mask')
    # Tractography ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    tractography = create_tractography()
    this_workflow.connect(datasource, 'outputnode.dwi_image', tractography, 'inputnode.dwi')
    this_workflow.connect(datasource, 'outputnode.bvals', tractography, 'inputnode.bvals')
    this_workflow.connect(datasource, 'outputnode.bvecs', tractography, 'inputnode.bvecs')
    this_workflow.connect(datasource, ('outputnode.bvals', number_of_b0_for_bvals),
                          tractography, 'inputnode.number_b0')
    # Connectome Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    connectome_mapper = create_connectome()
    this_workflow.connect(fs_atlas, 'outputnode.node_info', connectome_mapper, 'inputnode.node_info')
    this_workflow.connect(tractography, 'outputnode.spline_filtered', connectome_mapper, 'inputnode.track_file')
    this_workflow.connect(apply_registration, 'outputnode.roi_file', connectome_mapper, 'inputnode.roi_file')
    this_workflow.connect(tractography, 'outputnode.mask_file', connectome_mapper, 'inputnode.wm_file')
    # Write out data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    datawriter = create_datawriter(base_dir)
    this_workflow.connect(datasource, 'outputnode.subject_id', datawriter, 'inputnode.subject_id')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_image', datawriter, 'inputnode.t1_to_b0_img')
    this_workflow.connect(registration, 'outputnode.struct_to_dwi_mat', datawriter, 'inputnode.t1_to_b0_mat')
    this_workflow.connect(fs_atlas, 'outputnode.roi_img', datawriter, 'inputnode.roi')
    this_workflow.connect(fs_atlas, 'outputnode.wm_img', datawriter, 'inputnode.fs_mask')
    this_workflow.connect(apply_registration, 'outputnode.roi_file', datawriter, 'inputnode.roi_to_b0')
    this_workflow.connect(apply_registration, 'outputnode.wm_mask', datawriter, 'inputnode.fs_mask_to_b0')
    this_workflow.connect(tractography, 'outputnode.mask_file', datawriter, 'inputnode.fa_mask')
    this_workflow.connect(tractography, 'outputnode.streamline', datawriter, 'inputnode.streamline')
    this_workflow.connect(tractography, 'outputnode.spline_filtered', datawriter, 'inputnode.spline_filtered')
    this_workflow.connect(connectome_mapper, 'outputnode.network_file', datawriter, 'inputnode.network_file')
    # Return this workflow ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    return this_workflow
    

def create_datasource(base_dir, subject_ids, name='datasource'):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    # Create our workflow so we can connect as we go ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this_workflow = pe.Workflow(name=name)
    # Grab data for all subject ids from base directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['base_dir',
                                                                 'subject_id']),
                         name='inputnode')
    inputnode.iterables = ('subject_id', subject_ids)
    import nipype.interfaces.io as nio
    import os.path
    datagrab = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                                 outfields=['struct_image',
                                                            'dwi_image',
                                                            'bvals',
                                                            'bvecs']), name='datagrab')
    datagrab.inputs.base_directory = os.path.abspath(base_dir)    
    datagrab.inputs.template = ''
    datagrab.inputs.field_template = dict(struct_image='%s/NIFTI/T1.nii.gz',
                                          dwi_image='%s/NIFTI/DTI.nii.gz',
                                          bvals='%s/bvals.txt',
                                          bvecs='%s/bvecs.txt')
    datagrab.inputs.template_args = dict(struct_image=[['subject_id']],
                                         dwi_image=[['subject_id']],
                                         bvals=[['subject_id']],
                                         bvecs=[['subject_id']])
    this_workflow.connect(inputnode, 'subject_id', datagrab, 'subject_id')
    # Create output node
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['base_dir',
                                                                  'subject_id',
                                                                  'struct_image',
                                                                  'dwi_image',
                                                                  'bvals',
                                                                  'bvecs']), name='outputnode')
    this_workflow.connect(inputnode, 'base_dir', outputnode, 'base_dir')    
    this_workflow.connect(inputnode, 'subject_id', outputnode, 'subject_id')
    this_workflow.connect(datagrab, 'struct_image', outputnode, 'struct_image')
    this_workflow.connect(datagrab, 'dwi_image', outputnode, 'dwi_image')
    this_workflow.connect(datagrab, 'bvals', outputnode, 'bvals')
    this_workflow.connect(datagrab, 'bvecs', outputnode, 'bvecs')
    # returnt this workflow
    return this_workflow

def create_datawriter(base_dir, name='datawriter'):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.io as nio
    import os.path
    # create this workflow
    this_workflow = pe.Workflow(name=name)
    # inputnode
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['subject_id',
                                                                 't1_to_b0_img',
                                                                 't1_to_b0_mat',
                                                                 'roi',
                                                                 'fs_mask',
                                                                 'roi_to_b0',
                                                                 'fs_mask_to_b0',
                                                                 'fa_mask',
                                                                 'streamline',
                                                                 'spline_filtered',
                                                                 'network_file'
                                                             ]), name='inputnode')
    # datasink
    datasink = pe.Node(nio.DataSink(parameterization=False), name='datasink')
    datasink.inputs.base_directory = os.path.abspath(base_dir)
    datasink.inputs.substitutions = [('DTI_roi.nii.gz', 'DTI_first.nii.gz'),
                                     ('DTI_roi_out.nii.gz', 'Diffusion_b0_resampled.nii.gz'),
                                     ('T1_flirt.nii.gz', 'T1-TO-b0.nii.gz'),
                                     ('T1_flirt.mat', 'T1-TO-b0.mat'),
                                     ('roi.nii.gz', 'ROI.nii.gz')]
    # connect the pipes
    this_workflow.connect(inputnode, 'subject_id', datasink, 'container')
    this_workflow.connect(inputnode, 't1_to_b0_img', datasink, 'NIFTI.@t1_to_b0_img')
    this_workflow.connect(inputnode, 't1_to_b0_mat', datasink, 'NIFTI.transformations')
    this_workflow.connect(inputnode, 'roi', datasink, 'CMP.fs_output.HR.@roi')
    this_workflow.connect(inputnode, 'fs_mask', datasink, 'CMP.fs_output.HR.@fs_mask')            
    this_workflow.connect(inputnode, 'roi_to_b0', datasink, 'CMP.fs_output.HR__registered-TO-b0.@roi')
    this_workflow.connect(inputnode, 'fs_mask', datasink, 'CMP.fs_output.HR__registered-TO-b0.@fs_mask')
    this_workflow.connect(inputnode, 'fa_mask', datasink, 'CMP.fs_output.HR__registered-TO-b0.@fa_mask')
    this_workflow.connect(inputnode, 'streamline', datasink, 'CMP.fibers.@streamline')
    this_workflow.connect(inputnode, 'spline_filtered', datasink, 'CMP.fibers.@spline_filtered')
    this_workflow.connect(inputnode, 'network_file', datasink, 'CMP.fibers.matrices')
    # return this workflow
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
    extract_b0 = pe.Node(interface=fslutil.ExtractROI(t_min=0, t_size=1), name='extract_b0')
    # TODO: debug runtime error caused when specifying name of output roi file
    # extract_b0.inputs.roi_file = 'DTI_first.nii.gz'
    #
    # Resample b0 to voxel size of 1mm
    import nipype.interfaces.freesurfer.preprocess as fspre
    resample_b0 = pe.Node(interface=fspre.Resample(),
                          name='resample_b0')
    resample_b0.inputs.voxel_size = tuple((1.,1.,1.))
    # resample_b0.inputs.resampled_file = 'Diffusion_b0_resampled.nii.gz'
    # Align structural image to b0
    flirt = pe.Node(interface=fslpre.FLIRT( dof=6, cost='mutualinfo',
                                            uses_qform=True), name="flirt")
    # flirt.inputs.out_file = 'T1-TO-b0.nii.gz'
    # flirt.inputs.out_matrix_file = 'T1-TO-b0.mat'
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
    this_workflow.connect(resample_b0, "resampled_file", flirt, "reference")
    this_workflow.connect(extract_b0, "roi_file", outputnode, "dwi_first")
    this_workflow.connect(resample_b0, "resampled_file", outputnode, "dwi_resampled")
    this_workflow.connect(flirt, "out_file", outputnode, "struct_to_dwi_image")
    this_workflow.connect(flirt, "out_matrix_file", outputnode, "struct_to_dwi_mat")    
    # Return workflow
    return this_workflow

def create_tractography(name='tractography'):
    """Create a workflow designed to reconstruct diffusion image so
    that it is ready for fiber tracking.

    Example
    -------
    >>> import tn_image_processing.workflows as pipes
    >>> reconstruction = pipes.create_reconstruction()
    >>> reconstruction.inputnode.diff_volume = ['path/to/dwi/volume']
    >>> reconstruction.inputnode.gradient_table = ['path/to/gradient/table']
    >>> reconstruction.inputnode.number_b0 = 1 # number of b-zeros
    >>> reconstruction.run()

    Inputs::

      [mandatory]
      inputnode.dwi
      inputnode.bvals
      inputnode.bvecs
      inputnode.number_b0

      [optional]
      inputnode.mask

    Outputs::

      outputnode.mask_file
      outputnode.track_file
      outputnode.smoothed_track_file

    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.diffusion_toolkit as dtk
    # create this workflow
    this_workflow = pe.Workflow(name=name)
    # Get input
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['dwi', 'bvals', 'bvecs', 'number_b0'])
                        , name='inputnode')
    # Eddy Corrections
    import nipype.workflows.dmri.fsl.dti
    eddy_correct = nipype.workflows.dmri.fsl.dti.create_eddy_correct_pipeline("eddy_correct")
    eddy_correct.inputs.inputnode.ref_num = 0
    this_workflow.connect(inputnode, 'dwi', eddy_correct, 'inputnode.in_file')
    # Diffusion recon
    recon = pe.Node(interface=dtk.dti.DTIRecon(), name='recon')
    recon.inputs.oblique_correction = True
    recon.inputs.output_type = 'nii.gz'
    this_workflow.connect(eddy_correct, 'outputnode.eddy_corrected', recon, 'DWI')
    this_workflow.connect(inputnode, 'bvals', recon, 'bvals')
    this_workflow.connect(inputnode, 'bvecs', recon, 'bvecs')
    # Fiber
    tracker = pe.Node(interface=dtk.dti.DTITracker(), name='tracker')
    tracker.inputs.input_type = 'nii.gz'
    tracker.inputs.random_seed = 32
    tracker.inputs.mask1_threshold = 0.2
    tracker.inputs.output_file = 'streamline.trk'
    tracker.inputs.output_mask = 'wm_mask.nii.gz'
    this_workflow.connect(recon, 'FA', tracker, 'mask1_file')
    this_workflow.connect(recon, 'tensor', tracker, 'tensor_file')
    # Spline filtering
    spline_filter = pe.Node(interface=dtk.postproc.SplineFilter(), name='spline_filter')
    spline_filter.inputs.step_length = 1.0
    this_workflow.connect(tracker, 'track_file', spline_filter, 'track_file')
    # Output
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['mask_file',
                                                                  'streamline',
                                                                  'spline_filtered'])
                         , name='outputnode')
    this_workflow.connect(tracker, 'mask_file', outputnode, 'mask_file')
    this_workflow.connect(tracker, 'track_file', outputnode, 'streamline')
    this_workflow.connect(spline_filter, 'smoothed_track_file', outputnode, 'spline_filtered')
    # return workflow
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
    this_workflow.connect(inputnode, ('subject_id', freesurfer_dir, _base_dir), recon_all, 'subjects_dir')
    this_workflow.connect(inputnode, 'structural_image', recon_all, 'T1_files')
    # Output node
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['subjects_dir',
                                                                  'subject_id']),
                         name='outputnode')
    this_workflow.connect(recon_all, 'subjects_dir', outputnode, 'subjects_dir')
    this_workflow.connect(recon_all, 'subject_id', outputnode, 'subject_id')

    # Return workflow
    return this_workflow

def create_fs_atlas(name='fsatlas'):
    """Generate Freesurfer atlas given a fs directory that has been
    processed with recon-all

    Inputs::

      [mandatory]
      inputnode.freesurfer_dir: directory successfully processed with recon-all
      
    Outputs::

      node_info: node information associated with freesurfer atlas
      roi_img: roi image extracted from freesurfer
      wm_img: wm image extracted from freesurfer
    """
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import muscip.interfaces.atlas as atlas
    # Create this workflow
    this_workflow = pe.Workflow(name=name)
    # Setup input node
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['freesurfer_dir']),
                        name='inputnode')
    # Setup atlas extract node
    atlasExtract = pe.Node(interface=atlas.FreesurferAtlas(), name='atlasExtract')
    this_workflow.connect(inputnode, 'freesurfer_dir', atlasExtract, 'freesurfer_dir')
    # Setup output node
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['node_info',
                                                                  'roi_img',
                                                                  'wm_img']), name='outputnode')
    this_workflow.connect(atlasExtract, 'node_info', outputnode, 'node_info')
    this_workflow.connect(atlasExtract, 'roi_img', outputnode, 'roi_img')
    this_workflow.connect(atlasExtract, 'wm_img', outputnode, 'wm_img')
    # Return this workflow
    return this_workflow

def create_apply_registration(name="registerToBZero"):
    """Apply the given transformation to the list of input files"""
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.preprocess as fslpre
    # Create this workfow
    this_workflow = pe.Workflow(name=name)
    # Create inputnode
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['reference',
                                                                 'roi_file',
                                                                 'wm_mask',
                                                                 'xform']), name='inputnode')
    # Co-register ROI - b0
    register_ROI = pe.Node(interface=fslpre.ApplyXfm(), name='coregisterROI')
    register_ROI.inputs.output_type = 'NIFTI_GZ'
    # register_ROI.inputs.out_file = 'ROI_to_b0.nii.gz'
    register_ROI.inputs.apply_xfm = True
    register_ROI.inputs.interp = 'nearestneighbour'
    this_workflow.connect(inputnode, 'roi_file', register_ROI, 'in_file')
    this_workflow.connect(inputnode, 'reference', register_ROI, 'reference')
    this_workflow.connect(inputnode, 'xform', register_ROI, 'in_matrix_file')
    # Co-register wm_mask - b0
    register_WM = pe.Node(interface=fslpre.ApplyXfm(), name='coregisterWM')
    register_WM.inputs.output_type = 'NIFTI_GZ'
    # register_WM.inputs.out_file = 'fsmask_to_b0.nii.gz'
    register_WM.inputs.apply_xfm = True
    register_WM.inputs.interp = 'nearestneighbour'
    this_workflow.connect(inputnode, 'wm_mask', register_WM, 'in_file')
    this_workflow.connect(inputnode, 'reference', register_WM, 'reference')
    this_workflow.connect(inputnode, 'xform', register_WM, 'in_matrix_file')
    # Create output node
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['roi_file', 'wm_mask']),
                         name='outputnode')
    this_workflow.connect(register_ROI, 'out_file', outputnode, 'roi_file')
    this_workflow.connect(register_WM, 'out_file', outputnode, 'wm_mask')
    # Return this workflow
    return this_workflow

def create_connectome(name='connectome'):
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    # Create this workflow
    this_workflow = pe.Workflow(name=name)
    # Inputnode
    inputnode = pe.Node(interface=util.IdentityInterface(fields=['node_info',
                                                                 'roi_file',
                                                                 'track_file',
                                                                 'wm_file']), name='inputnode')
    # Connectome Generator
    import muscip.interfaces.connectome as ctome
    mapper = pe.Node(interface=ctome.ConnectomeGenerator(), name='mapper')
    mapper.inputs.network_file = 'connectome'
    this_workflow.connect(inputnode, 'node_info', mapper, 'node_info_file')
    this_workflow.connect(inputnode, 'roi_file', mapper, 'roi_file')
    this_workflow.connect(inputnode, 'track_file', mapper, 'track_file')
    this_workflow.connect(inputnode, 'wm_file', mapper, 'wm_file')
    # Outputnode
    outputnode = pe.Node(interface=util.IdentityInterface(fields=['network_file']), name='outputnode')
    this_workflow.connect(mapper, 'network_file', outputnode, 'network_file')
    # Return this workflow
    return this_workflow
    
def freesurfer_dir(subject_id, base_dir):
    """Return a freesurfer dir for given params"""
    import os.path
    return os.path.join(os.path.abspath(base_dir),subject_id)

def real_freesurfer_dir(subject_id, base_dir):
    """Return a freesurfer dir for given params"""
    import os.path
    return os.path.join(os.path.abspath(base_dir),subject_id,'FREESURFER')
    
def number_of_b0_for_bvals(bvals):
    """Return the number of b0's in the given b-value file."""
    import numpy
    # read in bvals
    try:
        read_bvals = numpy.loadtxt(bvals)
    except Exception as e:
        print e
        return None
    return len(numpy.where(read_bvals==0))
