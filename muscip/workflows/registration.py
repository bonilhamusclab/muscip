"""Collection of registration workflows written for my purposes."""

import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.fsl.preprocess as fslpre
import nipype.interfaces.fsl.util as fslutil

def create_register_structural_to_diff(name='structural_to_diff'):
    """Co-register images using an affine transformation.

    Example
    -------

    >>> from tn_image_processing.workflows import registration
    >>> registration = registration.create_affine_registration()
    >>> registration.inputs.inputnode.structural = ['sub
    >>> registration.inputs.inputnode.dwi = ['DTI.nii.gz']

    Inputs::

        [Mandatory]
        inputnode.structural_list: list of hi-resolution structural
                                   images
        inputnode.dwi_list: list of 4dim diffusion volumes

    Outputs::

        outputnode.struct_to_b0: structural image registered to b0
                                 space
        outputnode.struct_to_b0_mat: affine matrix used to register
                                     structural image to b0 space

    """
    

    # Define the inputnode
    inputnode = pe.Node(interface=util.IdentityInterface(
        fields=["structural_list", "dwi_list"]), name="inputnode")

    # Get first volume from 4d diffusion image (we assume this is b0)
    extract_b0 = pe.MapNode(interface=fslutil.ExtractROI(
        roi_file='DTI_first.nii.gz', t_min=0, t_size=1),
        iterfield=['in_file'], name='extract_b0')

    # Align structural image to b0
    flirt = pe.MapNode(
        interface=fslpre.FLIRT(
            dof=6, cost='mutualinfo',
            usesqform=True,
            out_file="T1_to_b0.nii.gz",
            out_matrix_file="T1_to_b0.mat"),
        iterfield=['in_file', 'ref_file'],
        name="flirt")

    # Define this workflow
    structural_to_diff = pe.Workflow(name=name)

    # Connect components of this workflow
    structural_to_diff.connect([
        (inputnode, extract_b0, [("dwi", "in_file")]),
        (inputnode, flirt, [("structural", "in_file")]),
        (extract_b0, flirt, [("roi_file", "ref_file")]),
    ])

    # Define the outputnode
    outputnode = pe.Node(interface=util.IdentityInterface(
        fields=['struct_to_b0', 'struct_to_b0_mat']))

    # Connect the output
    structural_to_diff.connect([
        (flirt, outputnode, [('out_file', 'struct_to_b0'),
                             ('out_matrix_file', 'struct_to_b0_mat')]),
    ])

    return structural_to_diff
