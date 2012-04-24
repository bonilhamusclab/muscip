import os
from nipype.interfaces.base import (BaseInterface, Directory, File,
                                    TraitedSpec)
from nipype import logging
iflogger = logging.getLogger('interface')

class FreesurferAtlasInputSpec(TraitedSpec):
    freesurfer_dir = Directory(exists=True, mandatory=True,
                               desc="Freesurfer directory after recon-all")
    node_info = File('node_info.graphml', usedefault=True, desc='Output name for node info file')
    roi_img = File('roi.nii.gz', usedefault=True, desc='Output name for roi image')
    wm_img = File('wm.nii.gz', usedefault=True, desc='Output name for wm image')    

class FreesurferAtlasOutputSpec(TraitedSpec):
    node_info = File(desc="Node information stored in graphML file")
    roi_img = File(desc="Generated ROI image for the given Freesurfer Dir")
    wm_img = File(desc="Generated WM image for the given Freesurfer Dir")

class FreesurferAtlas(BaseInterface):
    """Generate Freesurfer Atlas and related files given a processed
    Freesurfer directory *1

    Example
    -------
    >>> import muscip.interfaces.atlas as matlas
    >>> fsatlas = matlas.FreesurferAtlas()
    >>> fsatlas.inputs = 'path/to/processed/freesurfer/dir'
    >>> fsatlas.run() # doctest: +SKIP

    """

    input_spec = FreesurferAtlasInputSpec
    output_spec = FreesurferAtlasOutputSpec

    def _run_interface(self, runtime):
        from ..atlas import freesurfer as fs
        import networkx as nx
        import nibabel as nib
        iflogger.info('Extracting Freesurfer Atlas...')
        my_fs_atlas = fs.load(self.inputs.freesurfer_dir)
        iflogger.info('Writing Freesurfer Atlas Files...')
        nib.save(my_fs_atlas.roi_img, self.inputs.roi_img)
        nib.save(my_fs_atlas.wm_mask, self.inputs.wm_img)
        nx.write_graphml(my_fs_atlas.node_info, self.inputs.node_info)
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['node_info'] = os.path.abspath(self.inputs.node_info)
        outputs['roi_img'] = os.path.abspath(self.inputs.roi_img)
        outputs['wm_img'] = os.path.abspath(self.inputs.wm_img)
        return outputs
