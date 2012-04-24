import unittest
import muscip.atlas.freesurfer as fs

class Test_Freesurfer(unittest.TestCase):

    def setUp(self):
        import os.path as op
        self.freesurfer_dir = op.join(op.dirname(fs.__file__),'tests/data/freesurfer_dir')
        self.cmtk_dir = op.join(op.dirname(fs.__file__),'tests/data/cmtk_dir')
        self.cmtk_rois_path = op.join(self.cmtk_dir, 'freesurfer_roi.nii.gz')
        self.cmtk_wm_mask_path = op.join(self.cmtk_dir, 'freesurfer_wm.nii.gz')

    def tearDown(self):
        pass

    def test_freesurfer_rois_eq_to_cmtk_freesurfer_rois(self):
        my_fs = fs.load(self.freesurfer_dir)
        import nibabel
        cmtk_fs_rois = nibabel.load(self.cmtk_rois_path)
        assert my_fs.roi_img.get_data().all() == cmtk_fs_rois.get_data().all()

    def test_freesurfer_wm_mask_eq_to_cmtk_freesurfer_wm_mask(self):
        my_fs = fs.load(self.freesurfer_dir)
        import nibabel
        cmtk_wm_mask = nibabel.load(self.cmtk_wm_mask_path)
        assert my_fs.wm_mask.get_data().all() == cmtk_wm_mask.get_data().all()

    def test_freesurfer_has_node_info(self):
        my_fs = fs.load(self.freesurfer_dir)
        assert my_fs.node_info.node[str(82)]['name'] == 'Left-Amygdala'
