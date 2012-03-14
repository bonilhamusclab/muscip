import unittest
import fibers

FIBERS_DIR = 'test_data/fibers'
IMAGES_DIR = 'test_data/images'

# load trackvis files up front and share between all tests, otherwise,
# tests take too long to run
f_cmtk = fibers.read("%s/CMTK_scale33.trk" % FIBERS_DIR)  # CMTK
f_dtk = fibers.read("%s/DTK_RungeKutta.trk" % FIBERS_DIR) # DTK

class Test_TNImages(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_can_get_sum_of_inverse_lengths(self):
        # TODO: I need a better test for this, basically I need to
        # generate well-defined NIFTI images w/ more complexity than I
        # currently have and test the results
        pass
        
    def test_can_load_trackvis(self):
        global f_cmtk
        assert f_cmtk is not None

    def test_can_get_data(self):
        # test trackfile has a fiber count of 230814 according to
        # trackvis
        global f_cmtk
        assert len(f_cmtk.get_data()) == 230814

    def test_can_get_voxel_size(self):
        # test trackfile has a voxel size of [1,1,1]
        global f_cmtk
        assert (f_cmtk.get_voxel_size() == [1,1,1]).all()

    def test_can_get_dim(self):
        # test trackfile has dims of [208, 208, 120]
        global f_cmtk
        assert (f_cmtk.get_dim() == [208, 208, 120]).all()

    def test_can_get_voxel_spacing(self):
        # test trackfile has a voxel size of [2,2,2], so a given
        # vertex in voxel coordinates should be equal to half of the
        # same vertex in mm coordinates
        global f_dtk
        streamlines = f_dtk.get_data()
        mm_vertex = streamlines[0]['vertices_mm'][0]
        vox_vertex = streamlines[0]['vertices_vox'][0]
        assert ((mm_vertex / vox_vertex) == 2).all()

    def test_can_set_spacing(self):
        global f_cmtk
        f_cmtk.set_spacing('vox')
        assert f_cmtk.get_spacing() == 'vox'

    def test_can_map_vertices_to_roi(self):
        global f_cmtk
        f_cmtk.set_spacing('mm')
        import nibabel
        roi_img = nibabel.load("%s/CMTK_roi.nii.gz" % IMAGES_DIR)
        my_connectome = fibers.generate_connectome(f_cmtk, roi_img)
        import networkx
        expected_connectome = networkx.read_gpickle(
            "%s/CMTK_connectome_scale33.gpickle" % FIBERS_DIR)
        for i,j in my_connectome.edges_iter():
            assert my_connectome[i][j]['number_of_fibers'] == \
                   expected_connectome[i][j]['number_of_fibers']
            assert rms(my_connectome[i][j]['fiber_length_mean'],
                       expected_connectome[i][j]['fiber_length_mean']) < 10**-4

    def test_can_extract_scalars(self):
        global f_cmtk
        f_cmtk.set_spacing('mm')
        import nibabel
        roi_img = nibabel.load("%s/CMTK_roi.nii.gz" % IMAGES_DIR)
        my_connectome = fibers.generate_connectome(f_cmtk, roi_img)
        fa_img = nibabel.load("%s/CMTK_fa.nii.gz" % IMAGES_DIR)
        fibers.extract_scalars(f_cmtk, my_connectome, fa_img, 'fa', scale_factor=[2.,2.,2.])
        import networkx
        expected_connectome = networkx.read_gpickle(
            "%s/CMTK_connectome_scale33.gpickle" % FIBERS_DIR)
        for i,j in my_connectome.edges_iter():
            print i, j
            assert rms(my_connectome[i][j]['fa_mean'],
                       expected_connectome[i][j]['fa_mean']) < 10**-4

    def test_can_extract_hagmann_density(self):
        global f_cmtk
        f_cmtk.set_spacing('mm')
        import nibabel
        roi_img = nibabel.load("%s/CMTK_roi.nii.gz" % IMAGES_DIR)
        wm_img = nibabel.load("%s/CMTK_wm_1mm.nii.gz" % IMAGES_DIR)
        my_connectome = fibers.generate_connectome(f_cmtk, roi_img)
        fibers.extract_hagmann_density(my_connectome, roi_img, wm_img)
        for i,j in my_connectome.edges_iter():
            assert my_connectome[i][j]['hagmann_density'] > 0

    def test_can_get_cmat_for_key(self):
        global f_cmtk
        f_cmtk.set_spacing('mm')
        import nibabel
        roi_img = nibabel.load("%s/CMTK_roi.nii.gz" % IMAGES_DIR)
        my_connectome = fibers.generate_connectome(f_cmtk, roi_img)
        fiber_cmat = fibers.cmat_for_key(my_connectome, 'number_of_fibers')
        for i,j in my_connectome.edges_iter():
            expected_fibers = my_connectome[i][j]['number_of_fibers']
            assert fiber_cmat[i-1][j-1] == expected_fibers
            assert fiber_cmat[j-1][i-1] == expected_fibers


# helper functions
def rms(value1, value2):
    from math import sqrt
    return sqrt((value1-value2 )**2)
        
        
        

