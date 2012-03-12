import unittest
import fibers

FIBERS_DIR = 'test_data/fibers'
IMAGES_DIR = 'test_data/images/3dim_images'


f = fibers.read("%s/dti.trk" % FIBERS_DIR) # load once to avoid excess
                                           # run time


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
        global f
        assert f is not None

    def test_can_get_data(self):
        # test trackfile has a fiber count of 157434 according to
        # trackvis
        global f
        assert len(f.get_data()) == 157434

    def test_can_get_voxel_size(self):
        # test trackfile has a voxel size of [2,2,2]
        global f
        assert (f.get_voxel_size() == [2,2,2]).all()

    def test_can_get_dim(self):
        # test trackfile has dims of [104, 104, 60]
        global f
        assert (f.get_dim() == [104, 104, 60]).all()

    def test_can_get_voxel_spacing(self):
        # test trackfile has a voxel size of [2,2,2], so a given
        # vertex in voxel coordinates should be equal to half of the
        # same vertex in mm coordinates
        global f
        streamlines = f.get_data()
        mm_vertex = streamlines[0]['vertices_mm'][0]
        vox_vertex = streamlines[0]['vertices_vox'][0]
        assert ((mm_vertex / vox_vertex) == 2).all()
