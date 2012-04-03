import unittest
import connectome
import fibers
import nibabel

fib = fibers.read('test_data/fibers/CMTK_scale33.trk')
roi = nibabel.load('test_data/images/CMTK_roi.nii.gz')

class Test_TNConnectome(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_can_generate_connectome_from_fibers_and_rois(self):
        assert connectome.generate_connectome(fib, roi) is not None
