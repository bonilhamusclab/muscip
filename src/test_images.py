import unittest
import images

SUPPORT_DIR = 'test_data/images'

class Test_TNImages(unittest.TestCase):

    def setUp(self):
        self._i = images.TNImages()

    def tearDown(self):
        pass

    def test_can_instantiate_class(self):
        assert images.TNImages() is not None
        
    def test_can_load_images_from_directory(self):
        # should load images, and only images, from directory we have
        # placed four '.nii.gz' images in a directory in
        # tests/data/images, we should be able to load
        self._i.load_from_dir(SUPPORT_DIR)
        assert self._i.get_images() is not None
        assert len(self._i.get_images()) == 3

    def test_does_not_choke_on_non_image_file(self):
        error = False
        try:
            self._i.load_from_dir(SUPPORT_DIR)
        except:
            error = True
        assert not error
        assert len(self._i.get_images()) == 3

    def test_images_are_tn_images(self):
        self._i.load_from_dir(SUPPORT_DIR)
        test_img = self._i.get_images()[0]
        assert type(test_img) is images.TNImage

class Test_TNImage(unittest.TestCase):

    def setUp(self):
        self._i = images.TNImage()

    def tearDown(self):
        pass

    def test_can_instantiate_class(self):
        assert images.TNImage() is not None

    def test_can_instantiate_class_from_filename(self):
        from os.path import join
        assert images.TNImage(join(SUPPORT_DIR,'noise.nii.gz')) is not None

    def test_can_get_dz(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert (my_image.dz() == 0).all()

    def test_dz_is_accurate(self):
        from os.path import join
        test_img = join(SUPPORT_DIR,'noise.nii.gz')
        my_obj = images.TNImage(filename=test_img)
        mean_z = my_obj.z_mean()
        dz = my_obj.dz()
        assert (mean_z[1,:] - mean_z[0,:] == dz[1,:]).all()
        assert (mean_z[0,:] - mean_z[-1,:] == dz[0,:]).all()

    def test_can_get_z_mean(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert my_image.z_mean() is not None

    def test_z_mean_is_accurate(self):
        from nibabel import load
        from os.path import join
        import numpy as np
        test_img = join(SUPPORT_DIR,'noise.nii.gz')
        my_obj = images.TNImage(filename=test_img)
        data = load(test_img).get_data()
        flat = data[:,:,0,0].ravel()
        assert np.mean(flat) == my_obj.z_mean()[0][0]
