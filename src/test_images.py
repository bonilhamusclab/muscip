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

    def test_does_not_choke_on_non_image_file(self):
        error = False
        try:
            self._i.load_from_dir(SUPPORT_DIR)
        except:
            error = True
        assert not error

    def test_images_are_tn_images(self):
        self._i.load_from_dir(SUPPORT_DIR)
        test_img = self._i.get_images()[0]
        assert type(test_img) is images.TNImage

    def test_can_get_mean_delta_z_dz(self):
        from os.path import join
        self._i.load_from_dir(join(SUPPORT_DIR,'4dim_images'))
        assert self._i.mean_delta_z_dz() is not None

    def test_mean_delta_z_dz_is_expected_shape(self):
        from os.path import join
        self._i.load_from_dir(join(SUPPORT_DIR,'4dim_images'))
        result = self._i.mean_delta_z_dz()
        assert result.shape == (10,10)

    def test_mean_delta_z_dz_3_dim(self):
        from os.path import join
        self._i.load_from_dir(join(SUPPORT_DIR,'3dim_images'))
        result = self._i.mean_delta_z_dz()
        assert result.shape == (10,)

    def test_can_get_mean_delta_z_dv(self):
        from os.path import join
        self._i.load_from_dir(join(SUPPORT_DIR,'4dim_images'))
        assert self._i.mean_delta_z_dv() is not None

    def test_mean_delta_z_dv_is_expected_shape(self):
        from os.path import join
        self._i.load_from_dir(join(SUPPORT_DIR,'4dim_images'))
        result = self._i.mean_delta_z_dv()
        assert result.shape == (10,10)

    def test_can_load_from_directory_with_subjects_and_b_files(self):
        from os.path import join
        basepath = join(SUPPORT_DIR,'subject_dir')
        from os import listdir
        subjects = listdir(basepath)
        self._i.load_from_dir(basepath, subject_dirs = subjects,
                              image_filename = 'rdata.nii.gz',
                              bvals_filename = 'bvals',
                              bvecs_filename = 'bvecs')
        # get the images we found
        images = None
        try:
            images = self._i.get_images()
        except Exception as e:
            print e
        # make sure we got something back
        assert images is not None
        assert len(images) >= 1

    def test_can_load_from_directory_with_params_and_get_num_of_b_zeros(self):
        from os.path import join
        basepath = join(SUPPORT_DIR,'subject_dir')
        from os import listdir
        subjects = listdir(basepath)
        self._i.load_from_dir(basepath, subject_dirs = subjects,
                              image_filename = 'rdata.nii.gz',
                              bvals_filename = 'bvals',
                              bvecs_filename = 'bvecs')
        # each of these images should have 2 bvals
        for image in self._i.get_images():
            assert image.number_of_b_zeros() == 2

            
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

    def test_can_get_delta_z_dz(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert (my_image.delta_z_dz() == 0).all()

    def test_can_get_delta_z_dz_for_3dim(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'3dim_images/3dim.nii.gz'))
        assert my_image.delta_z_dz() is not None
        
    def test_delta_z_dz_is_accurate(self):
        from os.path import join
        test_img = join(SUPPORT_DIR,'noise.nii.gz')
        my_obj = images.TNImage(filename=test_img)
        mean_z = my_obj.z_mean()
        delta_z_dz = my_obj.delta_z_dz()
        assert (mean_z[1,:] - mean_z[0,:] == delta_z_dz[1,:]).all()
        assert (mean_z[0,:] - mean_z[-1,:] == delta_z_dz[0,:]).all()

    def test_can_get_z_mean(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert my_image.z_mean() is not None

    def test_can_get_z_mean_for_3dim(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'3dim_images/3dim.nii.gz'))
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

    def test_can_get_v_mean(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert my_image.v_mean() is not None

    def test_v_mean_is_accurate(self):
        from nibabel import load
        from os.path import join
        import numpy as np
        test_img = join(SUPPORT_DIR,'noise.nii.gz')
        my_obj = images.TNImage(filename=test_img)
        data = load(test_img).get_data()
        flat = data[:,:,:,0].ravel()
        assert np.mean(flat) == my_obj.v_mean()[0]
    
    def test_can_get_delta_z_dv(self):
        from os.path import join
        my_image = images.TNImage(join(SUPPORT_DIR,'one.nii.gz'))
        assert my_image.delta_z_dv() is not None

    def test_delta_z_dv_is_accurate(self):
        from os.path import join
        test_img = join(SUPPORT_DIR,'noise.nii.gz')
        my_obj = images.TNImage(filename=test_img)
        mean_z = my_obj.z_mean()
        delta_z_dv = my_obj.delta_z_dv()
        assert (mean_z[:,1] - mean_z[:,0] == delta_z_dv[:,1]).all()
        assert (mean_z[:,0] - mean_z[:,-1] == delta_z_dv[:,0]).all()
        
    def test_can_set_number_of_b_zeros(self):
        my_obj = images.TNImage()
        my_obj.number_of_b_zeros = 5
        assert my_obj.number_of_b_zeros == 5

    def test_access_number_of_b_zeros_without_set_no_error(self):
        my_obj = images.TNImage()
        error = False
        try:
            my_obj.number_of_b_zeros()
        except Exception as e:
            print e
            error = True
        assert not error

    def test_can_get_number_of_b_zeros(self):
        subject_dir = 'test_data/images/subject_dir/SUB_0001/'
        my_obj = images.TNImage(filename = subject_dir + 'rdata.nii.gz',
                                bvecs_file = subject_dir + 'bvecs',
                                bvals_file = subject_dir + 'bvals')
        # this particular subject has 2 b0's
        assert my_obj.number_of_b_zeros() == 2

    def test_can_get_wm_neighbors_for_value(self):
        # for a test we will provide a 3dim white-matter mask with one
        # zero value in the xyz-center of the cube...  this
        # configuration should yield 26 neighbors (entire cube
        # centered around center voxel, less the voxel itself).  The
        # ROI image will consist of all ones, being of the same shape
        # as the mask.
        roiImage = images.TNImage(filename = 'test_data/images/3dim_images/ones.nii.gz')
        maskImage = images.TNImage(filename = 'test_data/images/3dim_images/wm_mask_only_center.nii.gz')
        assert len(images.get_wm_neighbors_for_value(roiImage,maskImage,1)) == 26

    def test_can_get_wm_neighbors(self):
        # use the assumptions of the wm_neighbors_for_value test but
        # return the results in a plural context, even though these
        # test data will only yield one value
        roiImage = images.TNImage(filename = 'test_data/images/3dim_images/ones.nii.gz')
        maskImage = images.TNImage(filename = 'test_data/images/3dim_images/wm_mask_only_center.nii.gz')
        result = images.get_wm_neighbors(roiImage, maskImage)
        assert len(result) == 1
        assert result[1] == 26
