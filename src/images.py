import nibabel
import numpy as np

class TNImages(object):
    """Container class for images"""

    def __init__(self):
        self._images = list()

    def get_images(self):
        return self._images
        
    def load_from_dir(self,
                      theDirectory,
                      subject_dirs = None,
                      image_filename = None,
                      bvals_filename = None,
                      bvecs_filename = None
    ):
        from os import listdir
        from os.path import join
        theFiles = list()
        if subject_dirs is None:
            theFiles = listdir(theDirectory)
            for theFile in theFiles:
                try:
                    myBval = None
                    if bvals_filename is not None:
                        myBval = join(theDirectory,bvals_filename)
                    myBvec = None
                    if bvecs_filename is not None:
                        myBvec = join(theDirectory,bvecs_filename)
                    self._images.append(TNImage(filename=join(theDirectory,theFile),
                                                bvals_file=myBval,
                                                bvecs_file=myBvec))
                except:
                    continue
        elif subject_dirs is not None:
            for sdir in subject_dirs:
                target = join(theDirectory,sdir)
                if image_filename is not None:
                    try:
                        myFilename = join(target,image_filename)
                        myBval = None
                        if bvals_filename is not None:
                            myBval = join(target,bvals_filename)
                        myBvec = None
                        if bvecs_filename is not None:
                            myBvec = join(target,bvecs_filename)
                        self._images.append(TNImage(filename=myFilename,
                                                    bvals_file=myBval,
                                                    bvecs_file=myBvec))
                    except:
                        pass
                else:
                    for theFile in listdir(target):
                        try:
                            myFilename = join(target,theFile)
                            myBval = None
                            if bvals_filename is not None:
                                myBval = join(target,bvals_filename)
                            myBvec = None
                            if bvecs_filename is not None:
                                myBvec = join(target,bvecs_filename)
                            self._images.append(TNImage(filename=myFilename,
                                                        bvals_file=myBval,
                                                        bvecs_file=myBvec))
                        except:
                            pass

    def mean_delta_z_dv(self):
        expected_shape = None
        deltas_by_subject = None
        images = self.get_images()
        for image_idx, image in enumerate(self.get_images()):
            # get our expected shape
            if expected_shape is None:
                expected_shape = image.shape
                assert len(expected_shape) >= 4
            # make sure shape of current image matches expected shape
            assert image.shape == expected_shape
            # size our results array if we haven't already done so
            if deltas_by_subject == None:
                # our shape should match the dimensions 3 and up from
                # our expected shape, plus the additional dimension of
                # subject which we will append as the first index
                stack_shape = (len(images),) + expected_shape[2::]
                deltas_by_subject = np.empty(stack_shape)
            # throw those data on the stack (sounds like a geeky rap
            # refrain)
            deltas_by_subject[image_idx] = image.delta_z_dv()
        # return the mean (through subject axis)
        return deltas_by_subject.mean(0)
                
    def mean_delta_z_dz(self):
        expected_shape = None
        deltas_by_subject = None
        images = self.get_images()
        for image_idx, image in enumerate(self.get_images()):
            # get our expected shape
            if expected_shape is None:
                expected_shape = image.shape
                assert len(expected_shape) >= 3
            # make sure shape of current image matches expected shape
            assert image.shape == expected_shape
            # size our results array if we haven't already done so
            if deltas_by_subject == None:
                # our shape should match the dimensions 3 and up from
                # our expected shape, plus the additional dimension of
                # subject which we will append as the first index
                stack_shape = (len(images),) + expected_shape[2::]
                deltas_by_subject = np.empty(stack_shape)
            # throw those data on the stack (sounds like a geeky rap
            # refrain)
            deltas_by_subject[image_idx] = image.delta_z_dz()
        # return the mean (through subject axis)
        return deltas_by_subject.mean(0)

class TNImage(object):
    """My custom image class"""

    def __init__(self,
                 filename = None,
                 name = None,
                 bvecs_file = None,
                 bvals_file = None    ):
        """
        Keyword arguments:

        filename: full path to image file, if provided will load image
        name: image name as string, for convenience
        bvecs: full path to bvec file, if provided will load bvecs
        bvals: full path to bval file, if provided will load bvals
        """
        if filename:
            self._base = nibabel.load(filename)
            self.shape = self._base.shape
        if name:
            self.name = name
        if bvecs_file:
            try:
                from numpy import loadtxt
                self.bvecs = loadtxt(bvecs_file)
            except:
                self.bvecs = None
        if bvals_file:
            try:
                from numpy import loadtxt
                self.bvals = loadtxt(bvals_file)
            except:
                self.bvals = None

    def delta_z_dv(self):
        """Return the difference of the means of the z planes as we
        move through volumes in the positive direction.

        This will return a 2 or more dimension numpy array where the
        values represent mean(z,v,...) - mean(z,v-1,...). For the
        special case where v=0, the array is "wrapped" and the value
        represents mean(z,v_0,...) - mean(z,v_last,...).

        Input must have 4 or more dimensions.

        """
        # make sure we are dealing with at least 3 dimensions
        assert len(self._base.shape) >= 4
        # get mean z
        v0 = self.z_mean()
        # create v1 where last volume is inserted before first volume,
        # then old last volume is cut
        from numpy import insert
        v1 = insert(v0,0,v0[:,-1],1)
        v1 = v1[:,0:-1]
        # return the difference
        return v0 - v1
        
    def delta_z_dz(self):
        """Return the difference of the means of the z planes as we
        move through slices in the positive direction in terms of
        z-slice and volume.

        This will return a 1 or more dimension numpy array where the
        values represent, mean(z,(v),...) - mean(z-1,(v),...). For the
        secial case where z=0, the array is "wrapped" and the value
        represents mean(z_0,(v),...)-mean(z_(last),(v),...).

        Input must have 3 or more dimensions.

        """
        # make sure we are dealing with at least 3 dimensions
        assert len(self._base.shape) >= 3
        # get mean value for each (z,v) pair
        z0 = self.z_mean()
        # create z1 where last row is inserted before first row, then
        # old last row is cut
        from numpy import insert
        z1 = insert(z0,0,z0[-1],0)
        z1 = z1[0:-1]
        # return the difference
        return z0 - z1

    def number_of_b_zeros(self):
        """Return the number of b zeros if applicable, else return
        None

        """
        if hasattr(self, 'bvals') and self.bvals is not None:
            return len(self.bvals) - np.count_nonzero(self.bvals)
        else:
            return None
            
    def v_mean(self):
        """Return the mean across the x,y,z volume for each
        permutation of remaining dimensions.

        Input must be at least a 4-dimensional image

        Output will be a numpy array with dimensionality reduced by 3,
        so a 4d input yields a 1d output, a 5d input yields a 2d
        output, etc. Remaining shape will match preserved
        dimensions.

        """
        assert len(self._base.shape) >= 4
        return self._axes_mean((0,1,2))

    def z_mean(self):
        """Return the mean across the x,y plane for each permutation
        of remaining dimensions.

        Input must be at least a 3-dimensional image

        Output will be a numpy array with dimensionaliy reduced by
        two, so a 3d input yields a 1d output, a 4d input yields a 2d
        output, etc. Remaining shape will match preserved
        dimensions.

        """
        assert len(self._base.shape) >= 3
        return self._axes_mean((0,1))

    def _axes_mean(self, axes):
        """Function adapted from:

        http://stackoverflow.com/questions/8731517/
        numpy-computing-mean-and-std-over-multiple-
        non-consecutive-axes-2nd-attempt

        Returns mean over arbitrary, possibly non-consecutive axes.

        """
        inshape = self._base.shape
        nonaxes = tuple([i for i in range(len(inshape)) if i not in set(axes)])
        permuted = self._base.get_data().transpose(nonaxes + axes)
        newshape = tuple(inshape[i] for i in nonaxes) + (-1,)
        reshaped = permuted.reshape(newshape)
        outshape = newshape[:-1]
        return reshaped.mean(axis=-1).reshape(outshape)

# module level functions
def get_wm_neighbors(ROI_img, WM_img):
    ROI_data = ROI_img._base.get_data()
    WM_data = WM_img._base.get_data()
    result = dict()
    for idx, value in np.ndenumerate(ROI_data):
        if value == 0:
            continue
        elif voxel_has_neighbor_with_value(WM_data, np.asarray(idx), 1):
            if value not in result:
                result[value] = 1
            else:
                result[value] += 1
    return result
    
def get_wm_neighbors_for_value(ROI_img, WM_img, value):
    """Return the indices of all voxels in the ROI with the given
    value, which are adjacent to a white matter voxel (a voxel with a
    value of 1 in the WM image)

    ROI image and WM image must be of equal shape
    """
    assert ROI_img._base.shape == WM_img._base.shape
    ROI_data = ROI_img._base.get_data()
    WM_data = WM_img._base.get_data()
    results = list()
    test_idxs = np.argwhere(ROI_data==value)
    for idx in test_idxs:
        if voxel_has_neighbor_with_value(WM_data, idx, 1):
            results.append(idx)
    return results

def voxel_has_neighbor_with_value(test_img_data, voxel_idx, target_value):
    search_space = [-1,0,1]
    me = tuple(voxel_idx)
    x0 = voxel_idx[0]
    y0 = voxel_idx[1]
    z0 = voxel_idx[2]
    for dx in search_space:
        for dy in search_space:
            for dz in search_space:
                test_idx = tuple([x0+dx, y0+dy, z0+dz])
                # if test_idx is equal to original, let's move on to
                # the next
                if test_idx == me:
                    continue
                # if any of our xyz coordinates (for test idx) are
                # negative, let's move on to the next (this would mean
                # we are wrapping around a dimension, like pac-man
                # when he goes through one side of the screen and ends
                # up on the opposite)
                if -1 in test_idx:
                    continue
                # if we've made it here we are ready to test for the
                # target value for this particular voxel
                try: # use try to handle indexing errors (our algo
                     # will search beyond the max index for edge
                     # cases)
                    if test_img_data[test_idx] == target_value:
                        return True # we return true and exit as soon as
                except:
                    pass # nothing to do, this is normal for edge cases
    # any match is found if we made it through the previous loop, then
    # we did not find a match
    return False
            