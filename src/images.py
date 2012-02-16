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
        
        theFiles = listdir(theDirectory)
        for theFile in theFiles:
            try:
                self._images.append(TNImage(filename=join(theDirectory,theFile)))
            except:
                continue

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

        