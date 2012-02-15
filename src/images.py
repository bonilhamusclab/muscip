import nibabel

class TNImages(object):
    """Container class for images"""

    def __init__(self):
        self._images = list()

    def get_images(self):
        return self._images
        
    def load_from_dir(self, theDirectory):
        from os import listdir
        from os.path import join
        theFiles = listdir(theDirectory)
        for theFile in theFiles:
            try:
                self._images.append(TNImage(filename=join(theDirectory,theFile)))
            except:
                continue

class TNImage(object):
    """My custom image class"""

    def __init__(self,
                 filename=None):
        if filename:
            self._base = nibabel.load(filename)

    def dz(self):
        """Return the difference of the means of the z planes as we
        move through slices in the positive direction in terms of
        z-slice and volume.

        This will return a 2-dimensional numpy array where the values
        represent, mean(z,v) - mean(z-1,v). For the secial case where
        z=0, the array is "wrapped" and the value represents
        mean(z_0,v)-mean(z_(last),v).

        Input must be a 4-dimensional image.

        """
        # make sure we are dealing in 4 dimensions        
        assert len(self._base.shape) == 4
        # get mean value for each (z,v) pair
        z0 = self.z_mean()
        # create z1 where last row is inserted before first row, then
        # old last row is cut
        from numpy import insert
        z1 = insert(z0,0,z0[-1,:],0)
        z1 = z1[0:-1,:]
        # return the difference
        return z0 - z1
        
    def z_mean(self):
        """Return the mean across the x,y plane for each permutation
        of remaining dimensions.

        Input must be at least a 3-dimensional image

        Output will be a numpy array with dimensionaliy reduced by
        two, so a 3d input yields a 1d output, a 4d input yields a 2d
        output, etc. Remaining shape will match preserved
        dimensions.

        """
        inshape = self._base.shape
        assert len(inshape) >= 3
        axes = (0,1)
        nonaxes = tuple(range(2,len(inshape)))
        permuted = self._base.get_data().transpose(nonaxes + axes)
        newshape = tuple(inshape[i] for i in nonaxes) + (-1,)
        reshaped = permuted.reshape(newshape)
        outshape = newshape[:-1]
        return reshaped.mean(axis=-1).reshape(outshape)
        

        
        