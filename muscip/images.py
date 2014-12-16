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

    def get_data(self):
        """Pass this call up to nibabel and return the data as
        expected from nibabel

        """
        return self._base.get_data()

    def get_header(self):
        return self._base.get_header()
        
    def get_volume(self, volume_idx):
        """Return a specified volume belonging to a 4dim image. This
        can be useful when the image is large and we want to avoid
        loading the whole image into memory.

        """
        from nibabel.volumeutils import allopen, array_from_file

        def get_fileobj(nibimage):
            fileobj = allopen(nibimage.get_filename())
            return fileobj
            
        dtype = self._base.get_data_dtype()
        shape = self._base.get_header().get_data_shape()[:-1]
        bytes_per_volume = np.prod(shape) * dtype.itemsize
        offset = self._base.get_header().get_data_offset() + (volume_idx * bytes_per_volume)
        data = array_from_file(shape, dtype, get_fileobj(self._base), offset)
        return data
        
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
def dilate_rois(roi_img, iterations=1, mask=None, output=None):
    """Dilate an roi image, by first treating the roi image as an
    signle-value image, dilate using procedure from fslmaths, then for
    every voxel which now exists in the dilated image, but not in the
    original, poll the neighborhood of new voxel in the domain of the
    original roi image, taking votes for labels that decrease in
    weight as distance from center of new voxel increases.

    Input::

      roi_img: roi image consiting of integer labeled voxels, which to
      dilate

      iterations: number of times to perform dilation; 2 iterations
      would consist of dilating once, then dilating the dilated rois
      once more, default=1

      mask: Binary mask image that must be same dimensions as
      roi_img. If provided, this will limit both the search and target
      voxels to non-zero voxels in the mask.

      output: if provided, save the result to filename provided, else
      return result as image.

    """
    affine = None
    header = None
    if isinstance(roi_img, basestring):
        roi_img = nibabel.load(roi_img)
    if isinstance(roi_img, nibabel.spatialimages.SpatialImage):
        affine = roi_img.get_affine()
        header = roi_img.get_header()
        dil_data = roi_img.get_data()
    else:
        raise Exception("Do not understand ROI image format.")
    if mask is not None:
        if isinstance(mask, basestring):
            mask = nibabel.load(mask)
        mask_data = mask.get_data()
        dil_data = (mask_data > 0) * dil_data
    roi_labels = np.unique(dil_data[dil_data > 0])
    for i in range(0,iterations):
        roi_data = dil_data.copy()
        for label in roi_labels:
            voxels_to_set = get_non_labeled_neighbors_for_value(roi_data, label)
            for voxel in voxels_to_set:
                # if another roi has already claimed this voxel, zero
                # it back out, because we don't want to bother with
                # trying to figure out which labe to assign
                if dil_data[voxel] > 0 and dil_data[voxel] != label:
                    dil_data[voxel] = 0
                # else, this is the first time we are visiting this
                # voxel let's assign the value
                else:
                    dil_data[voxel] = label
    # generate output
    out_img = nibabel.Nifti1Image(dil_data, affine, header)
    if output is None:
        return out_img
    else:
        nibabel.save(out_img, output)
        return None

def erode_rois(roi_img, iterations=1, mask=None, output=None):
    """Erode an image by zero'ing every voxel that is adjacent to a
    voxel having a value not equal to that of itself.

    Input::

      roi_img: roi image consiting of integer labeled voxels, which to
      erode

      iterations: number of times to perform dilation; 2 iterations
      would consist of eroding once, then eroding the eroded rois
      once more, default=1

      mask: Binary mask image that must be same dimensions as
      roi_img. If provided, this will limit the voxels considerd for
      eroision to the non-zero voxels in the mask

      output: if provided, save the result to filename provided, else
      return result as image.

    """
    affine = None
    header = None
    if isinstance(roi_img, basestring):
        roi_img = nibabel.load(roi_img)
    if isinstance(roi_img, nibabel.spatialimages.SpatialImage):
        affine = roi_img.get_affine()
        header = roi_img.get_header()
        roi_data = roi_img.get_data()
    else:
        raise Exception("Do not understand ROI image format.")
    if mask is not None:
        if isinstance(mask, basestring):
            mask = nibabel.load(mask)
        mask_data = mask.get_data()
    else:
        mask_data = np.ones(roi_data.shape, dtype=np.int8)
    for i in range(0,iterations):
        # for each iteration, bookmark a start point using our roi data
        from copy import copy
        start_data = copy(roi_data)
        for x in range(0, roi_data.shape[0]):
            for y in range(0, roi_data.shape[1]):
                for z in range(0, roi_data.shape[2]):
                    # if the voxel is non zero and in the inclusion mask
                    voxel_value = roi_data[x,y,z]
                    if voxel_value != 0 and mask_data[x,y,z] > 0:
                        # if the voxel has a neighbor not equal to itself,
                        # then it must be zero'd... off with its bits!
                        # ...first get valid parameters for our local
                        # neighborhood
                        low = np.array([x,y,z]) - 1
                        high = np.array([x,y,z]) + 2 # plus two, because
                        # top of range is
                        # non-inclusive
                        # handle any cases where the neighborhood range
                        # extends beyond the bounds of original image
                        # space
                        for i in xrange(3):
                            if low[i] < 0:
                                low[i] = 0
                            if high[i] > roi_data.shape[i]:
                                high[i] = roi_data.shape[i]
                        local_neighborhood = start_data[low[0]:high[0],low[1]:high[1],low[2]:high[2]]
                        # check the neighborhood for outcasts
                        if (local_neighborhood != voxel_value).any():
                            roi_data[x,y,z] = 0 # off with its bits!
    # generate output
    out_img = nibabel.Nifti1Image(roi_data, affine, header)
    if output is None:
        return out_img
    else:
        nibabel.save(out_img, output)
        return None

def generate_masks_from_roi_list(roi_list=None, fa_map=None, threshold=0.2, flip_x=False, flip_y=False, units='VOX', outdir=None):
    """For each ROI in the list, create a mask where every voxel
    inside the spherical radius of the ROI is included in the mask if
    the value for FA in the subject's FA map is greater than the given
    threshold. In this way we can provide large spherical ROIs in
    order to deal with less than perfect registration, yet still
    restrict the ROIs to meaningful anatomical areas by thresholding
    FA.

    Input::
    
      roi_list - list of spherical ROIs following the form:
      name,x,y,z,r where name represents name of the mask, xyz are
      voxel coordinates and r is the radius of the sphere in
      millimeters
      
      fa_map - subject specific FA map in the same relative space as
      ROIs are defined

      threshold - required FA value for voxel to be included in mask
      (default = 0.2)

      pixdims - pixel dimensions of the fa_map; this enables us to
      provide the radius in millimeters in the case that the 

      flip_x - if needed flip the x axis, I've encountered this before
      when defining in mricron that the x-axis needs to be flipped
      (default = False)

      flip_y - trackvis flips Anterior/Posterior, so if coordinates
      were defined using Trackvis, we need to flip our y coordinate
      (default = False)

      units - units in which coordinates and radii are defined, must
      be one of { 'VOX' | 'MM' } (default = 'VOX')

      outdir - directory in which to write masks

    """
    from math import ceil, floor, sqrt

    roi_file = open(roi_list, 'rb')
    fa_img = nibabel.load(fa_map)
    fa_data = fa_img.get_data()
    voxdims = np.asarray(fa_img.get_header().get_zooms())
    min_voxdim = voxdims.min()
    if units == 'MM':
        scale_factor = ceil(1.0 / min_voxdim)
    else:
        scale_factor = 1.0

    def _euclidian_dist(index1, index2):
        x1, y1, z1 = index1
        x2, y2, z2 = index2
        return sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )

    for entry in roi_file.readlines():
        name, x, y, z, r = entry.split(',')
        x = int( floor( float(x) * scale_factor ) )
        if flip_x:
            # flip x if required, I've needed to do this in the past
            # after defining coordinates in mricron
            x = int( voxdims[0] * scale_factor - x )
        y = int( floor( float(y) * scale_factor ) )
        if flip_y:
            # trackvis flips Anterior/Posterior, so if coordinates
            # were defined using Trackvis, we need to flip our y
            # coordinate
            y = int( voxdims[1] * scale_factor - y )
        z = int( floor( float(z) * scale_factor ) )
        r =  float(r) * scale_factor
        offset = int( ceil( r ) )
        mask_data = np.zeros(fa_img.shape)
        print "Processing ROI: %s" % name

        for x_hat in range(x-offset, x+offset):
            for y_hat in range(y-offset, y+offset):
                for z_hat in range(z-offset, z+offset):
                    if _euclidian_dist( (x,y,z), (x_hat,y_hat,z_hat) ) <= r:
                        if fa_data[x_hat,y_hat,z_hat] >= threshold:
                            mask_data[x_hat,y_hat,z_hat] = 1

        mask_img = nibabel.Nifti1Image(mask_data, fa_img.get_affine())
        nibabel.save(mask_img, '%s/%s.nii.gz' % (outdir, name))
    roi_file.close()

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

def get_non_labeled_neighbors_for_value(roi_img, value):
    """Return indices of all non-labeled voxels neighboring an roi
    with a given value.

    """
    if isinstance(roi_img, basestring):
        data = nibabel.load(roi_img).get_data()
    if isinstance(roi_img, TNImage):
        data = roi_img._base.get_data()
    if isinstance(roi_img, np.ndarray):
        data = roi_img
    results = list()
    test_idxs = np.argwhere(data == value)
    for test_idx in test_idxs:
        for found_idx in neighbors_with_value(data, test_idx, 0):
            results.append(found_idx)
    return results

def neighbors_with_value(img, voxel, value):
    """Return the indices of voxels neighboring a given voxel, having
    a given value.

    """
    if isinstance(img, basestring):
        data = nibabel.load(img).get_data()
    if isinstance(img, TNImage):
        data = img._base.get_data()
    if isinstance(img, np.ndarray):
        data = img
    results = list()
    me = tuple(voxel)
    x0 = voxel[0]
    y0 = voxel[1]
    z0 = voxel[2]
    
    def test_voxel(voxel):
        test_idx = voxel
        # if test_idx is equal to original, let’s move on to
        # the next
        if test_idx == me:
            return
        # if any of our xyz coordinates (for test idx) are
        # negative, let’s move on to the next (this would mean
        # we are wrapping around a dimension, like pac-man
        # when he goes through one side of the screen and ends
        # up on the opposite)
        if -1 in test_idx:
            return
        # if we’ve made it here we are ready to test for the
        # target value for this particular voxel
        try: # use try to handle indexing errors (our algo
            # will search beyond the max index for edge
            # cases)
            if data[test_idx] == value:
                results.append(test_idx)
        except:
            pass # nothing to do, this is normal for edge cases

    test_voxel(tuple([x0-1, y0, z0]))
    test_voxel(tuple([x0+1, y0, z0]))
    test_voxel(tuple([x0, y0-1, z0]))
    test_voxel(tuple([x0, y0+1, z0]))
    test_voxel(tuple([x0, y0, z0-1]))
    test_voxel(tuple([x0, y0, z0+1]))

    return results

def parcellate_mask(mask_img, number_of_parcels, basename=None, outdir=None):
    """Take a given mask and parcel it into a given number of parcels,
    so that each parcel contains an (almost) equal number of non-zero
    voxels. This is done to facilitate parallel processing of the type
    where each voxel is processed independently and an image can be
    divided into multiple images, and the output of each can be summed
    to obtain the equivalent output of the whole.

    Inputs::
    
      mask_img - an image containing zero and non-zero voxels to be
      dispersed amongst n-images

      number_of_parcels - number of images which to parcel out
      non-zero voxels

      basename - basename of the output images, if not provided,
      basename of original image will be used; incremented integers
      will be appended to each image

      outdir - directory in which to place images, if not provided
      will use directory at point

    """
    # load mask img
    if isinstance(mask_img, basestring):
        orig_img = nibabel.load(mask_img)
        orig_data = orig_img.get_data()
        orig_filename = mask_img
        aff = orig_img.get_affine()
    if isinstance(mask_img, TNImage):
        orig_data = mask_img._base.get_data()
        orig_filename = mask_img._base.get_filename()
        aff = orig_img._base.get_affine()
    if isinstance(mask_img, np.ndarray):
        orig_data = img
        aff = np.ones((4,4), dtype=np.int8)
    # we assume the mask to be a three dimensional image... if this is
    # not true, lets throw an exception
    if len( orig_data.shape ) != 3:
        raise Exception("For the time being, this function only operates on 3dim images.")
    # if no basename has been provided, get the base name of the mask
    # image
    if basename is None:
        import os.path as op
        basename = op.basename(mask_filename).split('.', 1)[0]
    # if no outdir has been provided use the directory at point
    if outdir is None:
        outdir = op.abspath('.')
    # round up all the non-zero voxels, and while we're at it, wez
    # better cout'em
    nonzero_voxels = np.where(orig_data != 0)
    nonzero_voxel_count = len(nonzero_voxels[0])
    # create our shiny new data (empty for now)
    new_data = np.zeros(orig_data.shape + tuple([number_of_parcels]))
    # for each index in our nonzero indices
    idx = 0
    while idx < nonzero_voxel_count:
        parcel = 0
        while parcel < number_of_parcels and idx < nonzero_voxel_count:
            x,y,z = nonzero_voxels[0][idx], nonzero_voxels[1][idx], nonzero_voxels[2][idx]
            new_data[x,y,z,parcel] = orig_data[x,y,z]
            parcel += 1
            idx += 1
    # for each of our parcels, output an image
    num_digits = len( str( number_of_parcels + 1 ) ) # determine the
                                                     # number of
                                                     # digits for our
                                                     # incrementer so
                                                     # that we can
                                                     # zero pad nicely
    for i in range( 1, number_of_parcels + 1 ):
        suffix = str(i).zfill(num_digits)
        import os.path as op
        outfile = '%s_%s.nii.gz' % ( op.join( outdir, basename ),suffix )
        outdata = new_data[...,i-1]
        new_img = nibabel.Nifti1Image(outdata, aff)
        nibabel.save(new_img, outfile)
    
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

def surface_area_for_rois(roiImage, maskImage):
    """Return an array where arr[idx] contains the count of voxels
    belonging to ROI.label == idx, that have a white matter neighbor
    as per the mask image.

    roiImage - nibabel image of ROI where voxel values represent label
    of ROI

    maskImage - nibabel image of white matter where white matter
    voxels have values of 1 and non white matter voxels have value of
    zero

    """
    ROI_data = roiImage.get_data()
    WM_data = maskImage.get_data()
    # any voxel that is classified both as ROI and WM, for this case
    # should be considered WM only, which means that if ROIs are
    # dilated, we will respect the original border between grey and
    # white matter rather than the new border resulting from dilation
    ROI_voxels = np.argwhere(WM_data != 0)
    for idx in ROI_voxels:
        if ROI_data[tuple(idx)] != 0:
            ROI_data[tuple(idx)] = 0
    # now, go ahead and get the surface area for each ROI
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

def voxel_count_by_region(roi_image, low=1, high=None):
    """Given an roi image, return a dictionary that has region labels
    for keys and voxel counts for values."""
    if high is None:
        high = int(roi_image.max())
    result = dict()
    for i in range(low, high+1):
        result[i] = len(np.where(roi_image==i)[0])
    return result
