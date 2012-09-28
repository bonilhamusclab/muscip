import numpy, os, shutil
from ..connectome import TNConnectome
from ...images import TNImage
from ...fibers import length_of_streamline, read, transform_streamline_by_aff

__DEFAULT_MAX_FIBER_LENGTH__ = 200.0
__DEFAULT_MIN_FIBER_LENGTH__ = 20.0
__FIBERS_FILENAME__ = 'fibers.trk'
__FIBERS_TO_ROI_AFFINE_FILENAME__ = 'fibers-to-roi-affine.mat'
__TYPE__ = 'DTK'
__WM_IMAGE_PREFIX__ = 'wm'

class TNDtkConnectome(TNConnectome):
    
    def __init__(self,
                 fibers=None,
                 fibers_to_roi_affine=None,                 
                 filename=None,
                 max_fiber_length=None,
                 min_fiber_length=None,
                 roi_image=None,
                 wm_image=None
                 ):
        TNConnectome.__init__(self, filename=filename, roi_image=roi_image)
        self.type = __TYPE__
        if fibers is not None:
            self.fibers = fibers
        if fibers_to_roi_affine is not None:
            self.fibers_to_roi_affine = fibers_to_roi_affine
        if wm_image is not None:
            self.wm_image = wm_image
        if max_fiber_length is not None:
            self.max_fiber_length = max_fiber_length
        if min_fiber_length is not None:
            self.min_fiber_length = min_fiber_length

    def add_scalar(self, img, name, aff=None):
        # load image
        import nibabel
        scalar_data = nibabel.load(img).get_data()
        # if affine is not provided, use inverse of fibers-to-roi-affine
        if aff is None:
            aff = numpy.linalg.inv(self.fibers_to_roi_affine)
        # for every edge calculate the average scalar value for all
        # streamlines
        for i,j in self.network.edges_iter():
            accum_outside = 0.0
            for streamline in self.network[i][j]['streamlines']:
                streamline = transform_streamline_by_aff(streamline, aff, self.vox_dims)
                accum_inside = 0.0
                for idx in streamline:
                    try:
                        accum_inside += scalar_data[tuple(idx)]
                    except:
                        print "No data found at %s" % idx
                        continue
                accum_outside += accum_inside
            self.network[i][j][name] = accum_outside
        
    @property
    def fibers(self):
        try:
            if self._fibers is not None:
                return self._fibers
        except AttributeError:
            pass
        if self.fibers_path is not None:
            return read(self.fibers_path, format='trackvis')
        else:
            return None
        
    @fibers.setter
    def fibers(self, filename):
        """Take path to fibers file and set as our new fibers. This will
        overwrite previous fibers if one exists.
        
        """
        self._fibers = None
        self._fibers_file_path = filename

    @property
    def fibers_path(self):
        try:
            return self._fibers_file_path
        except:
            try:
                found_files = os.listdir(self._filename)
                for _file in found_files:
                    if _file == __FIBERS_FILENAME__:
                        return os.path.join(self._filename, _file)
                # else, we did not find fibers file in connectome dir
                return None
            except:
                return None

    @property
    def fibers_to_roi_affine(self):
        try:
            return self._fibers_to_roi_affine
        except AttributeError:
            try:
                return numpy.loadtxt(self.fibers_to_roi_affine_path)
            except:
                return numpy.eye(4) # return identity matrix if none exists

    @fibers_to_roi_affine.setter
    def fibers_to_roi_affine(self, new_affine):
        try:
            new_affine = numpy.loadtxt(new_affine)
            if not new_affine.shape == (4,4):
                raise Exception("Affine should be a 4x4 matrix.")
            self._fibers_to_roi_affine = new_affine
        except:
            raise Exception("Could not parse affine.")

    @property
    def fibers_to_roi_affine_path(self):
        try:
            return os.path.join(self._filename, __FIBERS_TO_ROI_AFFINE_FILENAME__)
        except:
            return None

    def generate_network(self, overwrite=False, min_length=20, max_length=240 ):
        """Generate network. Will not overwrite existing network
        unless overwrite is set to True.
        """
        if not overwrite and self.network.size() > 0:
            raise Exception("If you wish to overwrite existing \
                            network, you must explicity set overwrite=False")
        # else...
        print("...loading ROI and Fibers")
        roi_data = self.roi_image.get_data()
        streamline_counter = 0
        total_streamlines = self.fibers.number_of_streamlines
        for streamline in self.fibers.streamlines:
            streamline_counter += 1
            if streamline_counter % 10000 == 0:
                if total_streamlines is not None:
                    print "...%.2f of fibers processed" % (100 * (streamline_counter / float(total_streamlines)))
                else:
                    print "...%s fibers have been processed" % streamline_counter
            streamline = transform_streamline_by_aff(streamline, self.fibers_to_roi_affine)
            len_streamline = length_of_streamline(streamline)
            if len_streamline >= min_length and len_streamline <= max_length:
                i_value = int(roi_data[tuple(streamline[0])])
                j_value = int(roi_data[tuple(streamline[-1])])
                if i_value != 0 and j_value !=0 and i_value != j_value:
                    try:
                        self.network[i_value][j_value]['fiber_count'] += 1
                        # self.network[i_value][j_value]['streamlines'].append(streamline)
                    except KeyError:
                        self.network.add_edge(i_value, j_value)
                        self.network[i_value][j_value]['fiber_count'] = 1
                        # self.network[i_value][j_value]['streamlines'] = []
                        # self.network[i_value][j_value]['streamlines'].append(streamline)

    @property
    def max_fiber_length(self):
        try:
            return self._max_fiber_length
        except AttributeError:
            return __DEFAULT_MAX_FIBER_LENGTH__

    @max_fiber_length.setter
    def max_fiber_length(self, value):
        self._max_fiber_length = value

    @property
    def min_fiber_length(self):
        try:
            return self._min_fiber_length
        except AttributeError:
            return __DEFAULT_MIN_FIBER_LENGTH__

    @min_fiber_length.setter
    def min_fiber_length(self, value):
        self._min_fiber_length = value

    @property
    def scalars(self):
        try:
            if self._scalars is not None:
                return self._scalars
        except AttributeError:
            self._scalars = dict()
            return self._scalars
    
    @property
    def wm_image(self):
        try:
            return TNImage(filename=self.wm_image_path)
        except:
            return None

    @wm_image.setter
    def wm_image(self, filename):
        self._wm_image_path = filename

    @property
    def wm_image_path(self):
        try:
            return self._wm_image_path
        except AttributeError:
            try:
                found_files = os.listdir(self._filename)
                for _file in found_files:
                    if _file.split('.')[0] == __WM_IMAGE_PREFIX__:
                        return os.path.join(self._filename, _file)
                    # else, we did not find an wm file
                    return None
            except:
                return None

    def write(self, filename=None):
        """Write connectome to the given filename, do not need to
        provide filename if we are saving to same file from which we
        read.

        """
        # get filename if we need...
        if filename is None:
            filename = self._filename
        # if filename is still None, we have a problem...
        if filename is None:
            raise Exception("No filename has been provided.")
        # SUBCLASS SPECIFIC ITEMS
        ######################################################################
        # 1) manifest
        self.manifest['MAX_FIBER_LENGTH'] = self.max_fiber_length
        self.manifest['MIN_FIBER_LENGTH'] = self.min_fiber_length
        # 2) fibers
        if self.fibers is not None:
            fiber_dest_path = os.path.join(filename, __FIBERS_FILENAME__)
            if self.fibers_path != fiber_dest_path:
                self.fibers.write(fiber_dest_path)
        # 3) fibers to roi affine
        affine_dest_path = os.path.join(filename, __FIBERS_TO_ROI_AFFINE_FILENAME__)
        numpy.savetxt(affine_dest_path, self.fibers_to_roi_affine, fmt='%f')
        # 4) wm
        if self.wm_image is not None:
            # get the file extension of the fdt file, we will use this
            # later
            wm_file_ext = self.wm_image_path.split('.',1)[1]
            # generate standard path using extension
            wm_dest_path = os.path.join(filename, "%s.%s" %
                                         (__WM_IMAGE_PREFIX__, wm_file_ext))
            # if wm image path is not equal to our standard path,
            # then we need to copy the wm image to our standard path
            if self.wm_image_path != wm_dest_path:
                shutil.copyfile(self.wm_image_path, wm_dest_path)
        # SUPERCLASS WRITE
        #####################################################################
        TNConnectome.write(self, filename=filename)
