import numpy, os, shutil
from ..connectome import TNConnectome
from ...images import TNImage
from ...fibers import fiber_length, read_trackvis, transform_fiber_by_aff, write_trackvis

__DEFAULT_MAX_FIBER_LENGTH__ = 300.0
__DEFAULT_MIN_FIBER_LENGTH__ = 20.0
__FIBERS_TO_ROI_AFFINE_FILENAME__ = 'fibers-to-roi-affine.mat'
__FIBERS_FILENAME__ = 'fibers.trk'
__TYPE__ = 'DTK'
__WM_IMAGE_FILENAME__ = 'wm.nii.gz'

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
        # load inmage
        import nibabel
        scalar_data = nibabel.load(img).get_data()
        # if affine is not provided, use inverse of fibers-to-roi-affine
        if aff is None:
            aff = numpy.linalg.inv(self.fibers_to_roi_affine)
        # create a function to map fiber indices to voxel value in
        # corresponding image
        def value_at_point(img_data, pt):
            return img_data[tuple(pt)]
        # get roi data
        roi_data = self.roi_image.get_data()
        # for every edge calculate the average scalar value for all
        # fibers
        for fiber in self.filtered_fibers:
            fiber = transform_fiber_by_aff(fiber, aff)
            u = value_at_point(roi_data, fiber[0])
            v = value_at_point(roi_data, fiber[-1])
            inside_accum = 0.0
            for pt in fiber:
                inside_accum += value_at_point(scalar_data, pt)
            avg_value_for_fiber = inside_accum / float(len(fiber))
            try:
                self.network[u][v][name] += avg_value_for_fiber
            except KeyError:
                try:
                    self.network[u][v][name] = avg_value_for_fiber
                except:
                    print "Could not find edge %s,%s - WHAT IS GOING ON?" % (u,v)
                    raise
        # finally, for each edge, we need to divide accumulations by
        # fiber count to obtain average
        for u,v,data in self.network.edges_iter(data=True):
            self.network[u][v][name] = self.network[u][v][name] / \
                                       float(self.network[u][v]['fiber_count'])
        
    @property
    def fibers(self):
        try:
            if self._fibers is not None:
                return self._fibers
        except AttributeError:
            pass
        if self.fibers_path is not None:
            return read_trackvis(self.fibers_path)
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

    @property
    def filtered_fibers(self):
        if self.fibers.fibers is not None:
            fiber_ctx = 0
            filtered_fiber_ctx = 0
            filtered_fibers_indices = self.network.graph['filtered_fibers_indices']
            filtered_fiber_count = len(filtered_fibers_indices)
            for fiber in self.fibers.fibers:
                # TODO: fix this wierd offset after I fix in generate
                # network method... indexing should start at zero, but
                # for some odd reason I started at 1 (remember to
                # confront past-self about this!)
                if fiber_ctx + 1 == filtered_fibers_indices[filtered_fiber_ctx]:
                    if filtered_fiber_ctx < filtered_fiber_count - 1 - 1:
                        filtered_fiber_ctx += 1
                        yield fiber
                    else:
                        break
                fiber_ctx += 1
            
    def generate_network(self, overwrite=False, min_length=20, max_length=300 ):
        """Generate network. Will not overwrite existing network
        unless overwrite is set to True.
        """
        if not overwrite and self.network.size() > 0:
            raise Exception("If you wish to overwrite existing \
                            network, you must explicity set overwrite=False")
        # else...
        print("...loading ROI and Fibers")
        roi_data = self.roi_image.get_data()
        fiber_counter = 0
        total_fibers = self.fibers.number_of_fibers
        self.network.graph['filtered_fibers_indices'] = []
        for fiber in self.fibers.fibers:
            fiber_counter += 1
            if fiber_counter % 10000 == 0:
                if total_fibers is not None:
                    print "...%.2f%% of fibers processed" % (100 * (fiber_counter / float(total_fibers)))
                else:
                    print "...%s fibers have been processed" % fiber_counter
            # fiber = transform_fiber_by_aff(fiber, self.fibers_to_roi_affine)
            len_fiber = fiber_length(fiber)
            if len_fiber >= min_length and len_fiber <= max_length:
                i_value = int(roi_data[tuple(fiber[0])])
                j_value = int(roi_data[tuple(fiber[-1])])
                if i_value != 0 and j_value !=0 and i_value != j_value:
                    try:
                        self.network[i_value][j_value]['fiber_count'] += 1
                        self.network[i_value][j_value]['fiber_lengths'].append(len_fiber)
                        self.network.graph['filtered_fibers_indices'].append(fiber_counter)
                    except KeyError:
                        self.network.add_edge(i_value, j_value)
                        self.network[i_value][j_value]['fiber_count'] = 1
                        self.network[i_value][j_value]['fiber_lengths'] = []
                        self.network[i_value][j_value]['fiber_lengths'].append(len_fiber)
                        self.network.graph['filtered_fibers_indices'].append(fiber_counter)
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

    def populate_hagmann_density(self):
        """Add an entry for hagman density for each entry. Need to add
        complete description and reference. This method requires the
        connectome to have both roi and wm image defined.

        """
        # check that we have required data to proceed
        if self.wm_image is None or self.roi_image is None:
            print "Need to have defined ROI and WM before density can be extracted."
            return
        # define method for summing the inverse lengths of all
        # streamlines for a given edge
        def inverse_sum(streamlines):
            inverse_streamlines = []
            for streamline in streamlines:
                inverse_streamlines.append( 1.0 / streamline )
            from math import fsum
            return fsum(inverse_streamlines)
        # get surface areas for ROIs
        from ...images import surface_area_for_rois
        surface_area = surface_area_for_rois(self.roi_image, self.wm_image)
        # make sure surface areas are well formed (no zero values) by
        # adding a small epsilon value to any zeros
        epsilon = 0.000000001
        for n in self.network.nodes():
            try:
                surface_area[n]
            except KeyError:
                surface_area[n] = epsilon
        # for every edge in network...
        for i,j in self.network.edges_iter():
            # calculate hagmann density and add to data structure
            hd = ((2.0 / (surface_area[i] + surface_area[j])) * \
                  inverse_sum(self.network[i][j]['fiber_lengths']))
            self.network[i][j]['hagmann_density'] = hd
                
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
            if self.wm_image_path is not None:
                return TNImage(filename=self.wm_image_path)
            else:
                return None
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
                    if _file == __WM_IMAGE_FILENAME__:
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
        # if filename represents a directory, create if it does not exist
        if not os.path.exists(filename):
            os.mkdir(filename)
        # SUBCLASS SPECIFIC ITEMS
        ######################################################################
        # 1) manifest
        self.manifest['MAX_FIBER_LENGTH'] = self.max_fiber_length
        self.manifest['MIN_FIBER_LENGTH'] = self.min_fiber_length
        # 2) fibers
        if self.fibers is not None:
            fiber_dest_path = os.path.join(filename, __FIBERS_FILENAME__)
            if self.fibers_path != fiber_dest_path and os.path.exists(self.fibers_path):
                shutil.copyfile(self.fibers_path, fiber_dest_path)
        # 3) fibers to roi affine
        affine_dest_path = os.path.join(filename, __FIBERS_TO_ROI_AFFINE_FILENAME__)
        numpy.savetxt(affine_dest_path, self.fibers_to_roi_affine, fmt='%f')
        # 4) wm
        if self.wm_image is not None:
            # get the file extension of the fdt file, we will use this
            # later
            # generate standard path using extension
            wm_dest_path = os.path.join(filename, __WM_IMAGE_FILENAME__)
            # if wm image path is not equal to our standard path,
            # then we need to copy the wm image to our standard path
            if self.wm_image_path != wm_dest_path and os.path.exists(self.wm_image_path):
                shutil.copyfile(self.wm_image_path, wm_dest_path)
        # SUPERCLASS WRITE
        #####################################################################
        TNConnectome.write(self, filename=filename)

    def write_filtered_fibers(self, filename):
        """Write a trackvis file populated by fibers from the
        connectome.

        """
        hdr = self.fibers.hdr
        points_space = self.fibers.spacing
        streamlines = self.filtered_fibers
        write_trackvis(filename, hdr_mapping=hdr,
                       points_space=points_space, streamlines=streamlines)

