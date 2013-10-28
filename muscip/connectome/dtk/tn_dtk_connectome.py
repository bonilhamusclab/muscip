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
        if float(self.version) < 3.0:
            print "Upgrading connectome to 3.0... re-generating network"
            self.generate_network(overwrite=True)
            self.version = 3.0

    def add_scalar(self, img, name):
        # load inmage
        import nibabel
        scalar_data = nibabel.load(img).get_data()
        # create a function to map fiber indices to voxel value in
        # corresponding image
        def value_at_point(img_data, pt):
            return img_data[tuple(pt)]

        roi_data = self.roi_image.get_data() # get roi data

        # for every edge calculate the average scalar value for all
        # fibers
        for fiber in self.filtered_fibers:
            u = int(value_at_point(roi_data, fiber[0]))
            v = int(value_at_point(roi_data, fiber[-1]))
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
        
    def centroids_for_region_pairs(self, min_numfib=10, avg_length_dict=None, group_connectome=None):
        """Return a dictionary where streamlines between region a and
        region b are represented as a single list of points
        representing the centroid.

        PARAMS
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        min_numfib - create a centroid for a region pair if and only
                     if the fiber count for said pair is greater than
                     or equal than this number.

        avg_lengtgh_dict - A dictionary containing the average length
                           for all possible pairs (u,v). This may be
                           useful when one wants to constrain the
                           centroids for multiple subjects to the same
                           number of points (which is based on
                           length). If both an avg_length and
                           group_connectome is provided, the
                           avg_length will take precedence.

        group_connectome - if defined, the average length will be
                           taken from this connectome rather than the
                           one we are currently operating upon. this
                           is useful in group analysis when we might
                           like to force all the centroids of
                           individuals to have the same number of
                           points, and we do so by feeding our
                           centroid generating function the same
                           length factor. If both an avg_length and
                           group_connectome is provided, the
                           avg_length will take precedence.

        """

        from ...fibers import centroid_for_tracks

        # start with an empty dictionary for our centroids
        centroids = dict()

        for u,v,data in self.network.edges_iter(data=True):
            print "Generating centroid for regions [ %s -- %s ]" % (u,v)
            if data['fiber_count'] >= min_numfib:
                # if an avg_length dictionary was provided, use avg
                # length defined for u,v
                if avg_length_dict is not None:
                    if u in avg_length_dict.keys():
                        avg_length = avg_length_dict[u][v]
                    else:
                        avg_length = avg_length_dict[v][u]
                    thisCentroid = centroid_for_tracks(self.tracks_for_regions(u,v),
                                                       length_factor=avg_length)
                # else, if a group connectome was provided, use the
                # average of fibers for u,v as the length factor for
                # the centroid generation
                elif group_connectome is not None:
                    avg_length = numpy.mean(group_connectome.network[u][v]['fiber_lengths'])
                    thisCentroid = centroid_for_tracks(self.tracks_for_regions(u,v),
                                                       length_factor=avg_length)
                # otherwise use the default, which is the average
                # fiber length of this connectome
                else:
                    thisCentroid = centroid_for_tracks(self.tracks_for_regions(u,v))

                try:
                    centroids[u][v] = thisCentroid
                except KeyError:
                    centroids[u] = dict()
                    centroids[u][v] = thisCentroid
        return centroids

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
    def filtered_fibers(self, points_space='voxel'):
        """Return fibers counted in the connectome in one of three
        spaces: (voxel, voxmm, rasmm) default: voxel"""
        for u,v,data in self.network.edges_iter(data=True):
            for fiber in data['fibers']:
                if points_space == 'voxel':
                    yield fiber
                elif points_space == 'voxmm':
                    yield fiber * self.fibers.voxel_size
                elif points_space == 'rasmm':
                    yield transform_fiber_by_aff(fiber, self.fibers.vox_to_ras)
                else:
                    raise Exception("%s is not a valid points space")

    def generate_network(self, overwrite=False):
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
                if total_fibers is not None and total_fibers != 0:
                    print "...%.2f%% of fibers processed" % (100 * (fiber_counter / float(total_fibers)))
                else:
                    print "...%s fibers have been processed" % fiber_counter
            # fiber = transform_fiber_by_aff(fiber, self.fibers_to_roi_affine)
            len_fiber = fiber_length(fiber, self.fibers.voxel_size)
            try:
                if len_fiber >= self.min_fiber_length and len_fiber <= self.max_fiber_length:
                    i = tuple(fiber[0])        # index of first point
                    j = tuple(fiber[-1])       # index of last point
                    i_value = int(roi_data[i]) # get label under i
                    j_value = int(roi_data[j]) # get label under j
                    if i_value != 0 and j_value !=0 and i_value != j_value:
                        try:
                            self.network[i_value][j_value]['fiber_count'] += 1
                            self.network[i_value][j_value]['fiber_lengths'].append(len_fiber)
                            self.network[i_value][j_value]['fibers'].append(fiber)
                            
                        except KeyError:
                            self.network.add_edge(i_value, j_value)
                            self.network[i_value][j_value]['fiber_count'] = 1
                            self.network[i_value][j_value]['fiber_lengths'] = []
                            self.network[i_value][j_value]['fiber_lengths'].append(len_fiber)
                            self.network[i_value][j_value]['fibers'] = []
                            self.network[i_value][j_value]['fibers'].append(fiber)

            except IndexError:
                print "i: %s, j: %s" % (i,j)
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
        # get surface areas for ROIs
        from ...images import surface_area_for_rois
        surface_area = surface_area_for_rois(self.roi_image, self.wm_image)
        # make sure surface areas are well formed (no zero values) by
        # adding treating any zero value as a default value (1)
        default_surface_area = 1
        for n in self.network.nodes():
            try:
                surface_area[n]
            except KeyError:
                surface_area[n] = default_surface_area
        # for every edge in network...
        for i,j in self.network.edges_iter():
            # calculate hagmann density and add to data structure
            combined_surface_area = surface_area[i] + surface_area[j]
            inverse_lengths = 1. / numpy.asarray(self.network[i][j]['fiber_lengths'])
            self.network[i][j]['hagmann_density'] = ( 2. / combined_surface_area ) * inverse_lengths.sum()
                
    @property
    def scalars(self):
        try:
            if self._scalars is not None:
                return self._scalars
        except AttributeError:
            self._scalars = dict()
            return self._scalars
    
    def tracks_for_regions(self, u, v):
        """Return a list of fibers for a given u,v pair. If no fibers
        exist for this pair, return an empty list"""
        try:
            return self.network[u][v]['fibers']
        except:
            return []

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

