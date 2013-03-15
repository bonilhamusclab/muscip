import networkx, numpy, os, shutil
from ..connectome import TNConnectome
from ...images import TNImage

__FDT_IMAGE_PREFIX__ = 'fdt'
__WM_IMAGE_FILENAME__ = 'wm.nii.gz'
__TYPE__ = 'PROBTRACKX'

class TNProbtrackxConnectome(TNConnectome):

    def __init__(self,
                 filename=None,
                 fdt_image=None,
                 roi_image=None,
                 wm_image=None
                 ):
        TNConnectome.__init__(self, filename=filename, roi_image=roi_image)
        self.type = __TYPE__
        if fdt_image is not None:
            self.fdt_image = fdt_image
        if wm_image is not None:
            self.wm_image = wm_image

    @property
    def fdt_image(self):
        """Return 4-dimensional fiber distribution image. This will
        most-likely have been merged from the output of probtrackx.

        """
        if self.fdt_image_path is not None:
            try:
                return TNImage(filename=self.fdt_image_path)
            except Exception, e:
                raise e
        else:
            print "Connectome has no registered fdt image."
            return None

    @fdt_image.setter
    def fdt_image(self, filename):
        """Take path to fdt image and set as our new fdt image. This
        will overwrite previous fdt_image if one exists.

        """
        self._fdt_image_path = filename

    @property
    def fdt_image_path(self):
        """Return the path to our fdt image if we have one, if not,
        return None.

        """
        try:
            return self._fdt_image_path
        except:
            try:
                found_files = os.listdir(self._filename)
                for _file in found_files:
                    if _file.split('.')[0] == __FDT_IMAGE_PREFIX__:
                        return os.path.join(self._filename,_file)
                # didn't find in connectome dir
                return None
            except:
                return None

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

    def generate_network(self, overwrite=False, matrix=None):
        """Generate network. Will not overwrite existing network
        unless overwrite is set to True.

        Input::

          overwrite: overwrite existing network
        
          matrix: If an adjacency matrix is provided, simply transfer
          the contents of the matrix to the network, rather than
          exctract all info from the fdt image. This can be done when
          matrix is pre-computed, possibly on a cluster.

        """
        if not overwrite and self.network.size() > 0:
            raise Exception("If you wish to overwrite existing \
                            network, you must explicity set overwrite=False")
        else:
            # zero-out network so that we start with clean slate
            self.network = networkx.Graph()
        if matrix is not None:
            self._generate_network_from_matrix(matrix)
            return
        # else, we are generating network from fdt and roi images...
        print("...loading ROI and Fiber Distribution images")
        roi_data = self.roi_image.get_data()
        rois = numpy.unique(roi_data[numpy.where(roi_data != 0)])
        num_rois = rois.size
        print("...%s number of regions found in atlas" % num_rois)
        seeds = numpy.arange(0,self.fdt_image.shape[-1])
        num_seeds = seeds.size
        print("...%s seed regions found in fdt image" % num_seeds)
        for seed in seeds:
            fdt_volume = self.fdt_image.get_volume(seed)
            seed = seed + 1 # seeds count from 0, roi labels from 1
            print("...counting hits for seed region: %s" % seed)            
            for target in rois:
                seed = int(seed)
                target = int(target)
                if seed == target:
                    continue # we will not consider self-loops
                # else
                # get fiber count
                fiber_count = fdt_volume[numpy.where(roi_data == target)].sum()
                # and add fiber count to edge
                self._add_edge_value_for_key(seed,target,'fiber_count',fiber_count)
            
    def populate_density(self):
        """Add an entry for density for each edge. TODO: add complete
        description and reference. This method requires the connectome
        to have both roi and wm image defined."""
        # check that both wm and roi are defined
        if self.wm_image is None or self.roi_image is None:
            print "Need to have defined ROI and WM before density can be extracted."
            return
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
            density = self.network[i][j]['fiber_count'] / \
                float( surface_area[i] + surface_area[j] )
            self.network[i][j]['density'] = density

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
        # first let the parent object write
        TNConnectome.write(self, filename=filename)
        # then write our additions
        # FDT IMAGE (if needed) #############################################
        if self.fdt_image_path is not None:
            # get the file extension of the fdt file, we will use this
            # later
            fdt_file_ext = self.fdt_image_path.split('.',1)[1]
            # generate standard path using extension
            fdt_dest_path = os.path.join(filename, "%s.%s" %
                                         (__FDT_IMAGE_PREFIX__, fdt_file_ext))
            # if fdt image path is not equal to our standard path,
            # then we need to copy the fdt image to our standard path
            if self.fdt_image_path != fdt_dest_path:
                shutil.copyfile(self.fdt_image_path, fdt_dest_path)
        # WM IMAGE (if needed) 
        if self.wm_image is not None:
            # get the file extension of the fdt file, we will use this
            # later
            # generate standard path using extension
            wm_dest_path = os.path.join(filename, __WM_IMAGE_FILENAME__)
            # if wm image path is not equal to our standard path,
            # then we need to copy the wm image to our standard path
            if self.wm_image_path != wm_dest_path and os.path.exists(self.wm_image_path):
                shutil.copyfile(self.wm_image_path, wm_dest_path)


    def _add_edge_value_for_key(self,m,n,key,value):
        # 1) Edge exists
        try:
            # This is the case where we have already recorded i -> j
            # and now must figure the average of i -> j and j ->
            # i. Because this graph type is non-directed, setting ij
            # also effectively sets ji
            self.network[m][n][key] += value
            self.network[m][n][key] /= 2.0
        # 2) Edge does not exist
        except KeyError:
            # In this case we can simply add the edge along
            # with the given fiber count
            self.network.add_edge(m, n, {key: value})

    def _generate_network_from_matrix(self, matrix):
        if len(matrix.shape) != 2 or matrix.shape[0] != matrix.shape[1]:
            raise Exception("Matrix must be an NxN adjancency matrix")
        elements = matrix.shape[0]
        for m in range(0,elements):
            for n in range(0,elements):
                if m == n:
                    continue
                if matrix[m][n] > 0:
                    # add (to) edge
                    self._add_edge_value_for_key(m+1,n+1,'fiber_count',matrix[m][n])
