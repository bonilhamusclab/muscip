import networkx, numpy, os, shutil
from ..connectome import TNConnectome
from ...images import TNImage

__FDT_IMAGE_PREFIX__ = 'fdt'
__TYPE__ = 'PROBTRACKX'

class TNProbtrackxConnectome(TNConnectome):

    def __init__(self,
                 filename=None,
                 fdt_image=None,
                 roi_image=None ):
        TNConnectome.__init__(self, filename=filename, roi_image=roi_image)
        self.type = __TYPE__
        if fdt_image is not None:
            self.fdt_image = fdt_image

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

    def generate_network(self, overwrite=False):
        """Generate network. Will not overwrite existing network
        unless overwrite is set to True.

        """
        if not overwrite and self.network.size() > 0:
            raise Exception("If you wish to overwrite existing \
                            network, you must explicity set overwrite=False")
        # else...
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
                # 1) Edge exists
                try:
                    # This is the case where we have already recorded
                    # i -> j and now must figure the average of i -> j
                    # and j -> i. Because this graph type is
                    # non-directed, setting ij also effectively sets
                    # ji
                    self.network[seed][target]['fiber_count'] += fiber_count
                    self.network[seed][target]['fiber_count'] /= 2.0
                # 2) Edge does not exist
                except KeyError:
                    # In this case we can simply add the edge along
                    # with the given fiber count
                    self.network.add_edge(seed, target, {'fiber_count': fiber_count})
            
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
