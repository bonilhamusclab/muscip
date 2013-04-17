from .fibers import TNFibers

__FORMAT__ = 'dke'

class TNDKEFibers(TNFibers):

    def __init__(self, **kwargs):
        TNFibers.__init__(self, **kwargs)
        self._format = __FORMAT__
        self._fiber_key = kwargs.get('fiber_key', None)
        self._matfile = kwargs.get('matfile', None)
        self.ref_image = kwargs.get('ref_image', None)
        self._store_fibers = kwargs.get('store_fibers', True)
        # fibers must be set last, because it may depend on above
        # params
        self._fibers = kwargs.get('fibers', None)

    @property
    def fibers(self):
        # define a function that yields one fiber at a time for all
        # fibers
        def fibers_generator(all_fibers):
            for fiber in all_fibers:
                yield fiber
        # if we already have fibers loaded, then simply feed those
        # saved fibers into our fiber generator
        if self._fibers is not None:
            return fibers_generator(self._fibers)
        # if we have not returned anything yet, let's try to return
        # fibers from the matlab file
        if self._matfile is not None:
            # load the matlab file
            import scipy.io
            tmp_raw = scipy.io.loadmat(self._matfile)
            # get our fiber key
            if self._fiber_key is not None:
                my_fiber_key = self._fiber_key
            else:
                my_fiber_key = self._guess_fiber_key(tmp_raw)
            # use this information to get the fibers... in order to do
            # so, we assume that stored under the fiber key, is a list
            # of length N, where each entry is a list of length M,
            # where each entry is a 3-tuple:
            # 
            # fibers --> single fiber --> list of points as [ x y z ]
            ##########################################################
            # store fibers if we have been told to do so, discard if
            # not... in either case return fibers as a generator
            if self._store_fibers:
                self._fibers = tmp_raw[my_fiber_key][0]
                return fibers_generator(self._fibers)
            else:
                return fibers_generator(tmp_raw[my_fiber_key][0])

    @fibers.setter
    def fibers(self, new_fibers):
        """List of streamlines, where each streamline is a list of
        3-tuples representing xyz-coordinates."""
        self._fibers = new_fibers

    @property
    def fiber_key(self):
        """String representing the name the fiber data structure was
        saved under matlab."""
        return self._fiber_key

    @fiber_key.setter
    def fiber_key(self, new_fiber_key):
        """String representing the name the fiber data structure was
        saved under matlab."""
        self._fiber_key = new_fiber_key

    @property
    def matfile(self):
        """Path to matlab file containing the fiber information."""
        return self._matfile

    @matfile.setter
    def matfile(self, new_matfile):
        """Path to matlab file containing the fiber information."""
        self._matfile = new_matfile

    @property
    def ref_image(self):
        """Image which defines the space in which these fibers are
        defined. Examples may be the B0 (nodif), FA, etc. Basically
        anything in diffusion space. May be set either from a nibabel
        image object, or a filepath."""
        return self._ref_image

    @ref_image.setter
    def ref_image(self, new_ref_image):
        """Image which defines the space in which these fibers are
        defined. Examples may be the B0 (nodif), FA, etc. Basically
        anything in diffusion space. May be set either from a nibabel
        image object, or a filepath."""
        if type(new_ref_image) == str:
            import nibabel
            self._ref_image = nibabel.load(new_ref_image)
        else:
            self._ref_image

    @property
    def store_fibers(self):
        """Boolean indicating whether fibers should be kept in memory,
        or discarded after reading."""
        return self._store_fibers

    @store_fibers.setter
    def store_fibers(self, new_store_fibers):
        """Boolean indicating whether fibers should be kept in memory,
        or discarded after reading."""
        self._new_store_fibers = new_store_fibers

    def write_to_trackvis(self, filename):
        """Write fibers to given filename in trackvis format."""
        # make sure that we have a reference image, because we're
        # gonna need one
        if self.ref_image is None:
            raise Exception('Reference image is needed to write to ' + 
                            'trackvis format.')
        # now let's get to writing our trackvis file
        import nibabel.trackvis
        # create an empty trackvis header
        hdr = nibabel.trackvis.empty_header()
        # apply reference image's affine to the header
        nibabel.trackvis.aff_to_hdr(self.ref_image.get_affine(),
                                    hdr, pos_vox=True,
                                    set_order=True)
        # apply reference image's dim info to the header
        hdr['dim'] = self.ref_image.shape
        # zero the origin to prevent offset
        hdr['origin'] = [0,0,0]
        # now we can actually write the file
        from .tn_trackvis_fibers import write_trackvis
        write_trackvis(filename,
                       hdr_mapping=hdr,
                       points_space='voxel',
                       streamlines=self.fibers)

    def _guess_fiber_key(self, matlab_dict):
        all_keys = matlab_dict.keys()
        # we will make an assumption that the key does not begin with
        # an underscore
        all_keys = [key for key in all_keys if key[0] != '_']
        # if we have only one key left, we will assume that is the key
        # we are looking for
        if len(all_keys) == 1: return all_keys[0]
        # otherwise, lets further assume that the key may contain the
        # word track
        all_keys = [key for key in all_keys if 'track' in key.lower()]
        # again, if we have only one key left, we will assume we have
        # found our golden ticket
        if len(all_keys) == 1: return all_keys[0]
        # STILL NO!? Okay, let's make one more assumption that the key
        # should also contain the string 'filter'
        all_keys = [key for key in all_keys if 'filter' in key.lower()]
        # let's perform the same check, but this time, if we have not
        # found our match we must give up :(
        if len(all_keys) == 1:
            return all_keys[0]
        else:
            raise Exception('I could not guess the field name for ' +
                            'your fiber structure')
            return None

#######################################################################

def read_dke(matfile, ref_image, fiber_key=None):
    from .tn_dke_fibers import TNDKEFibers
    return TNDKEFibers(matfile=matfile,
                       ref_image=ref_image,
                       fiber_key=fiber_key)
