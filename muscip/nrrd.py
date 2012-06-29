

class BaseNrrdData(object):

    valid_formats = ["raw", "txt", "text", "ascii", "hex", "gz",
                     "gzip"]

    def __init__(self, data, format='gzip'):
        if not format in BaseNrrdData.valid_formats:
            raise Exception("%s is not a valid data format" % format)
        self._data = data
        self._format = format

    @property
    def data(self):
        return self._data

    @property
    def format(self):
        return self._format

    def write_to_file(self, fileobj, dtype):
        raw = self.data.astype(dtype)
        if self.format in ['gzip', 'gz']:
            self._write_gzip_to_file(fileobj, raw)
        if self.format in ['raw', 'txt', 'text', 'ascii']:
            self._write_txt_to_file(fileobj, raw)
        if self.format == 'hex':
            self._write_hex_to_file(fileobj, raw)

    def _write_gzip_to_file(self, fileobj, raw):
        import gzip
        f = gzip.GzipFile(mode='wb', fileobj=fileobj)
        try:
            f.write(raw)
        except Exception as e:
            raise e

    def _write_txt_to_file(self, fileobj, raw):
        try:
            fileobj.write(raw)
        except Exception as e:
            raise e
        
    def _write_hex_to_file(self, raw):
        ##TODO: implement hex write
        pass


class BaseNrrdHeader(dict):
    """Nrrd header"""

    mandatory_keys = ['type', 'dimension', 'space', 'sizes',
                      'space directions', 'endian',
                      'encoding', 'space origin']
    valid_spaces = ["right-anterior-superior",
                    "right-posterior-superior", "left-anterior-superior",
                    "left-posterior-superior"]
    valid_types = ["signed char", "int8", "int8_t", "uchar",
                   "unsigned char", "uint8", "uint8_t", "short",
                   "short int", "signed short", "signed short int",
                   "int16", "int16_t", "ushort", "unsigned short",
                   "unsigned short int", "uint16", "uint16_t", "int",
                   "signed int", "int32", "int32_t", "uint",
                   "unsigned int", "uint32", "uint32_t", "longlong",
                   "long long", "long long int", "signed long long",
                   "signed long long int", "int64", "int64_t",
                   "ulonglong", "unsigned long long",
                   "unsigned long long int", "uint64", "uint64_t",
                   "float", "double", "block"]

    def __init__(self):
        self._load_defaults()

    def __str__(self):
        return self.contents()

    def is_valid(self):
        for key in BaseNrrdHeader.mandatory_keys:
            if self[key] is None:
                return False
        return True

    def write_to_file(self, fileobj):
        try:
            fileobj.write("%s\n" % self.banner)
            fileobj.write("type: %s\n" % self.type)
            fileobj.write("dimension: %s\n" % self.dimension)
            fileobj.write("space: %s\n" % self.space)
            fileobj.write("sizes:")
            for item in self.sizes:
                fileobj.write(" %s" % item)
            fileobj.write("\n")
            fileobj.write("endian: %s\n" % self.endian)
            fileobj.write("encoding: %s\n" % self.encoding)
            fileobj.write("space origin: (%s,%s,%s)\n" %
                          (self.space_origin[0],
                           self.space_origin[1],
                           self.space_origin[2]))
            fileobj.write("space directions:")
            for tpl in self.space_directions:
                if tpl is None:
                    fileobj.write(" none")
                else:
                    fileobj.write(" (%s,%s,%s)" % (tpl[0],tpl[1],tpl[2]))
            fileobj.write("\n")
        except Exception as e:
            raise e

    @property
    def banner(self):
        return self._value_for('banner')

    @banner.setter
    def banner(self, value):
        ##TODO: type checking
        self['banner'] = value
        
    @property
    def type(self):
        return self._value_for('type')
        
    @type.setter
    def type(self, value):
        if value in BaseNrrdHeader.valid_types:
            self['type'] = value
        else:
            raise Exception("%s is not a valid type")

    @property
    def dimension(self):
        return self._value_for('dimension')

    @dimension.setter
    def dimension(self, value):
        try:
            self['dimension'] = int(value)
        except Exception as e:
            raise e

    @property
    def space(self):
        return self._value_for('space')

    @space.setter
    def space(self, value):
        """Can use full name (right-anterior-superior) or
        capitalized abbrv. (RAS)

        """
        abbrv = dict({'RAS':'right-anterior-superior',
                      'RPS':'right-posterior-superior',
                      'LAS':'left-anterior-superior',
                      'LPS':'left-posterior-superior'})
        if value in abbrv.keys():
            value = abbrv[value]
        if value in BaseNrrdHeader.valid_spaces:
            self['space'] = value
        else:
            raise Exception("%s is not a valid space")

    @property
    def sizes(self):
        return self._value_for('sizes')

    @sizes.setter
    def sizes(self, value):
        try:
            newlist = list()
            for item in value:
                newlist.append(item)
            self['sizes'] = newlist
        except Exception as e:
            raise e
        if self.dimension is not None and \
           self.dimension != len(self.sizes):
            dim_provided = len(self.sizes)
            raise Exception("Sizes should contain %s values, equal \
                            to the specified dimension" % dim_provided)

    @property
    def space_directions(self):
        return self._value_for('space_directions')

    @space_directions.setter
    def space_directions(self, value):
        """Expects a list of 3-tuples with a length equal to the
        number of dimensions. Last item may be None in the case of a
        time-series or diffusion volume.

        """
        try:
            for tpl in value:
                if not (tpl is None or len(tpl)==3):
                    raise Exception("Each tuple in list must contain \
                                    3-elements")
            if self.dimension:
                if not self.dimension == len(value):
                    raise Exception("Number of entries should be equal \
                                    to the number of dimensions")
            self['space_directions'] = value
        except Exception as e:
            raise e

    @property
    def endian(self):
        return self._value_for('endian')

    @endian.setter
    def endian(self, value):
        possible_vals = ('little', 'big')
        if value not in possible_vals:
            raise Exception("Value must be one of: %s" % \
                            possible_vals)
        self['endian'] = value

    @property
    def encoding(self):
        return self._value_for('encoding')

    @encoding.setter
    def encoding(self, value):
        possible_values = ["raw", "txt", "text", "ascii", "hex", "gz",
                           "gzip", "bz2", "bzip2"]
        if value not in possible_values:
            raise Exception("Value must be one of the following: %s" % \
                            possible_values)
        self['encoding'] = value

    @property
    def space_origin(self):
        return self._value_for('space_origin')

    @space_origin.setter
    def space_origin(self, value):
        """3-tuple representing space origin in terms of space
        directions.

        """
        ##TODO:
        # need to check to assure validity of space origin, maybe len
        # of space origin should be equal to non-nil space direction
        # entries
        self['space_origin'] = value

    def _load_defaults(self):
        version = 4
        self.banner = "NRRD%04d\n# Complete NRRD file format specification at:\n" \
                      "# http://teem.sourceforge.net/nrrd/format.html" % version
        self.endian = 'little'
        
    def _value_for(self, key):
        """Return value for the given key, or None if key does not
        exist.

        """
        try:
            return self[key]
        except KeyError:
            return None
        except Exception as e:
            raise e
        

class BaseNrrdImage(object):
    """Nrrd image with header and data."""

    @property
    def header(self):
        return self._header

    @property
    def data(self):
        return self._data
        
    def __init__(self, header=BaseNrrdHeader(), data=None):
        self._header = header
        self._data = BaseNrrdData(data)

    def write(self, filename):
        try:
            f = open(filename, 'wb')
            self.header.write_to_file(f)
            f.write('\n')
            self.data.write_to_file(f, self.header.type)
        except Exception as e:
            raise e
        finally:
            f.close()


class DwiNrrdHeader(BaseNrrdHeader):
    """Nrrd DWI header"""

    def __init__(self, bvecs=None, bval=None):
        super(DwiNrrdHeader, self).__init__()
        self.bvecs = bvecs
        self.bval = bval
        self.kinds = ('space', 'space', 'space', 'list')

    def write_to_file(self, fileobj):
        try:
            BaseNrrdHeader.write_to_file(self, fileobj)
            fileobj.write("kinds: %s %s %s %s\n" % (self.kinds[0],
                                                    self.kinds[1],
                                                    self.kinds[2],
                                                    self.kinds[3]))
            fileobj.write("modality:=DWMRI\n")
            fileobj.write("DWMRI_b-value:=%s\n" % self.bval)
            for idx, vec in enumerate(self.bvecs):
                fileobj.write("DWMRI_gradient_%04d:=%s\t%s\t%s\n" %
                              (idx, vec[0], vec[1], vec[2]))
        except Exception as e:
            raise e

    @property
    def bvecs(self):
        return self._value_for('bvecs')

    @bvecs.setter
    def bvecs(self, value):
        """Expects a N x 3 array-like object"""
        ##TODO: validate bvecs on set
        self['bvecs'] = value

    @property
    def bval(self):
        return self._value_for('bval')
        
    @bval.setter
    def bval(self, value):
        ##TODO: validate?
        self['bval'] = int(value)

    @property
    def kinds(self):
        return self._value_for('kinds')

    @kinds.setter
    def kinds(self, value):
        ##TODO: assert number of entries match dimension
        self['kinds'] = value
        
        
class DwiNrrdImage(BaseNrrdImage):
    """Nrrd DWI image."""

    def __init__(self, data=None, bvecs=None, bval=None):
        super(DwiNrrdImage, self).__init__( data=data,
                                            header=DwiNrrdHeader(bvecs=bvecs,
                                                                 bval=bval))
                                 
        
def load_3d_nifti(filename):
    """Read a nifti file and convert to nrrd object.

    Inputs::

    in_filename - nifti file to be converted

    Example
    -------
    >>> import muscip.images.nrrd as nrrd
    >>> nrrd.load_3d_nifti('/path/to/file.nii.gz')
    >>> nrrd.write('/path/to/new/file.nrrd')

    """
    import nibabel
    # get nifti
    nifti = nibabel.load(filename)
    # create base nrrd image
    nrrd = BaseNrrdImage(nifti.get_data())
    nrrd.header.type = str(nifti.get_data_dtype())
    nrrd.header.dimension = 3
    nrrd.header.space = 'RAS'
    nrrd.header.sizes = nifti.shape
    nrrd.header.encoding = 'gzip'
    aff = nifti.get_affine()
    nrrd.header.space_origin = aff[:,3][0:3]
    nrrd.header.space_directions = (aff[0:3,0], aff[0:3,1], aff[0:3,2])
    # return the new nrrd object
    return nrrd

def load_dwi_nifti(filename, bvecs, bval):
    import nibabel
    # get nifti
    nifti = nibabel.load(filename)
    # create base nrrd image
    nrrd = DwiNrrdImage(data=nifti.get_data(), bvecs=bvecs,
                        bval=bval)
    nrrd.header.type = str(nifti.get_data_dtype())
    nrrd.header.dimension = 4
    nrrd.header.space = 'RAS'
    nrrd.header.sizes = nifti.shape
    nrrd.header.encoding = 'gzip'
    aff = nifti.get_affine()
    nrrd.header.space_origin = aff[:,3][0:3]
    nrrd.header.space_directions = (aff[0:3,0], aff[0:3,1], aff[0:3,2], None)
    return nrrd
