
class TNFibers(object):
    """Container class for fibers"""

    def __init__(self,
                 data = None,
                 format = 'trackvis',
                 header = None,
                 spacing = 'mm'
    ):
        self.data = data
        self.format = format
        self.header = header
        self._load_alternative_spacings()
        
    def get_data(self):
        return self.data

    def get_dim(self):
        if self.format == 'trackvis':
            if self.header:
                return self.get_value_for_header_key('dim')
        
    def get_header(self):
        return self.header
        
    def get_data(self):
        if self.format == 'trackvis':
            return self.data

    def get_voxel_size(self):
        if self.format == 'trackvis':
            if self.header is not None:
                return self.get_value_for_header_key('voxel_size')

    def get_value_for_header_key(self, key):
        if self.format == 'trackvis':
            try:
                return self.header[self._trackvis_header_fields()[key]]
            except:
                return None

    def _load_alternative_spacings(self):
        """Populate the data with both voxel and mm based spacing
        values.

        The outcome is that within data, there will be keys for
        both.

        """
        voxel_size = self.get_voxel_size()
        if voxel_size is not None:
            if not self.data[0].has_key('vertices_vox'):
                for idx in self.data:
                    self.data[idx]['vertices_vox'] = self.data[idx]['vertices_mm'] / voxel_size
            if not self.data[0].has_key('vertices_mm'):
                for idx in self.data:
                    self.data[idx]['vertices_mm'] = self.data[idx]['vertices_vox'] * voxel_size
        
    def _trackvis_header_fields(self):
        return {
            'id_string': 0,
            'dim': 1,
            'voxel_size': 2,
            'origin': 3,
            'n_scalars': 4,
            'scalar_name': 5,
            'n_properties': 6,
            'property_name': 7,
            'vox_to_ras': 8,
            'reserved': 9,
            'voxel_order': 10,
            'pad2': 11,
            'image_orientation_patient': 12,
            'pad1': 13,
            'invert_x': 14,
            'invert_y': 15,
            'invert_z': 16,
            'swap_xy': 17,
            'sway_yz': 18,
            'swap_zx': 19,
            'n_count': 20,
            'version': 21,
            'hdr_size': 22
        }
                
        
# Module Level Functions
def sum_of_inverse_fiber_lengths(roi_img, endpoints, lengths):
    """For a given ROI atlas and set of fibers, return an adjacency
    matrix where each component(i,j) holds the sum of inverse fiber
    lengths for all fibers running between region(i) and region(j).

    Inputs:

    roi_img - a nibabel image, or subclass of nibabel image, such as a
    TNImage, that contains the ROIs encoded with integer labels

    endpoints - a numpy ndarray that contains the two, xyz endpoints
    for each fiber

    lengths - a numpy ndarray that contains a length value for each
    fiber (must be ordered same as endpoint array so that endpoint)

    """
    import numpy as np
    roi_data = roi_img.get_data()
    n = len(np.unique(roi_data))
    result = np.zeros((n,n))
    for idx, ends in enumerate(endpoints):
        endpoint_i, endpoint_j = ends
        ix,iy,iz = (endpoint_i[0],endpoint_i[1],endpoint_i[2])
        jx,jy,jz = (endpoint_j[0],endpoint_j[1],endpoint_j[2])
        try:
            i = roi_data[ix,iy,iz]
            j = roi_data[jx,jy,jz]
            l = lengths[idx]
            if i > 0 and j > 0 and l > 0:
                try:
                    result[i][j] += 1.0 / l
                except Exception as e:
                    print e
                    raise
        except IndexError:
            print "Warning: one of the following endpoints is out of " + \
            "range: %s - %s" % (endpoint_i, endpoint_j)
        except Exception as e:
            print e
            raise
    return result

def read(file, format='trackvis'):
    """Load a TNFibers object from file and return"""
    accepted_formats = ['trackvis']
    if format not in accepted_formats:
        print "%s is not an accepted format" % format
        raise
    if format == 'trackvis':
        import nibabel.trackvis
        import numpy
        read_data = nibabel.trackvis.read(file)
        my_data = dict()
        for idx, rec in enumerate(read_data[0]):
            my_data[idx] = {'vertices_mm': rec[0]}
        my_header = read_data[1].tolist()
        return TNFibers(data=my_data, format='trackvis', header=my_header)
