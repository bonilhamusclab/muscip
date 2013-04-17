
class TNFibers(object):
    """Container class for fibers"""

    def __init__(self, **kwargs):
        self._filename = kwargs.get('filename', None)
        self._fibers = kwargs.get('fibers', None)

    @property
    def fibers(self):
        raise Exception("Method <fibers> is not implemented by \
                        subclass: %s" % type(self))

    @property
    def filename(self):
        try:
            return self._filename
        except AttributeError:
            return None
            
    @property
    def format(self):
        try:
            return self._format
        except:
            raise Exception("Subclass: %s does not properly encode \
                            format" % type(self))
            
    @property
    def number_of_fibers(self):
        raise Exception("Method <number_of_fibers> is not implemented by \
                        sublcass: %s" % type(self))
    @property
    def shape(self):
        raise Exception("Method <shape> is not implemented by \
                        subclass: %s" % type(self))

    @property
    def spacing(self):
        raise Exception("Method <spacing> is not implemented by \
                        subclass: %s" % type(self))

    @property
    def voxel_size(self):
        raise Exception("Method <voxel_size> is not implemented by \
                        subclass: %s" % type(self))
    
    def write(self, filename):
        raise Exception("Method <write> is not implemented by \
                        subclass: %s" % type(self))
        
################################################################################        
# Module Level Functions -------------------------------------------------------
################################################################################
def fiber_length(fiber, vox_dims=[1.,1.,1.]):
    """Return the length of streamline as determined by it's
    vertices

    '"""
    segments = []
    idx = 0
    vertices = fiber
    while idx < len(vertices) - 1:
        start = vertices[idx]
        end = vertices[idx+1]
        from math import sqrt
        distance = sqrt( ( (end[0]-start[0]) / vox_dims[0])**2 + \
                         ( (end[1]-start[1]) / vox_dims[1])**2 + \
                         ( (end[2]-start[2]) / vox_dims[2])**2 )
        segments.append(distance)
        idx += 1
    from math import fsum
    return fsum(segments)

def transform_fiber_by_aff(fiber, aff):
    import numpy
    xFiber = []
    for idx in fiber:
        xIdx = numpy.dot(aff, numpy.append(idx , 1))[0:3]
        xFiber.append(xIdx)
    return xFiber
    
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
