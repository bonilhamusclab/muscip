from table import *

class TNFibers(object):
    """Container class for fibers"""

    def __init__(self,
                 affine = None,
                 dims = None,
                 spacing = None,
                 streamlines = None,
                 voxel_size = None
    ):
        # create an hdf5 table in a temp directory
        self.set_affine(affine)
        self.set_dims(dims)
        self.set_spacing(spacing)
        self.set_voxel_size(voxel_size)
        self.load_streamlines(streamlines)

    @property
    def affine(self):
        try:
            return self._affine
        except:
            return None

    @property
    def dims(self):
        try:
            return self._dims
        except:
            return None

    @property
    def spacing(self):
        try:
            return self._spacing
        except:
            return None

    @property
    def streamlines(self):
        ##TODO: implement streamlines generator
        return None

    @property
    def voxel_size(self):
        try:
            return self._voxel_size
        except:
            return None

    @affine.setter
    def affine(self, new_affine):
        try:
            from numpy import asarray
            new_value = asarray(new_affine)
            if new_value.shape == (4,4):
                self._affine = new_value
            else:
                raise Exception("Affine should be a 4x4 matrix, but instead was %s" % new_value.shape)
        except Exception, e:
            raise e

    @dims.setter
    def dims(self, new_dims):
        try:
            from numpy import asarray
            new_value = asarray(new_dims)
            if new_value.shape == (1,3):
                self._dims = new_value
            else:
                raise Exception("Dims should be a 1x3 matrix, but instead was %s" % new_value.shape)
        except Exception, e:
            raise e

    @spacing.setter
    def spacing(self, new_spacing):
        try:
            accepted_spacings = ['voxel', 'voxelmm', 'rasmm']
            if new_spacing in accepted_spacings:
                self._spacing = new_spacing
            else:
                raise Exception("Valid spacing options include: %s, but %s was provided."
                                % (accepted_spacings, new_spacing))
        except Exception, e:
            raise e

    @voxel_size.setter
    def voxel_size(self, new_voxel_size):
        try:
            from numpy import asmatrix
            new_value = asmatrix(new_voxel_size)
            if new_value.shape == (1,1):
                self._voxel_size = new_value * [1,1,1]
            if new_value.shape == (1,3):
                self._voxel_size = new_value
            else:
                raise Exception("Voxel size must be provided as a 1x1 or 1x3 matrix, not %s." %
                                new_value.shape)
        except Exception, e:
            raise e
        
    def length_of_streamline(self, streamline_key):
        ##TODO: impolement length of streamline
        pass
        
    def write(self, filename, format='trackvis'):
        if format == 'trackvis':
            ##TODO: implement trackvis write
            pass
        if format == 'fib'
           ##TODO: implement fib write
                
    def _trackvis_header(self):
        

        
################################################################################        
# Module Level Functions -------------------------------------------------------
################################################################################
    
def length_of_streamline(streamline, vox_dims=[1.,1.,1.]):
    """Return the length of streamline as determined by it's
    vertices

    '"""
    segments = []
    idx = 0
    vertices = streamline['vertices_mm']
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

def read(file, format='trackvis'):
    """Load a TNFibers object from file and return"""
    accepted_formats = ['trackvis','fiber_tools','mitk']
    if format not in accepted_formats:
        print "%s is not an accepted format" % format
        raise
    # TRACKVIS
    if format == 'trackvis':
        import nibabel.trackvis
        read_data = nibabel.trackvis.read(file)
        my_data = dict()
        for idx, rec in enumerate(read_data[0]):
            my_data[idx] = {'vertices_mm': rec[0]}
        return TNFibers(data=my_data, format='trackvis', header=read_data[1])
    # FIBERTOOLS
    if format == 'fiber_tools':
        from scipy.io import loadmat
        import nibabel.trackvis
        all_data = loadmat(file)
        fib_data = all_data['curveSegCell']
        aff = all_data['hMatrix']
        my_data = dict()
        my_header = nibabel.trackvis.empty_header()
        nibabel.trackvis.aff_to_hdr(aff, my_header)
        for idx, fiber in enumerate(fib_data):
            my_data[idx] = {'vertices_vox': fiber[0][:,[1,0,2]]}
        return TNFibers(data=my_data, format='trackvis', header=my_header, spacing='vox')
    # MITK
    if format == 'mitk':
        import numpy as np
        # create our new data structures
        import nibabel.trackvis        
        my_data = dict()
        my_header = nibabel.trackvis.empty_header()
        # setup for vtk parsing
        import vtk                
        reader=vtk.vtkPolyDataReader()
        reader.SetFileName(file)
        reader.Update()
        lines = reader.GetOutput().GetLines()
        points = reader.GetOutput().GetPoints()
        ptr = vtk.vtkIdList()
        lines.InitTraversal()
        def read_streamlines():
            while lines.GetNextCell(ptr):
                vertices = list()
                for i in range(0,ptr.GetNumberOfIds()):
                    vertices.append(points.GetPoint(ptr.GetId(i)))
                    yield np.asarray(vertices, dtype=np.float32)
        return TNFibers(streamlines=read_streamines(), header=my_header, spacing='vox')

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
