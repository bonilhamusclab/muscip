import tables

class TNFibers(object):
    """Container class for fibers"""

    def __init__(self,
                 affine = None,
                 dims = None,
                 h5_file = None,
                 spacing = None,
                 streamlines = None,
                 voxel_size = None
    ):
        if h5_file is None:
            self._create_h5_file()
        else:
            self._load_h5_file(h5_file)
        self.affine = affine
        self.dims = dims
        self.spacing = spacing
        self.voxel_size = voxel_size
        self._write_info()
        if streamlines is not None:
            self.load_streamlines(streamlines)

    def __del__(self):
        try:
            self._info.flush()
            self._points.flush()
            self._streamlines.flush()
            self._hd5.flush()
            self._hd5.close()
            if self._tmpdir is not None:
                from shutil import rmtree
                rmtree(self._tmpdir)
        except Exception, e:
            raise e
            
    @property
    def affine(self):
        try:
            return self._info[0]['affine']
        except:
            return None

    @property
    def dims(self):
        try:
            return self._info[0]['dims']
        except:
            return None

    @property
    def header(self, format='trackvis'):
        try:
            if format=='trackvis':
                import nibabel.trackvis
                header = nibabel.trackvis.empty_header()
                header['vox_to_ras'] = self.affine
                header['dim'] = self.dims
                header['voxel_size'] = self.voxel_size
                return header
            ## format not handled, raise exception
            raise Exception("Format %s is not supported" % format)
        except Exception, e:
            raise e
        
    @property
    def spacing(self):
        try:
            return self._info[0]['spacing']
            self._info.flush()
        except:
            return None

    @property
    def streamlines(self):
        idx = 0
        while idx <= self._points[:]['streamline_idx'].max():
            streamline = []
            for point in self._points.where('streamline_idx == idx'):
                streamline.append(([point['x'],point['y'],point['z']],None,None))
            yield streamline
            idx += 1

    @property
    def voxel_size(self):
        try:
            return self._info[0]['voxel_size']
            self._info.flush()
        except:
            return None

    @affine.setter
    def affine(self, new_affine):
        try:
            if self._info.nrows == 0:
                new_info = self._info.row
                new_info['affine'] =  new_affine
                new_info.append()
            else:
                self._info[0]['affine'] = new_affine
        except Exception, e:
            raise e
        self._flush_data()
        
    @dims.setter
    def dims(self, new_dims):
        try:
            if self._info.nrows == 0:
                new_info = self._info.row
                new_info['dims'] =  new_dims
                new_info.append()
            else:
                self._info[0]['dims'] = new_dims
        except Exception, e:
            raise e
        self._flush_data()

    @spacing.setter
    def spacing(self, new_spacing):
        try:
            if self._info.nrows == 0:
                new_info = {'spacing': new_spacing}
                self._info.append(new_info)
            else:
                self._info[0]['spacing'] = new_spacing
        except Exception, e:
            raise e
        self._flush_data()

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
        
    def length_of_streamline(self, streamline_idx):
        try:
            return self._streamlines[streamline_idx]['length']
        except Exception, e:
            raise e

    def load_streamlines(self, streamlines):
        try:
            self._points, self._streamlines = self._create_ytable()
            for idx, entry in enumerate(streamlines):
                streamline = entry[0]
                streamline_data = dict()
                streamline_data['idx'] = idx
                from numpy import floor
                streamline_data['start_voxel'] = floor(streamline[0])
                streamline_data['end_voxel'] = floor(streamline[-1])
                if self.spacing == 'voxel':
                    streamline_data['length'] = length_of_streamline(streamline, self.voxel_size)
                else:
                    streamline_data['length'] = length_of_streamline(streamline)
                self._streamlines.append(streamline_data)
                for point in streamline:
                    point_data = dict()
                    point_data['streamline_idx'] = idx
                    point_data['x'] = point[0]
                    point_data['y'] = point[1]
                    point_data['z'] = point[0]
                    self._points.append(point_data)
            self._points.flush()
            self._streamlines.flush()
            self._hd5.flush()
        except Exception, e:
            raise e
        
    def write(self, filename, format='trackvis'):
        if format == 'trackvis':
            import nibabel.trackvis
            nibabel.trackvis.write(filename, self.streamlines(),
                                   hdr_mapping=self.header(format='trackvis'),
                                   point_space = self.spacing)

    def _create_h5_file(self):
        try:
            import tempfile
            self._tmpdir = tempfile.mkdtemp()
            from os.path import join
            self._hd5 = tables.openFile(join(self._tmpdir,'fibers.h5'), mode='a')
            self._info = self._hd5.createTable(self._hd5.root, 'info',
                                               _TNFiber_Info)
            self._points = self._hd5.createTable(self._hd5.root,
                                                 'points', _TNFiber_Points)
            self._streamlines = self._hd5.createTable(self._hd5.root,
                                                      'streamlines',
                                                      _TNFiber_Streamlines)
            self._flush_data()
        except Exception, e:
            raise e

    def _flush_data(self):
        self._info.flush()
        self._points.flush()
        self._streamlines.flush()
        self._hd5.flush()

    def _load_h5_file(self, h5_file):
        try:
            self._hd5 = tables.openFile(h5_file, mode='a')
            self._info = self._hd5.getNode('/info')
            self._points = self._hd5.getNode('/points')
            self._streamlines = self._hd5.getNode('/streamlines')
        except Exception, e:
            raise e

    def _write_info(self):
        # if we haven't yet written info
        if self._info.nrows is None:
            # create and add
            info = dict()
            info['affine'] = self.affine
            info['dims'] = self.dims
            info['spacing'] = self.spacing
            info['voxel_size'] = self.voxel_size
            self._info.append(info)
        # if we have existing data in our hd5 file...
        else:
            self._info[0]['affine'] = self.affine
            self._info[0]['dims'] = self.dims
            self._info[0]['spacing'] = self.spacing
            self._info[0]['voxel_size'] = self.voxel_size
        self._flush_data()
        
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
        distance = sqrt( ( (end[0]-start[0]) * vox_dims[0])**2 + \
                         ( (end[1]-start[1]) * vox_dims[1])**2 + \
                         ( (end[2]-start[2]) * vox_dims[2])**2 )
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
    ###########################################################################        
    # TRACKVIS
    if format == 'trackvis':
        import nibabel.trackvis
        data, hdr = nibabel.trackvis.read(file, as_generator=True,
                                          points_space='voxel')
        def streamlines(data):
            for record in data:
                yield record[0]
        return TNFibers(affine = hdr['vox_to_ras'],
                        dims = hdr['dim'],
                        spacing = 'voxel',
                        streamlines = streamlines(data),
                        voxel_size = hdr['voxel_size'])
    ###########################################################################
    # FIBERTOOLS
    # if format == 'fiber_tools':
    #     from scipy.io import loadmat
    #     import nibabel.trackvis
    #     all_data = loadmat(file)
    #     fib_data = all_data['curveSegCell']
    #     aff = all_data['hMatrix']
    #     my_data = dict()
    #     my_header = nibabel.trackvis.empty_header()
    #     nibabel.trackvis.aff_to_hdr(aff, my_header)
    #     for idx, fiber in enumerate(fib_data):
    #         my_data[idx] = {'vertices_vox': fiber[0][:,[1,0,2]]}
    #     ##TODO: review and wrap-up loading
    ###########################################################################            
    # MITK
    # if format == 'mitk':
    #     import numpy as np
    #     # create our new data structures
    #     import nibabel.trackvis        
    #     my_data = dict()
    #     my_header = nibabel.trackvis.empty_header()
    #     # setup for vtk parsing
    #     import vtk                
    #     reader=vtk.vtkPolyDataReader()
    #     reader.SetFileName(file)
    #     reader.Update()
    #     lines = reader.GetOutput().GetLines()
    #     points = reader.GetOutput().GetPoints()
    #     ptr = vtk.vtkIdList()
    #     lines.InitTraversal()
    #     def read_streamlines():
    #         while lines.GetNextCell(ptr):
    #             vertices = list()
    #             for i in range(0,ptr.GetNumberOfIds()):
    #                 vertices.append(points.GetPoint(ptr.GetId(i)))
    #                 yield np.asarray(vertices, dtype=np.float32)
    #     ##TODO: review and wrap-up loading

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

class _TNFiber_Info(tables.IsDescription):
    """Description of internal info storage."""
    affine = tables.Float32Col(shape=(4,4))
    dims = tables.Float32Col(shape=3)
    spacing = tables.StringCol(5)
    voxel_size = tables.Float32Col(shape=3)

class _TNFiber_Points(tables.IsDescription):
    """Description of internal points storage... this forms the basis
    for streamlines.

    """
    streamline_idx = tables.Int32Col(pos=1)
    x = tables.Float32Col(pos=2)
    y = tables.Float32Col(pos=3)
    z = tables.Float32Col(pos=4)

class _TNFiber_Streamlines(tables.IsDescription):
    """Description of internal steramlines storage. This is mainly for
    convenience of grabbing, endpoints, as well as length.

    """
    idx = tables.Int32Col(pos=1)
    start = tables.Float32Col(shape=3,pos=2)
    end = tables.Float32Col(shape=3,pos=3)
    length = tables.Float32Col(pos=4)
