
class TNFibers(object):
    """Container class for fibers"""

    def __init__(self,
                 format = 'trackvis',
                 header = None,
                 spacing = 'mm',
                 streamlines = None
    ):
        self.format = format
        self.set_header(header)
        self.set_spacing(spacing)
        # self._load_alternative_spacings()
        self._streamlines = streamlines
        
    def get_data(self):
        if self.format == 'trackvis':
            return self.data

    def get_dim(self):
        return self.get_value_for_header_key('dim')
    
    def get_header(self):
        return self.header
        
    def get_spacing(self):
        return self._spacing
        
    def get_voxel_size(self):
        if self.format == 'trackvis':
            if self.header is not None:
                return self.get_value_for_header_key('voxel_size')

    def get_value_for_header_key(self, key):
        if not self.header:
            return None
        if self.format == 'trackvis':
            try:
                return self.header[key]
            except:
                return None
        return None

    @property
    def number_of_streamlines(self):
        return self.get_value_for_header_key('n_count')
                
    def set_dims(self,x,y,z):
        self.set_value_for_header_key('dim',[x,y,z])
        
    def set_header(self, header):
        if not header:
            self.header = None
            return
        if hasattr(header, 'flags.writeable'):
            header.flags.writeable = True
        self.header = header

    def set_spacing(self, new_spacing):
        accepted_spacing = ['mm', 'vox']
        if new_spacing not in accepted_spacing:
            print "%s is not an accepted spacing" % new_spacing
        if new_spacing == 'mm':
            if self.format == 'trackvis':
                self._vertex_key = 'vertices_mm'
        if new_spacing == 'vox':
            if self.format == 'trackvis':
                self._vertex_key = 'vertices_vox'
        self._spacing = new_spacing

    def set_value_for_header_key(self, key, value):
        if not self.header:
            return None
        if self.format == 'trackvis':
            try:
                self.header[key] = value
            except Exception as e:
                raise e
            
    def set_voxel_size(self,x,y,z):
        self.set_value_for_header_key('voxel_size',[x,y,z])

    @property
    def streamlines(self):
        try:
            self._streamlines
            return self._streamlines
        except AttributeError:
            return None
            
    def write(self, filename):
        if self.format == 'trackvis':
            import nibabel.trackvis
            if self.get_spacing == 'vox':
                ps = 'voxel'
            else:
                ps = None
            if self.header is not None:
                nibabel.trackvis.write(filename, self.streamlines(),
                                       hdr_mapping=self.header, points_space=ps)
            else:
                nibabel.trackvis.write(filename, self.streamlines(), point_space=ps)
                
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
                

        
################################################################################        
# Module Level Functions -------------------------------------------------------
################################################################################

def cmat_for_key(connectome, key, number_of_nodes=None,
                 force_symmetric=True):
    """Return a N x N connection matrix for given connectome and
    key. The connection matrix is returned as a numpy ndarray.

    """

    # create our new shiny connection matrix
    import numpy
    if number_of_nodes is None:
        n = max(connectome.nodes())
    else:
        n = number_of_nodes
    new_cmat = numpy.zeros((n,n))

    # extract the value for key for every edge in the given connectome
    for i,j in connectome.edges_iter():
        new_cmat[i-1][j-1] = connectome[i][j][key]
        
    # do we need to do anything regarding symmetry?
    if force_symmetric and (new_cmat - new_cmat.T != 0).any():
        #...if one-sided (no information below diagonal)
        if (numpy.tril(new_cmat,-1) == 0).all():
            # project above diagonal onto below diagonal
            new_cmat += numpy.tril(new_cmat.T, -1)
        #...else, we will assume two-sided unequal
        else:
            # our solution will be to take the mean of each pair of
            # reflected indices
            new_cmat = (new_cmat + new_cmat.T ) / 2.0

    # return the cmat
    return new_cmat
        
    
def extract_hagmann_density(connectome, roi_img, wm_img):
    """Populate hagmann density for given connectome"""

    def inverse_sum(elements):
        inverse_elements = []
        for element in elements:
            inverse_elements.append( 1.0 / element )
        from math import fsum
        return fsum(inverse_elements)

    # get surface areas for ROIs
    import images
    surface_area = images.surface_area_for_rois(roi_img, wm_img)

    # for every edge...
    for i,j in connectome.edges_iter():
        calc_hd = ( ( 2.0 / ( surface_area[i] + surface_area[j] ) ) * \
                  inverse_sum( connectome[i][j]['streamlines_length'] ) )
        connectome[i][j]['hagmann_density'] = calc_hd
    
    
def extract_scalars(fibers, connectome, scalar_img, scalar_name, scale_factor=[1.,1.,1.]):
    """Extract scalar values for given fiber, img combination, and add
    to connectome.

    Return copy of updated connectome.

    """
    import numpy # we will need this for mean, std over arrays

    # get scalar img data
    scalar_data = scalar_img.get_data()

    # get vertex key
    if fibers.get_spacing() == 'mm':
        vertex_key = 'vertices_mm'
    elif fibers.get_spacing() == 'vox':
        vertex_key = 'vertices_vox'

    # generate out keys
    mean_key = "%s_mean" % scalar_name
    std_key = "%s_std" % scalar_name

    # record image path in connectome
    from os.path import abspath
    connectome.graph['%s_img' % scalar_name] = abspath(scalar_img.get_filename())

    # get streamlines
    streamlines = fibers.get_data()
    # for every edge in our connectome...
    for i,j in connectome.edges_iter():
        # start a new collection
        collected_values = []
        # for every streamline belonging to edge...
        for streamline in connectome[i][j]['streamlines']:
            # for every vertex belonging to streamline
            for vertex in streamlines[streamline][vertex_key]:
                voxel = [int(vertex[0] / scale_factor[0]),
                         int(vertex[1] / scale_factor[1]),
                         int(vertex[2] / scale_factor[2])]
                # try to get value for voxel (it is possible it is out
                # of bounds)
                try:
                    value = scalar_data[tuple(voxel)]
                except IndexError:
                    continue
                collected_values.append(value)

        # calculate aggregate values and write to connectom
        collected_values_as_array = numpy.asarray(collected_values)
        connectome[i][j][mean_key] = collected_values_as_array.mean()
        connectome[i][j][std_key] = collected_values_as_array.std()

    # return updated connectome
    return connectome
    
def generate_connectome(fibers, roi_img):
    """Return connectome as a networkx object"""
    import networkx
    # create empty graph
    connectome = networkx.Graph()
    from os.path import abspath
    connectome.graph['roi_img'] = abspath(roi_img.get_filename())
    # get our ROI data
    roi_data = roi_img.get_data()
    # get our vertex key
    if fibers.get_spacing() == 'mm':
        vertex_key = 'vertices_mm'
    elif fibers.get_spacing() == 'vox':
        vertex_key = 'vertices_vox'
    # TODO: should write an interater so that we don't need to load
    # all of this into memory
    streamlines = fibers.get_data()
    for streamline in streamlines:
        # get the endpoints of streamline in terms of voxel indices
        vertex_i = streamlines[streamline][vertex_key][0]
        vertex_j = streamlines[streamline][vertex_key][-1]
        voxel_i = [int(vertex_i[0]), int(vertex_i[1]),
                   int(vertex_i[2])]
        voxel_j = [int(vertex_j[0]), int(vertex_j[1]),
                   int(vertex_j[2])]
        # try and get value for voxel indices from roi data (it is
        # possible that we are outside of the bounds of the ROI, due
        # to the propegation of the tracks beyond the bounds of ROI
        # image in tracking)
        try:
            label_i = roi_data[tuple(voxel_i)]
            label_j = roi_data[tuple(voxel_j)]
        except IndexError:
            continue
        # if both endpoints reside in ROIs (both endpoints are
        # non-zero)...
        if label_i != 0 and label_j != 0:
            # ...then we need to add to the fiber count for the
            # specified edge
            try:
                connectome[label_i][label_j]['number_of_fibers'] += 1
                connectome[label_i][label_j]['streamlines'].append(streamline)
                connectome[label_i][label_j]['streamlines_length'].append(
                    length_of_streamline(streamlines[streamline]))
            # handle the case where the edge does not yet exist
            except KeyError:
                connectome.add_edge(label_i, label_j)
                connectome[label_i][label_j]['number_of_fibers'] = 1
                connectome[label_i][label_j]['streamlines'] = [streamline]
                connectome[label_i][label_j]['streamlines_length'] = [length_of_streamline(streamlines[streamline])]

    # calculate and store mean fiber lengths and std
    import numpy
    for i,j in connectome.edges_iter():
        connectome[i][j]['fiber_length_mean'] = numpy.asarray(connectome[i][j]['streamlines_length']).mean()
        connectome[i][j]['fiber_length_std'] = numpy.asarray(connectome[i][j]['streamlines_length']).std()

    # return our results
    return connectome
    
def length_of_streamline(streamline, vox_dims=[1.,1.,1.]):
    """Return the length of streamline as determined by it's
    vertices

    '"""
    segments = []
    idx = 0
    vertices = streamline
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
        data, hdr = nibabel.trackvis.read(file, as_generator=True, points_space='voxel')
        def streamlines(data):
            for record in data:
                yield record[0]
        return TNFibers(format = 'trackvis',
                        header = hdr,
                        spacing = 'mm',
                        streamlines = streamlines(data))
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
        streamline_idx = 0
        # setup for vtk parsing
        import vtk                
        reader=vtk.vtkPolyDataReader()
        reader.SetFileName(file)
        reader.Update()
        lines = reader.GetOutput().GetLines()
        points = reader.GetOutput().GetPoints()
        ptr = vtk.vtkIdList()
        lines.InitTraversal()
        while lines.GetNextCell(ptr):
            vertices = list()
            for i in range(0,ptr.GetNumberOfIds()):
                vertices.append(points.GetPoint(ptr.GetId(i)))
            my_data[streamline_idx] = {'vertices_vox':
                                       np.asarray(vertices, dtype=np.float32)}
            streamline_idx += 1
        return TNFibers(data=my_data, format='trackvis', header=my_header, spacing='vox')

def transform_streamline_by_aff(streamline, aff):
    import numpy
    xStreamline = []
    for idx in streamline:
        xIdx = numpy.dot(aff, numpy.append(idx , 1))[0:3]
        xStreamline.append(xIdx)
    return xStreamline
    
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
