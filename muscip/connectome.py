import networkx

class TNConnectome(networkx.Graph):
    """Connectome as defined by me!"""

    def __init__(self,
                 filename=None
    ):
        networkx.Graph.__init__(self)

    def populate_hagmann_density(self, WM_img, ROI_img):
        """Populate hagmann density."""

        def inverse_sum(elements):
            inverse_elements = []
            for element in elements:
                inverse_elements.append( 1.0 / element )
            from math import fsum
            return fsum(inverse_elements)

        # get surface areas for ROIs
        import images
        surface_area = images.surface_area_for_rois(ROI_img, WM_img)

        # for every edge...
        for i,j in self.edges_iter():
            calc_hd = ( ( 2.0 / ( surface_area[i] + surface_area[j] ) ) * \
                        inverse_sum( self[i][j]['streamlines_length'] ) )
            self[i][j]['hagmann_density'] = calc_hd

    def populate_node_info(self, ROI_img_data, node_info_file):
        """Populate node information such as ROI name and center of
        mass in the given connectome.

        Inputs::
        
          ROI_img - the ROI image data (as numpy ndarray) from which
                    the connectome was generated. We will use this to
                    generate centers of mass.
          node_info_file - a graphml file which contains node
                           information

        """
        import numpy as np
        # read node info file
        try:
            node_info = networkx.read_graphml(node_info_file)
        except:
            print "Could not read node_info_file!"
            return
        # add node information
        for label, data in node_info.nodes_iter(data=True):
            self.add_node(int(label), data)
            self.node[int(label)]['subject_position'] = tuple( np.mean ( np.where ( ROI_img_data == int(label) ), axis=1 ) )

# MODULE LEVEL FUNCTIONS
def generate_connectome(fib, roi_img, node_info=None):
    """Return a TNConnectome object

    Example
    -------
    import tn_image_processing as tnip
    import nibabel
    fibers = tnip.fibers.read('path/to/track.trk')
    roi = nibabel.load('path/to/roi.nii.gz')
    C = tnip.connectome.generate_connectome(fibers, roi)

    Input::

      [mandatory]
      fibers - a loaded fiber object from
               tn_image_processing.fibers
      roi_img - a loaded nibabel image

    """
    import os, fibers

    # create connectome object and store initializing properties
    connectome = TNConnectome()
    connectome.graph['roi_img'] = os.path.abspath(roi_img.get_filename())
    # get ROI data
    roi_data = roi_img.get_data()
    # load node info if provided
    if node_info:
        connectome.populate_node_info(roi_data, node_info)
    # get our vertex key
    if fib.get_spacing() == 'mm':
        vertex_key = 'vertices_mm'
    elif fib.get_spacing() == 'vox':
        vertex_key = 'vertices_vox'
    # TODO: should write an interater so that we don't need to load
    # all of this into memory
    streamlines = fib.get_data()
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
                    fibers.length_of_streamline(streamlines[streamline], vox_dims=fib.get_voxel_size()))
            # handle the case where the edge does not yet exist
            except KeyError:
                connectome.add_edge(label_i, label_j)
                connectome[label_i][label_j]['number_of_fibers'] = 1
                connectome[label_i][label_j]['streamlines'] = [streamline]
                connectome[label_i][label_j]['streamlines_length'] = \
                    [fibers.length_of_streamline(streamlines[streamline], vox_dims=fib.get_voxel_size())]

    # calculate and store mean fiber lengths and std
    import numpy
    for i,j in connectome.edges_iter():
        connectome[i][j]['fiber_length_mean'] = numpy.asarray(connectome[i][j]['streamlines_length']).mean()
        connectome[i][j]['fiber_length_std'] = numpy.asarray(connectome[i][j]['streamlines_length']).std()

    # return our results
    return connectome
    