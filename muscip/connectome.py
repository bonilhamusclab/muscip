import networkx


class TNConnectome(networkx.Graph):
    """Connectome as defined by me!"""

    def __init__(self):
        networkx.Graph.__init__(self)

    def edge_observations_for_key(self, key, include_info=True, node_name_key=None):
        """Return record for a given key in the form of a dictionary.

        Inputs::

          key: key for which to export - key must contain no more than
               a single value for each edge

          include_info: (Boolean) should export include info

          node_name_key: if provided, will use information stored in
                         the node key, to label edges, if not
                         provided, will use node labels - edges will
                         take on the form of nodeA-nodeB

        """
        try:
            if include_info:
                try:
                    from copy import copy
                    data = copy(self.graph['info'])
                except KeyError:
                    print "No global info found for connectome."
                    data = dict()                    
                except Exception, e:
                    print e
            else:
                data = dict()
                # for every edge, get data and put in result
            for a,b in self.edges_iter():
                if node_name_key is not None:
                    # get node names
                    nodeA = self.node[a][node_name_key]
                    nodeB = self.node[b][node_name_key]
                else:
                    nodeA = a
                    nodeB = b
                data["%s-%s" % (str(nodeA), str(nodeB))] = self[a][b][key]
            return data
        except Exception, e:
            print e
            return None

    def node_observations_for_key(self, key, include_info=True, node_name_key=None):
        """Return record for a given key in the form of a dictionary.

        Inputs::

          key: key for which to export - key must contain no more than
               a single value for each edge

          include_info: (Boolean) should export include info

          node_name_key: if provided, will use information stored in
                         the node key, to label nodes, if not
                         provided, will use node labels
        
        """
        try:
            if include_info:
                try:
                    from copy import copy
                    data = copy(self.graph['info'])
                except KeyError:
                    print "No global info found for connectome."
                    data = dict()
                except Exception, e:
                    print e
            else:
                data = dict()
                # for every edge, get data and put in result
            for a in self.nodes_iter():
                if node_name_key is not None:
                    # get node names
                    node_name = self.node[a][node_name_key]
                else:
                    node_name = a
                data[node_name] = self.node[a][key]
            return data
        except Exception, e:
            print e
            return None
        
    def get_info(self):
        """Return connectome info if any."""
        try:
            return self.graph['info']
        except KeyError:
            print "Connectome contains no associated info."
            return None
        except Exception, e:
            print e
            return None
            
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
        
    def set_info(self, info):
        """Set info for connectome. This will add a info dictionary
        and provided keys/values to the connectome.

        Inputs::

          info: a dictionary containing info to be stored

        """
        self.graph['info'] = info

    def write(self, filename):
        """Write connectom to the given filename as gpickle.

        Inputs::

          filename - path at which to save output
        
        """
        try:
            networkx.write_gpickle(self, filename)
        except Exception, e:
            print e
        
    def write_fibers(self, filename):
        """Write fibers as trackvis file to the give filename.

        Inputs::

          filename - path at which to save output
        
        """
        try:
            streamlines = []
            for i,j in self.edges_iter():
                for streamline in self[i][j]['streamlines']:
                    streamlines.append(streamline)
            from nibabel.trackvis import TrackvisFile
            TrackvisFile(streamlines).to_file(filename)
        except Exception, e:
            print e

# MODULE LEVEL FUNCTIONS
def edge_data_for_connectomes_as_df(C_list, edge_data_key, filename,
                                    include_info=True, node_name_key=None):
    """Write data frame as csv file to the given path.

    Inputs::

      C_list: list of connectome paths
      edge_data_key: extract values belonging to this key
      filename: write file to this path
      include_info: (Boolean) include info?; default=True
      node_name_key: if provided, use this key to extract node names
    
    """
    # get the union of all keys, we need this to be all inclusive
    all_info_keys = list()
    all_data_keys = list()
    for C_path in C_list:
        C = read_gpickle(C_path)
        if include_info:
            # add any info keys not yet in union
            try:
                for k in C.graph['info'].keys():
                    if k not in all_info_keys:
                        all_info_keys.append(k)
            except KeyError:
                print "Info requested, but no info exists in Connectome."
            except Exception, e:
                print e
        # add every edge key we find
        try:
            c_keys = C.edge_observations_for_key(edge_data_key,
                                                 include_info=False,
                                                 node_name_key=node_name_key).keys()
            for k in c_keys:
                if k not in all_data_keys:
                    all_data_keys.append(k)
        except Exception, e:
            print e
    all_info_keys.sort()
    all_data_keys.sort()
    all_keys = all_info_keys + all_data_keys
    print all_keys
    # export to csv file ~~~~~~~~~~~~~~~~~~~
    try:
        import csv
        import os.path as op
        fout = open(op.abspath(filename), 'wt')
        writer = csv.DictWriter(fout, fieldnames=all_keys)
        header = dict( (k,k) for k in all_keys)
        writer.writerow(header)
        for C_path in C_list:
            C = read_gpickle(C_path)
            data = C.edge_observations_for_key(edge_data_key,
                                               include_info=include_info,
                                               node_name_key=node_name_key)
            writer.writerow(data)
    except Exception, e:
        print e
    finally:
        fout.close()
    
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
    # get our vertex key - we will use voxel spacing to determine
    # intersection of ROIs, we are assuming that ROI atlas is in same
    # space as diffusion, as is best practice
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

def read_gpickle(filename):
    """Create connectome object from a given gpickle.

    Input::

      filename - path to networkx gpickle file
    """
    read_graph = networkx.read_gpickle(filename)
    connectome = TNConnectome()
    connectome.graph = read_graph.graph
    for i,data in read_graph.nodes_iter(data=True):
        connectome.add_node(i,data)
    for i,j,data in read_graph.edges_iter(data=True):
        connectome.add_edge(i,j,data)
    return connectome
