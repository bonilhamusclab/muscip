import datetime, networkx, os, pickle, shutil
from ..images import TNImage

__CLINICAL_INFO_FILENAME__ = 'clinical_info.pickle'
__MANIFEST_FILENAME__ = 'manifest.pickle'
__MODIFICATIONS_FILENAME__ = 'modifications.pickle'
__NETWORK_FILENAME__ = 'graph.gpickle'
__ROI_IMAGE_PREFIX__ = 'roi'
__VERSION__ = 2.0

class TNConnectome(object):
    """Connectome as defined by me!"""

    def __init__(self,
                 clinical_info = None,
                 filename = None,
                 network = None,
                 roi_image = None
    ):
        if filename is not None:
            self._filename = filename
        else:
            self._filename = None
            self.created = datetime.datetime.now()
            self.version = __VERSION__
        if clinical_info is not None:
            self.clinical_info = clinical_info
        if network is not None:
            self.network = network
        if roi_image is not None:
            self.roi_image = roi_image

    @property
    def created(self):
        """Return date of connectome creation."""
        try:
            return self.manifest['CREATED']
        except KeyError:
            print "Connectome contains no creation date."
            return None
        except Exception, e:
            raise e

    @created.setter
    def created(self, timestamp):
        if isinstance(timestamp, datetime.datetime):
            self.manifest['CREATED'] = timestamp
        else:
            raise Exception("Expected value of type: datetime.datetime")
            
    @property
    def clinical_info(self):
        """Return clinical info, if any, as dictionary."""
        # 1) clinical info is an object in memory                    
        try:
            return self._clinical_info
        except AttributeError:
            # 2) clinical info exists as a file at std path
            if self._filename is not None:
                clinfo_path = self.clinical_info_path
                if os.path.exists(clinfo_path):
                    try:
                        clinfo_file = open(clinfo_path, 'rb')
                        self._clinical_info = pickle.load(clinfo_file)
                        clinfo_file.close()
                    except:
                        raise Exception('Could not parse clinical info file.')
                # 3) we have a connectome directory in which no
                # clinical info exists, we need to generate and return
                # an empty dictionary
                else:
                    self._clinical_info = dict()
            # 4) we have no connectome directory, thus we have no
            # clinical info on disk
            else:
                self._clinical_info = dict()
            # return
            return self._clinical_info

    @clinical_info.setter
    def clinical_info(self, new_info):
        """Set info for connectome. This will add a info dictionary
        and provided keys/values to the connectome.

        Inputs::

          new_info: a dictionary, or filename pointing to pickled
                    dictionary, containing info to be stored

        """
        if isinstance(new_info, basestring):
            try:
                f = open(new_info, 'rb')
                new_info = pickle.load(new_info)
            finally:
                f.close()
        if isinstance(new_info, dict):
            self._clinical_info = new_info
        else:
            raise Exception('New info must be provided in the form of \
                            a dictionary, or filename of pickled dictionary.')

    def generate_connectome(self):
        """Place-holder method. This should be implemented by each
        sub-class.

        """
        raise Exception('Method: <generate_connectome> should be \
                        implemented by subclass but it is not.')
        
    @property
    def manifest(self):
        """Return the manifest if one exists. If it does not exist
        create an empty manifest and return.

        """
        try:
            # 1) existing manifest object in memory
            return self._manifest
        except AttributeError:
            # 2) manifest file stored on disk within connectome file
            if self._filename is not None:
                manifest_path = self.manifest_file_path
                if os.path.exists(manifest_path):
                    try:
                        manifest_file = open(manifest_path, 'rb')
                        self._manifest = pickle.load(manifest_file)
                    except:
                        raise Exception('Could not parse manifest file.')
                    finally:
                        manifest_file.close()
            # 3) no existing manifest (disk|mem), must create new
            else:
                self._manifest = dict()
            # return
            return self._manifest

    @property
    def manifest_file_path(self):
        """Return standard path to manifest file if one exists, else
        return None. Path will not exist unless filename is set.

        """
        if self._filename is not None:
            return os.path.join(self._filename, __MANIFEST_FILENAME__)
        else:
            return None
            
    @property
    def modifications(self):
        """Return list of modifications as a list of dictionaries w/
        the following fields: { timestamp | user | desc }.

        """
        # existing modifications object in memory
        try:
            return self._modifications
        except AttributeError:
            # modifications file stored on disk within connectome file
            if self._filename is not None:
                modpath = self.modifications_file_path
                if os.path.exists(modpath):
                    try:
                        modifications_file = open(modpath, 'rb')
                        self._modifications = pickle.load(modifications_file)
                    except:
                        raise Exception('Could not parse modifications file.')
                    finally:
                        modifications_file.close()
            # generate empty modifications list
            else:
                self._modifications = []
            # return our generated modifications object
            return self._modifications

    @property
    def modifications_file_path(self):
        """Return standard path to modifications file if one exists,
        else return None. Path will not exist unless filename is
        set.

        """
        if self._filename is not None:
            return os.path.join(self._filename,__MODIFICATIONS_FILENAME__)
        else:
            return None

    @property
    def network(self):
        """Return network as a networkx graph."""
        try:
            return self._network
        except AttributeError:
            try:
                self._network = networkx.read_gpickle(self.network_path)
            except:
                self._network = networkx.Graph()
            return self._network

    @network.setter
    def network(self, new_network):
        """Set the connectome's network to new value. This might be
        useful for doing something like null modeling, or assigning
        connectomes that were generated through some other
        interface.

        Input::

          new_network: this can be either a networkx.Graph object, or
                       a file name pointing to a pickled
                       networkx.Graph object.

        """
        if isinstance(new_network, basestring):
            try:
                new_network = networkx.read_gpickle(new_network)
            except:
                raise Exception("Could not load %s" % new_network)
        if isinstance(new_network, networkx.Graph):
            self._network = new_network
        else:
            raise Exception("Could not load %s" % new_network)
            
    @property
    def roi_image(self):
        """Return ROI image as TNImage."""
        if self.roi_image_path is not None:
            try:
                return TNImage(filename=self.roi_image_path)
            except Exception, e:
                raise e
        else:
            print "Connectome has no registered ROI image."
            return None

    @roi_image.setter
    def roi_image(self, filename):
        """Take path to roi image and set as our new roi image. This
        will overwrite previous roi_image if one exists.

        """
        self._roi_image_path = filename

    @property
    def subject_id(self):
        """Return subject id if it is set, otherwise return None."""
        try:
            return self.clinical_info['ID']
        except KeyError:
            return None
        except Exception, e:
            raise e

    @subject_id.setter
    def subject_id(self, new_id):
        """Set the subject id."""
        self.clinical_info['ID'] = new_id

    @property
    def type(self):
        """Return the type of the connectome as string. This could be
        one of the following {'PROBTRACKX'|'DTK'}. This is stored in a
        file w/in the connectome and is used to determine which
        subclass of connectome should handle reading.

        """
        try:
            return self.manifest['TYPE']
        except KeyError:
            print "Type could not be determined."
            return None

    @type.setter
    def type(self, new_type):
        """Set type to new value."""
        self.manifest['TYPE'] = new_type

    @property
    def version(self):
        """Return the version of the conectome. This will help in
        determining how to handle connectomes and how to deal with
        compatibility issues as the conenctome structure is developed
        further.

        """
        try:
            return self.manifest['VERSION']
        except KeyError:
            return None
        except Exception, e:
            raise e

    @version.setter
    def version(self, new_version):
        """Set version... this can only be done if a version does not
        already exist, so as to protect against version corruption.

        """
        if self.version is None:
            self.manifest['VERSION'] = new_version
        else:
            raise Exception("Cannot overwrite previous version.")

    @property
    def filename(self):
        """Internal attribute returning filename or None if filename
        has not been set.

        """
        try:
            return self._filename
        except:
            return None

    @property
    def clinical_info_path(self):
        """Return the path to our info dictionary if we have one, if
        not return None.

        """
        if self._filename is not None:
            return os.path.join(self._filename, __CLINICAL_INFO_FILENAME__)
        else:
            return None
            
    @property
    def network_path(self):
        """Return the path to our network if we have one, if not,
        return None.

        """
        try:
            return self._network_path
        except AttributeError:
            try:
                found_files = os.listdir(self._filename)
                for _file in found_files:
                    if _file == __NETWORK_FILENAME__:
                        return os.path.join(self._filename, _file)
                return None # we didn't find it
            except:
                return None

    @property
    def roi_image_path(self):
        """Return the path to our roi image if we have one, if not,
        return None.

        """
        try:
            return self._roi_image_path
        except:
            try:
                found_files = os.listdir(self._filename)
                for _file in found_files:
                    if _file.split('.')[0] == __ROI_IMAGE_PREFIX__:
                        return os.path.join(self._filename,_file)
                # didn't find in connectome dir
                return None
            except:
                return None
    
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
                    data = self.clinical_info
                except KeyError:
                    print "No clinical info found for connectome."
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

    def matrix_for_key(self, key, force_symmetric=True,
                       binarize=False, number_of_nodes=None,
                       zero_diagonal=True):
        """Return a NxN matrix for given connectome and key. The
        connection matrix is returned as a numpy ndarray.

        """
        ##TODO: Add functionality for fixed desnity
        # create our new shiny connection matrix
        import numpy
        if number_of_nodes is None:
            n = max(self.nodes())
        else:
            n = number_of_nodes
        new_cmat = numpy.zeros((n,n))

        # extract the value for key for every edge in the given connectome
        for i,j in self.edges_iter():
            new_cmat[i-1][j-1] = self[i][j][key]
                
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

        # if we need to binarize...
        if binarize:
            new_cmat = new_cmat >= 1

        # if we need to zero diagonal...
        if zero_diagonal:
            numpy.fill_diagonal(new_cmat,0)

        # return the cmat
        return new_cmat
    
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
                    data = copy(self.clinical_info)
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
            self.node[int(label)]['subject_position'] = tuple(np.mean
                                                              (np.where
                                                               (ROI_img_data == int(label)),
                                                               axis=1))

    def submatrix_for_key(self, submatrix_nodes, key):
        """Return a NxN matrix for key, where N ==
        len(submatrix_nodes). Submatrix nodes are first sorted, then
        metrics are extracted for each edge (i,j) in order according
        to sort.

        """
        import numpy
        submatrix_nodes.sort()        
        n = len(submatrix_nodes)
        new_cmat = numpy.zeros((n,n))
        for i in range(0,n):
            for j in range(0,n):
                node0 = submatrix_nodes[i]
                node1 = submatrix_nodes[j]
                if node0 == node1:
                    new_cmat[i][j] = 0
                else:
                    try:
                        new_cmat[i][j] = self[node0][node1][key]
                    except KeyError:
                        pass
                    except Exception, e:
                        raise e
        return new_cmat

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
        # CONNECTOME DIRECTORY ##############################################
        if not os.path.isdir(filename):
            os.mkdir(filename,0755)
        # GRAPH #############################################################
        networkx.write_gpickle(self.network,
                               os.path.join(filename,__NETWORK_FILENAME__))
        # VERSION and TYPE ##################################################        
        manifest_info = {'VERSION': self.version,
                         'TYPE': self.type,
                         'CREATED': self.created}
        try:
            manifest_file = open(os.path.join(filename,__MANIFEST_FILENAME__), 'wb')
            pickle.dump(manifest_info, manifest_file)
        except Exception, e:
            raise e
        finally:
            manifest_file.close()
        # CLINICAL INFO #####################################################
        try:
            clinfo_file = open(os.path.join(filename, __CLINICAL_INFO_FILENAME__), 'wb')
            pickle.dump(self.clinical_info, clinfo_file)
        except Exception, e:
            raise e
        finally:
            clinfo_file.close()
        # ROI IMAGE (if needed) #############################################
        if self.roi_image_path is not None:
            # get the file extension of the roi file, we will use this
            # later
            roi_file_ext = self.roi_image_path.split('.',1)[1]
            # generate standard path using extension
            roi_dest_path = os.path.join(filename, "%s.%s" %
                                         (__ROI_IMAGE_PREFIX__, roi_file_ext))
            # if roi image path is not equal to our standard path,
            # then we need to copy the roi image to our standard path
            if self.roi_image_path != roi_dest_path:
                shutil.copyfile(self.roi_image_path, roi_dest_path)
        
    def _log_modification(self, description):
        timestamp = datetime.datetime
        try:
            user = os.environ['USER']
        except KeyError:
            user = 'User could not be determined'
        self.modifications.append({'timestamp': timestamp,
                                   'user': user,
                                   'desc': description})

# MODULE LEVEL FUNCTIONS
def read(filename):
    """Read the connectome by sending it to the appropriate handler
    and return the connectome of approrpiate type.

    Input::

      filename - path to the connectome directory

    """
    # read the connectome manifest
    try:
        manifest_file = open(os.path.join(filename,__MANIFEST_FILENAME__), 'rb')
        manifest = pickle.load(manifest_file)
        _version = manifest['VERSION']
        _type = manifest['TYPE']
        manifest_file.close()
    except:
        raise Exception("Could not parse manifest belonging to %s." % filename)
    if _version < 2.0:
        raise Exception("Connectome must be version 2.0 or greater")
    if _type == 'DTK':
        from . import TNDtkConnectome
        return TNDtkConnectome(filename=filename)
    if _type == 'PROBTRACKX':
        from . import TNProbtrackxConnectome
        return TNProbtrackxConnectome(filename=filename)



####    
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
            label_i = int(roi_data[tuple(voxel_i)])
            label_j = int(roi_data[tuple(voxel_j)])
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

def generate_connectome_from_cmat(cmat,
                                  metric_key,
                                  roi_img = None,
                                  node_info = None):
    """Return a TNConnectome object

    Example
    -------
    import muscip.connectome as mcon
    C = mcon.connectome.generate_connectome_from_cmat(cmat, 'my_new_key')

    Input::

    [mandatory]

      cmat - connection matrix indicating some metric

      metric_key - key under which to store metric

    [optional]

      roi_img - a loaded nibabel image representing rois

      node_info - node information

    """
    C = TNConnectome()
    if roi_img is not None and node_info is not None:
        C.populate_node_info(roi_img.get_data(), node_info)
    for i in range(0,cmat.shape[0]):
        for j in range(0,cmat.shape[1]):
            C.add_edge(i,j)
            C[i][j][metric_key] = (cmat[i,j] + cmat[j,i]) / 2.0
    return C
