from . import connectome as con

class TNConnectomeGroup(object):
    """Group of connectomes with built-in functions for group
    analysis.

    """

    def __init__(self,
                 subject_dir=None,     # path to directory containing subjects
                 subjects=None,        # list of subjects
                 connectome_path=None  # path to append to
                                       # subject_dir/subject/... to
                                       # realise complete path to
                                       # connectome

    ):
        import os.path as op
        self.subject_dir = op.abspath(subject_dir)
        self.subjects = subjects
        self.connectome_path = connectome_path
        self._count = None
        self._info_keys = None
        self._metric_keys = None
        self._number_of_nodes = None

    @property
    def count(self):
        if self._count:
            return self._count
        else:
            ctx = 0
            for con in self.connectomes():
                ctx += 1
            self._count = ctx
            return self._count

    @property
    def info_keys(self):
        if self._info_keys:
            return self._info_keys
        else:
            self.load_info_keys()
            return self._info_keys
        
    @property
    def metric_keys(self):
        if self._metric_keys:
            return self._metric_keys
        else:
            self.load_metric_keys()
            return self._metric_keys

    @property
    def number_of_nodes(self):
        if self._number_of_nodes:
            return self._number_of_nodes
        else:
            self.guess_number_of_nodes()
            return self._number_of_nodes
            
    def connectomes(self):
        """Iterate and yield connectomes."""
        for subject in self.subjects:
            from os.path import join
            yield con.read(join(self.subject_dir, subject, self.connectome_path))

    def edge_list_for_key(self, key, outfile, clinical_info=False, sep=','):
        """Create an edge list for connectomes for a given key and
        output to a text file.

        Input::

          key - for which key do we want values?

          outfile - write to this filepath

          clinical_info - should each record contain all clinical
                          info?

          sep - charater to use as delimiter (defaults to comma)

        Example::

          aConnectomeGroup.edge_list_for_key('fiber_count', 'foo.csv')

        """
        # open outfile for writing
        hOutfile = open(outfile, 'wb')

        # if we are not writing clinical info, then the header should
        # be "subject.id,u,v,<key>"
        if not clinical_info:
            hdr = "subject.id%su%sv%s%s" % (sep,sep,sep,key)
        # otherwise, grab the first subject and get clinical
        # info... we assume that all subjects have the same clinical
        # info dictionary... if this is not true, then it may be best
        # to run this method with the clinical_info=False option
        else:
            hdr = ""
            C_exemplar = self.connectomes().next()
            # add an entry in header for every key in the clinical
            # info dictionary of the example connectome
            clinical_info_keys = C_exemplar.clinical_info.keys()
            for clinkey in clinical_info_keys:
                hdr += (clinkey + sep)
            # append our edge header stuff
            hdr += "u%sv%s%s" % (sep,sep,key)

        # write the header to the top of the outfile
        hOutfile.write(hdr + "\n")

        # for each connectome in group...
        for con in self.connectomes():
            # form a record prefix containing either only the subject
            # id, or the complete clinical info dictionary depending
            # on parameter
            if clinical_info:
                recordPrefix = ""
                for clinkey in clinical_info_keys:
                    recordPrefix += ( str(con.clinical_info[clinkey]) + sep)
            else:
                recordPrefix = (con.subject_id + sep)

            # get edge list for key
            edgeList = con.edge_list_for_key(key)
            # for each record in edgeList...
            for record in edgeList:
                # append record prefix and write this record
                thisRecord = recordPrefix + str(record['u']) + sep + \
                    str(record['v']) + sep + str(record[key])
                hOutfile.write(thisRecord + '\n')

        # close up outfile
        hOutfile.close()

    def export_to_matlab(self, filename,
                         fiber_lengths=False,
                         fiber_lengths_key='fiber_lengths',
                         number_of_nodes=None,
                         subnetwork_nodes=None):
        structure = dict()
        structure['data'] = list()
        exclude_keys = ['streamlines', 'streamlines_length', 'fiber_lengths', 'fibers']
        if number_of_nodes is None:
            number_of_nodes = self.number_of_nodes
        for connectome in self.connectomes():
            # grab info for connectome
            record = dict()
            for key in self.info_keys:
                try:
                    record[key] = connectome.clinical_info[key]
                except KeyError:
                    record[key] = None
            # grab metrics for connectome
            for key in self.metric_keys:
                if key in exclude_keys:
                    continue
                if subnetwork_nodes is None:
                    record[key] = connectome.matrix_for_key(key,
                                                            number_of_nodes=number_of_nodes)
                else:
                    record[key] = connectome.submatrix_for_key(subnetwork_nodes, key)
                # add fiber length histograms if requested
                if fiber_lengths:
                    if not fiber_lengths_key in self.metric_keys:
                        raise Exception("%s is not the correct fiber length key" % fiber_lengths_key)
                    record[fiber_lengths_key] = connectome.fiber_lengths
            structure['data'].append(record)
            from scipy.io import savemat
            savemat(filename, structure)
        
    def dataFrame(self, info_keys=None, metric_keys=None, number_of_nodes=None):
        try:
            df = dict()
            import numpy as np

            # if info keys is none...assume we want them all
            if not info_keys:
                info_keys = list()
                for C in self.connectomes():
                    try:
                        for key in C.clinical_info.keys():
                            if key not in info_keys:
                                info_keys.append(key)
                    except Exception as e:
                        raise e
                    
            for C in self.connectomes():
                for key in info_keys:
                    try:
                        if key not in df.keys():
                            df[key] = list()
                        df[key].append(C.clinical_info[key])
                    except KeyError:
                        continue
                    except Exception as e:
                        raise e
                for key in metric_keys:
                    try:
                        if key not in df.keys():
                            df[key] = list()
                        df[key].append(C.matrix_for_key(key))
                    except KeyError:
                        continue
                    except Exception as e:
                        raise e
            import pandas
            return pandas.DataFrame(df)
        except Exception as e:
            raise e
        
    def visualize_adj_matrix(self,
                             function=None,
                             group_key=None,
                             metric_key=None,
                             number_of_nodes=None,
                             zero_diagonal=True,
                             binarize=False,
                             fixed_density=None
    ):
        pass

    def load_info_keys(self):
        self._info_keys = list()
        for C in self.connectomes():
            try:
                for key in C.clinical_info:
                    if key not in self._info_keys:
                        self._info_keys.append(key)
            except KeyError:
                continue
            except Exception, e:
                raise e
        self._info_keys.sort()

    def load_metric_keys(self):
        self._metric_keys = list()
        try:
            for C in self.connectomes():
                for u,v,data in C.network.edges_iter(data=True):
                    for key in data.keys():
                        if key not in self._metric_keys:
                            self._metric_keys.append(key)
        except Exception, e:
            raise e
        self._metric_keys.sort()

    def guess_number_of_nodes(self):
        guess = 0
        for connectome in self.connectomes():
            max_node = max(connectome.network.nodes())
            if max_node > guess:
                guess = max_node
        return guess

    def qa_image(self, filename=None, show=True, size=(10,10), title=None):
        """Display or save a figure of all the log fiber count
        matrices in the connectome group. This can be used to quickly
        spot outlies where problems may exist.

        Input::

          filename - filename to save output image, image will not be
          saved to disk if filename is not provided

          show - show the figure in a window

          size - tuple in inches representing the size of the image

          title - title to place on the figure as a whole

        """
        import matplotlib.pyplot as plot, math, numpy as np, warnings
        FIBER_COUNT_KEY = "fiber_count"
        ncols = int(math.ceil(math.sqrt(self.count)))
        fig = plot.figure(figsize=size)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for i, c in enumerate(self.connectomes()):
                # get the adjacency matrix of the natural log of fiber counts
                mat = c.matrix_for_key(FIBER_COUNT_KEY)
                mat = np.log(mat)
                mat[np.where(mat < 0)] = 0
                ax = fig.add_subplot(ncols,ncols,i+1)
                ax.imshow(mat, interpolation="nearest")
                try:
                    ax.set_title(c.subject_id, size='x-small')
                except Exception, e:
                    pass
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
        if title:
            fig.suptitle(title, size='x-large')
        if filename:
            fig.savefig(filename, orientation='landscape')
        if show:
            plot.show()
        
    def set_info(self, csvfile, id_key, new_filename=None):
        
            import os.path as op
            # create the csv reader and grab data
            try:
                import csv
                read_info = dict()
                f = open(op.abspath(csvfile), 'rt')
                reader = csv.DictReader(f)
                for record in reader:
                    read_info[record[id_key]] = record
            except Exception, e:
                print e
            finally:
                f.close()
            # for every connectome entry, try and get info from the info we
            # previously read
            for idx, C in enumerate(self.connectomes()):
                # try to get info for subject
                try:
                    subject = self.subjects[idx]
                    new_info = read_info[subject]
                except KeyError:
                    print "No info found in csv file for: %s" % subject
                    continue
                # apply new info
                for key, value in new_info.iteritems():
                    C.clinical_info[key] = value
                # try to save connectome
                try:
                    if new_filename:
                        outpath = op.join(self.subject_dir, subject, new_filename)
                    else:
                        outpath = op.join(self.subject_dir, subject, self.connectome_path)
                    C.write(outpath)
                except Exception, e:
                    print "Could not save connectome file -- %s" % e
