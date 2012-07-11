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
        self._info_keys = None
        self._metric_keys = None

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
            
    def connectomes(self):
        """Iterate and yield connectomes."""
        try:
            for subject in self.subjects:
                from os.path import join
                yield con.read_gpickle(join(self.subject_dir,
                                            subject,
                                            self.connectome_path))
        except Exception as e:
            raise e

    def dataFrame(self, info_keys=None, metric_keys=None, number_of_nodes=None):
        try:
            df = dict()
            # if info keys is none...assume we want them all
            if not info_keys:
                info_keys=self.info_keys
            # if metric_keys is none...assume we want them all
            if not metric_keys:
                metric_keys=self.metric_keys
            # loop through connectomes and get values
            for C in self.connectomes():
                for key in info_keys:
                    if key not in df.keys():
                        df[key] = list()
                    df[key].append(C.graph['info'][key])
                for key in metric_keys:
                    excluded_keys = ('streamlines', 'streamlines_length')
                    if key not in df.keys():
                        df[key] = list()
                    df[key].append(C.matrix_for_key(key))
            # convert to pandas data frame
            import pandas
            return pandas.DataFrame(df)
        except Exception as e:
            print e

    def export_to_matlab_struct(self):
        struct = dict()
        struct['Sub'] = list()
        ##TODO: finish implementation
        pass

        
    def load_info_keys(self):
        self._info_keys = list()
        for C in self.connectomes():
            try:
                for key in C.graph['info']:
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
                for u,v,data in C.edges_iter(data=True):
                    for key in data.keys():
                        if key not in self._metric_keys:
                            self._metric_keys.append(key)
        except Exception, e:
            raise e
        self._metric_keys.sort()
        
    def visualize_adj_matrix(self,
                             function='mean',
                             group_key=None,
                             metric_key=None,
                             number_of_nodes=None,
                             zero_diagonal=True,
                             binarize=False,
                             fixed_density=None
    ):
        # check that function is valid
        valid_functions = ['tstat', 'mean', 'max', 'min']
        if function not in valid_functions:
            raise Exception('Function must be one of %s' %
                            valid_functions)
        # get data frame
        df = self.dataFrame(info_keys=(group_key,),
                            metric_keys=(metric_key,),
                            number_of_nodes=number_of_nodes)
        # function cases:
        results = dict()
        import numpy as np
        if function=='mean':
            results['title'] = "Mean %s" % metric_key
            results['subplots'] = list()
            matrices = df.groupby(group_key)[metric_key].apply(np.mean)
            for key in matrices.keys():
                results['subplots'].append({'title': key, 'data':
                                            matrices[key]})
        if function=='max':
            results['title'] = "Max %s" % metric_key
            results['subplots'] = list()
            matrices = df.groupby(group_key)[metric_key].apply(np.nanmax)
            for key in matrices.keys():
                results['subplots'].append({'title': key, 'data':
                                            matrices[key]})
        if function=='min':
            results['title'] = "Min %s" % metric_key
            results['subplots'] = list()
            matrices = df.groupby(group_key)[metric_key].apply(np.nanmin)
            for key in matrices.keys():
                results['subplots'].append({'title': key, 'data':
                                            matrices[key]})
        if function=='ttest':
            ##TODO: complete implementation
            pass
            # def ttest(series, first_group, second_group):
            #     a = series[first_group]
            #     b = series[second_group]
            #     t, _ = ttest_ind(a,b,axis=2)
            #     return t
                
            # results['title'] = "T-Values %s" % metric_key
            # results['subplots'] = list()
            # from itertools import combinations
            # comparisons = list(combinations(df.groupby(group_key)[metric_key].groups.keys(), 2))
            # for comparison in comparisons:
            #     results['subplot'].append({'title': "%s > %s" %
            #                                (comparison[0], comparison[1]),
            #                                'data': 

        # plot results
        import matplotlib.pyplot as plt
        plt.figure(1)
        plt.title(results['title'])
        # determine shape
        from math import ceil, floor, sqrt
        num_subplots = len(results['subplots'])
        ncols = ceil(sqrt(num_subplots))
        nrows = floor(sqrt(num_subplots))
        if nrows * ncols < num_subplots:
            nrows += 1
        # get min, max
        global_min = np.inf
        global_max = -np.inf
        for subplot in results['subplots']:
            if np.max(subplot['data']) > global_max:
                global_max = np.max(subplot['data'])
            if np.min(subplot['data']) < global_min:
                global_min = np.min(subplot['data'])
        for idx, subplot in enumerate(results['subplots']):
            plt.subplot(nrows,ncols,idx)
            plt.imshow(subplot['data'], interpolation='nearest',
                       vmax=global_max, vmin=global_min)
            plt.title(subplot['title'])
        plt.show()

    def set_info(self, csvfile, id_key, new_filename=None):
        
            import os.path as op
            from copy import copy
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
                C.set_info(new_info)
                # try to save connectome
                try:
                    if new_filename:
                        outpath = op.join(self.subject_dir, subject,
                                          op.dirname(self.connectome_path),
                                          new_filename)
                    else:
                        outpath = op.join(self.subject_dir, subject, self.connectome_path)
                    C.write(outpath)
                except Exception, e:
                    print "Could not save connectome file -- %s" % e
