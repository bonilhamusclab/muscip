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

    def visualize_adj_matrix(self,
                             function=None,
                             group_key=None,
                             metric_key=None,
                             number_of_nodes=None,
                             zero_diagonal=True,
                             binarize=False,
                             fixed_density=None,
    ):
        try:
            import numpy as np
            # create an empty numpy array to hold all matrices
            stack = np.zeros((len(self.subjects), number_of_nodes,
                              number_of_nodes))
            # get a count of groups
            groups = dict()
            for idx, C in enumerate(self.connectomes()):
                group = C.graph['info'][group_key]
                if group not in groups.keys():
                    groups[group] = [idx]
                else:
                    groups[group].append(idx)
                stack[idx,:,:] = C.matrix_for_key(metric_key,
                                                  binarize=binarize,
                                                  number_of_nodes=number_of_nodes,
                                                  zero_diagonal=zero_diagonal)
            # generate the result
            submatrices = list()
            for key in groups.keys():
                print key, groups[key]
                print stack[groups[key],:,:].shape
                submatrices.append(stack[groups[key],:,:].reshape((len(stack[groups[key]]),
                                                                   number_of_nodes ** 2)))
            result = apply(function, submatrices) #.reshape((number_of_nodes, number_of_nodes))
            print result
            # display the result
            from matplotlib.pyplot import imshow
            imshow(result, interpolation='nearest')

        except Exception as e:
            raise e
