import nibabel

class Freesurfer(object):
    """Given a directory that has been run through recon-all, this
    object loads and wraps freesurfer atlas and white-matter masks.

    Example
    -------
    >>> import muscip.atlas.freesurfer as fsa
    >>> fs_atlas = fsa.load('processed-fs-dir')
    >>> wm_mask = fs_atlas.wm_mask
    >>> rois = fs_atlas.roi_img

    Attributes:
      gm_lut - grey matter lookup table
      roi_img - the ROI atlas as a nibabel image
      wm_lut - white matter lookup table
      wm_mask - the white-matter mask as a nibabel image

    """

    def __init__(self, fs_dir=None):
        """Return a freesurfer atlas with labels as defined by
        look-up-tables

        """
        self.fs_dir = fs_dir
        self.gm_lut = _def_gm_lut()
        self.node_info = _def_node_info()
        self.roi_img = None
        self.wm_lut = _def_wm_lut()
        self.wm_mask = None

    def extract_stats(self):
        """Parse freesurfer stat files and return a dictionary with
        entries for each node as well as include entries for file-wide
        measures.

        """
        # create an empty dictionary, with entries for each node
        # w/ key and name, as well as a reverse lookup by name
        import networkx as nx
        struct = nx.Graph()
        label_by_name = dict()
        for label, data in self.node_info.nodes_iter(data=True):
            struct.add_node(label,{'name': data['name']})
            label_by_name[data['name']] = label
        # get lh-cortex
        self._parse_stats('lh.aparc.stats', 'ctx-lh-', struct, label_by_name)
        # get rh-cortex
        self._parse_stats('rh.aparc.stats', 'ctx-rh-', struct, label_by_name)
        # get sub-cortical
        self._parse_stats('aseg.stats', None, struct, label_by_name)
        # return result
        return struct

    def _parse_stats(self, stats_file, prefix, out_struct, label_lut):
        from os.path import join
        fhandle = open(join(self.fs_dir, 'stats', stats_file), 'rb')
        header = None
        header_label = 'ColHeaders'
        if prefix is None:
            prefix = '' # we are appending to this, so it must be an
                        # empty string...
        for line in fhandle.readlines():
            # we will reach header first, other processing depends
            # on succesfully parsing header for keys, so let us
            # first get that..
            if header is None:
                # if this is the header...
                try:
                    if line.split()[1] == header_label:
                        # ...get the keys from the header
                        header = line.split()[2:]
                        label_col_idx = header.index('StructName')
                # exception will occur if this is just an empty
                # comment line
                except:
                    continue
            # ...if header is not None, then we can assume it has
            # already been processed
            else:
                # check to make sure line is not a comment, if
                # it is, move on to next line
                if line.split()[0] == '#':
                    continue
                # otherwise, first term is prefix-less-name of region
                # to which we must find label - we must also handle
                # the case where we find a region in the stats file
                # that is not in our expected regions, as our expected
                # regions may be a subset of all regions
                try:
                    region_label = label_lut[prefix + line.split()[label_col_idx]]
                except KeyError:
                    continue
                except Exception, e:
                    raise e
                # ...and we need to get values for each stat label
                for idx, stat in enumerate(line.split()):
                    out_struct.node[region_label][header[idx]] = stat

                    
# MODULE LEVEL FUNCTIONS
###################################################
def load(freesurfer_dir, gm_lut=None, wm_lut=None):
    """Return an instance of our Freesurfer atlas class"""
    newAtlas = Freesurfer(fs_dir=freesurfer_dir)    
    # load new atlas starting from aparc+aseg.mgz
    #############################################
    # convert mgz files to nii.gz files
    from os.path import join
    aparc_aseg = _convert_mgz_to_niigz_and_load(join(freesurfer_dir,
                                                     'mri/aparc+aseg.mgz'))
    aparc_aseg_data = aparc_aseg.get_data()
    # load any custom luts
    ######################
    if gm_lut:
        newAtlas.gm_lut = gm_lut
    if wm_lut:
        newAtlas.wm_lut = wm_lut
    # map fs labels to expected labels
    ##################################
    import numpy as np
    roi_data = np.zeros(aparc_aseg_data.shape, dtype=np.uint8)
    for map_entry in newAtlas.gm_lut:
        roi_data[ aparc_aseg_data == map_entry[1]] = map_entry[0]
    wm_data = np.zeros(aparc_aseg_data.shape, dtype=np.uint8)
    for map_entry in newAtlas.wm_lut:
        wm_data[ aparc_aseg_data == map_entry] = 1
    # re-orient wm and roi images to orig
    #####################################
    newAtlas.roi_img = _reorient_fs_to_orig(
        nibabel.nifti1.Nifti1Image(roi_data, aparc_aseg.get_affine() ,
                                   aparc_aseg.get_header()), freesurfer_dir)
    newAtlas.wm_mask = _reorient_fs_to_orig(
        nibabel.nifti1.Nifti1Image(wm_data, aparc_aseg.get_affine() ,
                                   aparc_aseg.get_header()), freesurfer_dir)
    return newAtlas
    
def _convert_mgz_to_niigz_and_load(filepath):
    """Return the given mgz file as a nii loaded and wrapped in
    nibabel

    """
    import os.path as op
    import subprocess
    import tempfile
    out_filename = op.splitext(op.split(filepath)[1])[0] + '.nii.gz'
    tmpdir = tempfile.mkdtemp()
    tmpfile = op.join(tmpdir, out_filename)
    print "IN: %s\tOUT: %s" % (filepath, tmpfile)
    try:
        subprocess.check_call(["mri_convert", filepath, tmpfile])
        tmpimg = nibabel.load(tmpfile)
        return nibabel.nifti1.Nifti1Image(tmpimg.get_data(),
                                          tmpimg.get_affine(),
                                          header=tmpimg.get_header())
    except Exception as e:
        print e
        return None
    finally:
        from shutil import rmtree
        try:
            rmtree(tmpdir)
        except:
            pass

def _def_gm_lut():
    """Returnt the path to default grey matter lookup table"""
    try:
        import pkgutil
        data = pkgutil.get_data(__name__, 'data/freesurfer_gm_lut.txt')
    except ImportError:
        import pkg_resources
        data = pkg_resources.resource_string(__name__, 'data/freesurfer_gm_lut.txt')
    read_lut = []
    try:
        for line in data.split('\n'):
            if not line == "":
                in_line = line.split()
                read_lut.append([eval(in_line[0]), eval(in_line[1])])
    except Exception, e:
        print e
        return None
    return read_lut

def _def_node_info():
    """Return node info as networkx object with data for each node"""
    import networkx
    import os.path
    module_path = os.path.split(__file__)[0]
    node_info_path = os.path.join(module_path,'data/freesurfer_node_info.graphml')
    return networkx.read_graphml(node_info_path)

def _def_wm_lut():
    """Returnt the path to default white matter lookup table"""
    try:
        import pkgutil
        data = pkgutil.get_data(__name__, 'data/freesurfer_wm_lut.txt')
    except ImportError:
        import pkg_resources
        data = pkg_resources.resource_string(__name__, 'data/freesurfer_wm_lut.txt')
    read_lut = []
    try:
        for line in data.split('\n'):
            if not line == "":
                read_lut.append(eval(line))
    except Exception, e:
        print e
        return None
    return read_lut
    
def _reorient_fs_to_orig(input_img, freesurfer_dir):
    """Freesurfer outputs images in its own space, this utility
    re-samples to dims found in first image in orig folder.

    Input:
      input_img - nibabel image to convert
      freesurfer_dir - root freesurfer directory

    Output:
      nibabel image in original space

    """
    # find reference image in orig directory
    import os.path as op
    from os import listdir
    import subprocess
    import tempfile    
    search_path = op.join(freesurfer_dir, 'mri/orig')
    found_files = listdir(search_path)
    found_files.sort()
    reference = op.join(search_path, found_files[0])
    tmpdir = tempfile.mkdtemp()
    tmp_in = op.join(tmpdir, 'infile.nii.gz')
    tmp_out = op.join(tmpdir, 'outfile.nii.gz')
    try:
        nibabel.save(input_img, tmp_in)
        subprocess.check_call(["mri_convert", "-rl", reference, "-rt",
                               "nearest", tmp_in, "-nc", tmp_out])
        tmp_img = nibabel.load(tmp_out)
        return nibabel.nifti1.Nifti1Image(tmp_img.get_data(),
                                          tmp_img.get_affine(),
                                          header=tmp_img.get_header())
    except Exception as e:
        print e
        return None
    finally:
        from shutil import rmtree
        try:
            rmtree(tmpdir)
        except:
            pass
    
    