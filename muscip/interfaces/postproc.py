import os
from nipype.interfaces.base import (BaseInterface, File, TraitedSpec,
                                    isdefined)
from nipype.utils.filemanip import split_filename
# tn_image_processing will soon change to muscip
from .. import connectome
from .. import fibers
import nibabel
import networkx

from nipype import logging
iflogger = logging.getLogger('interface')

class ConnectomeGeneratorInputSpec(TraitedSpec):
    network_file = File(genfile=True, desc='NetworkX graph that will be generated')
    roi_file = File(exists=True, mandatory=True, desc='ROI volume')
    track_file = File(exists=True, mandatory=True, desc='Trackvis track file')
    wm_file = File(exists=True, mandatory=True, desc='White matter mask used in tractography')

class ConnectomeGeneratorOutputSpec(TraitedSpec):
    network_file = File(desc="NetworkX graph representing the connectome", exists=True)

class ConnectomeGenerator(BaseInterface):
    """Maps connectome and outputs the result as a NetworkX graph

    Example
    -------
    >>> import nipype.interfaces.muscip as muscip
    >>> ng = muscip.NetworkGenerator()
    >>> ng.inputs.network_filename = 'desired_output.pkl'
    >>> ng.inputs.roi_file = 'path/to/roi/volume'
    >>> ng.inputs.track_file = 'path/to/trackfile.trk'
    >>> ng.inputs.wm_file = 'path/to/wm_file'
    >>> ng.run() # doctest: +SKIP
    """

    input_spec = ConnectomeGeneratorInputSpec
    output_spec = ConnectomeGeneratorOutputSpec

    def _run_interface(self, runtime):
        if isdefined(self.inputs.network_file):
            path, name, _ = split_filename(self.inputs.network_file)
            self.inputs.network_file = os.path.abspath(name + '.pkl')
        else:
            self.inputs.network_filename = 'connectome.pkl'
        iflogger.info('Loading ROI and WM images...')
        roi = nibabel.load(self.inputs.roi_file)
        wm = nibabel.load(self.inputs.wm_file)
        iflogger.info('Loading .trk file...')
        fib = fibers.read(self.inputs.track_file)
        iflogger.info('Generating connectome...')
        con = connectome.generate_connectome(fib, roi)
        iflogger.info('Extracting Hagmann density...')
        con.populate_hagmann_density(wm, roi)
        iflogger.info('Saving output...')
        networkx.write_gpickle(con, os.path.abspath(self.inputs.network_file))
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['network_file'] = self.inputs.network_file
        return outputs
