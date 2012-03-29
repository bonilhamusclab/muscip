"""Dummy workflow to test behavior of nipype."""

import traits
from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine
)

class EchoInputSpec(CommandLineInputSpec):
    in_text = traits.trait_types.Str(desc = "In String", mandatory = True, position = 0, argstr="-e %s")

class EchoOutputSpec(TraitedSpec):
    out_text = traits.trait_types.Str(desc = "In String", mandatory = True, position = 0, argstr="-e %s")

class Echo(CommandLine):
    _cmd = 'echo'
    input_spec = EchoInputSpec
    output_spec = EchoOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_text'] = "This is what I said: %s" % self.inputs.in_text
        return outputs
