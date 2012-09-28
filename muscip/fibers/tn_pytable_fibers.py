import tables
from .fibers import TNFibers

__FORMAT__ = 'PyTables'

class TNPyTableFibers(TNFibers):
    """Fiber implementation built on top of PyTables to appropriately
    deal with large datasets.

    """
    def __init__(self, **kwargs):
        TNFibers.__init__(self, **kwargs)
        self._format = __FORMAT__
        # initialize tables

    @property
    def h5f(self):
        try:
            return self._h5f
        except AttributeError:
            if self.filename is not None:
                self._h5f = self.filename
            else:
                import tempfile
                _, self._h5f = tempfile.mkstemp()
                self._h5f.close()
            return self._h5f
        
