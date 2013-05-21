from nibabel import trackvis
from .fibers import TNFibers, fiber_length

__FORMAT__ = 'trackvis'

class TNTrackvisFibers(TNFibers):

    def __init__(self, **kwargs):
        TNFibers.__init__(self, **kwargs)
        self._format = __FORMAT__
        self._hdr = kwargs.get('hdr', None)

    @property
    def fibers(self):
        def fibers_generator(fibers):
            for rec in fibers:
                yield rec[0]
        data, hdr = trackvis.read(self._filename, as_generator=True,
                                  points_space='voxel')
        return fibers_generator(data)

    @property
    def hdr(self):
        if self._hdr is not None:
            return self._hdr
        else:
            return trackvis.empty_header()
            
    @property
    def number_of_fibers(self):
        try:
            return self.hdr['n_count']
        except:
            return None

    @property
    def shape(self):
        try:
            return self.hdr['dim']
        except:
            return None
    @property
    def spacing(self):
        return 'voxel' # for trackvis, we will always handle in voxel
                       # spacep

    @property
    def voxel_size(self):
        try:
            return self.hdr['voxel_size']
        except:
            return None
        
    def write(self, filename, **kwargs):
        endianness = kwargs.get('endianness', None)
        hdr_mapping = kwargs.get('hdr_mapping', self.hdr)
        points_space = kwargs.get('points_space', 'voxel')
        streamlines = kwargs.get('streamlines', self.fibers)
        def streamlines_gen(streamlines):
            for streamline in streamlines:
                yield (streamline, None, None)
        trackvis.write(filename, streamlines_gen(streamlines),
                       hdr_mapping=hdr_mapping, endianness=endianness,
                       points_space=points_space)

    def _trackvis_header_fields(self):
        return {
            'id_string',
            'dim',
            'voxel_size',
            'origin',
            'n_scalars',
            'scalar_name',
            'n_properties',
            'property_name',
            'vox_to_ras',
            'reserved',
            'voxel_order',
            'pad2',
            'image_orientation_patient',
            'pad1',
            'invert_x',
            'invert_y',
            'invert_z',
            'swap_xy',
            'sway_yz',
            'swap_zx',
            'n_count',
            'version',
            'hdr_size'
        }
        
################################################################################

def centroid_for_tracks(tracks, length_factor=None, point_factor=2.0):
    """Use the Quick Bundle algorithm to generate a centroid for a
    given group of tracks."""
    import dipy.segment.quickbundles as qb
    import numpy
    # if average length is none, then we need to calculate average
    # length from the fibers themselves so that we may multiply by our
    # point factor
    if length_factor is None:
        lengths = []
        for track in tracks:
            lengths.append(fiber_length(track))
        length_factor = numpy.mean(lengths)
    # use quick bundle with an outrageously large distance threshold
    # so that we are only returning one track cluster... this is not
    # how this tool is meant to be used, but for our simple purpose, i
    # believe this is what we need... also we will make the number of
    # points be a function of the track length
    Bundle = qb.QuickBundles(tracks, dist_thr=1000000,
                             pts=int(length_factor * point_factor))
    return Bundle.virtuals()[0]

def read_trackvis(filename):
    data, hdr = trackvis.read(filename, as_generator=True)
    return TNTrackvisFibers(filename=filename, hdr=hdr)

def write_trackvis(filename, **kwargs):
    endianness = kwargs.get('endianness', None)
    hdr_mapping = kwargs.get('hdr_mapping', None)
    points_space = kwargs.get('points_space', 'voxel')
    streamlines = kwargs.get('streamlines', [])
    def streamlines_gen(streamlines):
        for streamline in streamlines:
            yield (streamline, None, None)
    trackvis.write(filename, streamlines_gen(streamlines),
                   hdr_mapping=hdr_mapping, endianness=endianness,
                   points_space=points_space)
    
