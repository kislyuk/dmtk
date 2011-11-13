import logging
import numpy

from dmtk.io.cmph5 import Basemap

class ConsensusReferenceMapper:
    """
    Provide a mapping between reference coordinate space and consensus coordinate space for a given cmp.h5 and variants dataset.
    
    Usage:
    from dmtk.io import cmph5, GffIO
    cmp_fh = cmph5.factory.create(cmp_fn)
    variants_fh = GffIO.GffReader(gff_fn)
    mapper = ConsensusReferenceMapper(cmp_fh, variants_fh)
    for ref_info in cmp_fh.refInfoIterator():
        ref2cons = mapper.referenceToConsensus(ref_info)
        cons2ref = mapper.consensusToReference(ref_info)
        ref_pos = ref2cons[consensus_pos]
        consensus_pos = cons2ref[ref_pos]
    ref_names can be used to restrict the mapper to only given references:
    mapper = ConsensusReferenceMapper(cmp_fh, variants_fh, ref_names=['chrI', 'chrII'])
    gap_value governs how gaps are represented:
     - None: gaps in the target mapping are given with the coordinate of the preceding non-gap column
     - Defined: The value opposite a gap is set to gap_value. Note gap_value must be an integer (e.g. -1)
    """
    def __init__(self, cmp_filehandle, variants_gff_reader, ref_names=None, gap_value=None):
        self._consensus_dsets = {}
        self._ref2cons = {}
        self._cons2ref = {}
        self._ref_fullname2name = {}
        self._consensus_in_ref_space = {}
        self.gap_value = gap_value
        
        self._indels = {}
        for feature in variants_gff_reader:
            if feature.type == "insertion" or feature.type == "deletion":
                self._indels.setdefault(feature.seqid, [])
                self._indels[feature.seqid].append(feature)
        
        for ref_info in cmp_filehandle.refInfoIterator():
            if ref_names and ref_info.fullName not in ref_names:
                continue
            ref_group = cmp_filehandle.refGroupFromFullName(ref_info.fullName)
            self._ref_fullname2name[ref_info.fullName] = ref_group.name.lstrip("/") # check me
            if ref_group.hasConsensus:
                consensus_array = numpy.empty_like(ref_group.consensusCalls._dataset)
                ref_group.consensusCalls._dataset.read_direct(consensus_array)

#                import h5py
#                f=h5py.File("/home/UNIXHOME/stang/bugzilla/andreyBug/cons0.h5", 'r')
#                consensus_array = numpy.empty_like(f["/ref000001/Consensus/ConsensusCalls"])
#                f["/ref000001/Consensus/ConsensusCalls"].read_direct(consensus_array)
#                consensus_array = numpy.squeeze(numpy.transpose(consensus_array))
                
                self._consensus_dsets[ref_info.fullName] = "".join(Basemap[consensus_array])
                self._computeConsensusReferenceMapping(ref_info)
    
    def _computeConsensusReferenceMapping(self, ref_info):
        ref_name = ref_info.fullName
        consensus_dset = self._consensus_dsets[ref_name]
        ref2cons_diff = numpy.empty(ref_info.length, dtype=int)
        ref2cons_diff.fill(1)
        ref2cons_diff[0] = 0
        cons2ref_diff = numpy.empty(len(consensus_dset), dtype=int)
        cons2ref_diff.fill(1)
        cons2ref_diff[0] = 0
        
        my_ref_id = self._ref_fullname2name[ref_name]
        my_indels = self._indels[my_ref_id] if my_ref_id in self._indels else []
        my_indels = sorted(my_indels, key=lambda(f): f.start)
        cur_consensus_shift = 0
        
        for indel in my_indels:
            indel_length = int(indel.attributes["length"])
            if indel.type == "insertion":
                ref2cons_diff[indel.start-1] += indel_length
                for p in range(cur_consensus_shift+indel.start-1, cur_consensus_shift+indel.start+indel_length-1):
                    print indel, "\n\t", p, cur_consensus_shift, "Ref/cons array len:", len(ref2cons_diff), len(cons2ref_diff)
                    try:
                        cons2ref_diff[p] = 0
                    except IndexError:
                        logging.error("BUGCHECK: consensus generation error")
                        pass
                        print 'Failed:', indel, "\n\t", p, cur_consensus_shift, "Ref/cons array len:", len(ref2cons_diff), len(cons2ref_diff)
                cur_consensus_shift += indel_length
            else: # deletion
                cons2ref_diff[indel.start-1 + cur_consensus_shift] += indel_length
                for p in range(indel.start-1, indel.start+indel_length-1):
                    print indel, "\n\t", p, cur_consensus_shift, "Ref/cons array len:", len(ref2cons_diff), len(cons2ref_diff)
                    ref2cons_diff[p] = 0
                cur_consensus_shift -= indel_length
        
        self._ref2cons[ref_name] = numpy.cumsum(ref2cons_diff)
        self._cons2ref[ref_name] = numpy.cumsum(cons2ref_diff)
        
        if self.gap_value:
            self._ref2cons[ref_name][ref2cons_diff == 0] = self.gap_value
            self._cons2ref[ref_name][cons2ref_diff == 0] = self.gap_value

    def referenceToConsensus(self, ref_info):
        return self._ref2cons[ref_info.fullName]
        
    def consensusToReference(self, ref_info):
        return self._cons2ref[ref_info.fullName]

    def alignedReference(self, ref_info):
        return

    def alignedConsensus(self, ref_info):
        return
    
    def _computeConsensusInRefSpace(self, ref_info):
        cons_dset = self._consensus_dsets[ref_info.fullName]
        ref2cons = self.referenceToConsensus(ref_info)
        cons_in_ref_space = [cons_dset[ref2cons[i]] if ref2cons[i] != -1 else '-' for i in range(ref_info.length)]
        self._consensus_in_ref_space[ref_info.fullName] = "".join(cons_in_ref_space)

    def consensusInRefSpace(self, ref_info):
        ''' Return consensus as a string in reference space,
        i.e. same length as reference, with no insertions and gaps at deletions.
        '''
        if ref_info.fullName not in self._consensus_in_ref_space:
            self._computeConsensusInRefSpace(ref_info)
        return self._consensus_in_ref_space[ref_info.fullName]
