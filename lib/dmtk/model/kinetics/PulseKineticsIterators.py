"""
PulseKineticsUtils: Functions for loading and tabulating SMRT sequencing kinetic data,
and computing statistics for inferring DNA modifications

Author: Andrey Kislyuk
"""

import os, sys, itertools, logging, re, bz2, csv
import numpy, scipy, scipy.stats, random
from numpy import mean, median, std, var, sum, sqrt, sort, log, log10, nansum, nonzero, product, floor
from scipy.stats import ks_2samp, ttest_ind, ttest_rel
from numpy.random import randint
from multiprocessing import Pool, cpu_count

from dmtk.io import cmph5
from dmtk.io.cmph5 import CmpH5ColIterators


#class KineticContext:
#    """ Data structure for kinetic information at a read position. See instantiations below for field names. """
#    def __init__(self, **entries): self.__dict__.update(entries)
#    
#    def __str__(self): return self.__class__.__name__+': '+self.__dict__.__str__()


class KineticContextSet:
    """ Data structure for a set of observations of pulse kinetics at and around an aligned reference position. """
    def __init__(self, ref_pos, strand, aln_ctxt, extended_aln_ctxt, metrics,
                 alignment_ids=None, zmw_ids=None, subread_counts=None):
        self.ref_pos, self.strand, self.aln_ctxt, self.extended_aln_ctxt, self.metrics = ref_pos, strand, aln_ctxt, extended_aln_ctxt, metrics
        self.alignment_ids = alignment_ids
        self.zmw_ids = zmw_ids
        self.subread_counts = subread_counts

    def __str__(self): return self.__class__.__name__+': '+self.__dict__.__str__()


class AlignedPulseKineticsIterator:
    ''' Iterates over reference positions. For each position and strand, emits
    arrays of pulse metrics measured in reads aligned to those positions.
    
    The default iterator returns a tuple (reference position, metrics).
    Metrics is a dict keyed by metric, then strand, containing lists of metric values.
    '''
    def __init__(self, cmpH5, ref_info, ref_start=0, ref_end=None,
                 pulse_tables=['PulseWidth', 'IPD'], template_length=None,
                 pulse_metric_normalizer=None,
                 aln_context_filter=None,
                 target_coverage_ceiling = None,
                 upstream_enzyme_footprint = 2, downstream_enzyme_footprint = 5):
        
        self.pulse_tables, self.ref_info, self.template_length = pulse_tables, ref_info, template_length
        self.upstream_enzyme_footprint, self.downstream_enzyme_footprint = upstream_enzyme_footprint, downstream_enzyme_footprint # used for TransitTime calculation
        
        self.extended_metrics = set(['TransitTime', 'LocalErrorRate', 'IABD'])
        if 'TransitTime' in self.pulse_tables or 'IABD' in self.pulse_tables:
            if 'IPD' not in self.pulse_tables: self.pulse_tables.append('IPD')
            if 'PulseWidth' not in self.pulse_tables: self.pulse_tables.append('PulseWidth')
        
        self.pulse_metric_normalizer = pulse_metric_normalizer
        self.aln_context_filter = aln_context_filter
        if template_length is None:
            self.contexts = CmpH5ColIterators.ColumnIterator(cmpH5, ref_info, ref_start=ref_start, ref_end=ref_end, pulse_tables=pulse_tables, target_coverage_ceiling=target_coverage_ceiling)
        else:
            self.contexts = CmpH5ColIterators.RolledColumnIterator(cmpH5, ref_info, template_length, pulse_tables=pulse_tables)

    def __iter__(self):
        for ref_pos, aln_ctxts in self.contexts:
            # note: aln_ctxts is an array of tuples (AlignmentHit, position)
            if self.aln_context_filter:
                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+' and self.aln_context_filter.filter(c)],
                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-' and self.aln_context_filter.filter(c)]}
            else:
                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+'],
                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-']}
            metrics = {}
            for metric in self.pulse_tables:
                metrics[metric] = {}
                for strand in '+', '-':
                    metrics[metric][strand] = numpy.array([self._getMetricValue(metric, aln_hit, aln_pos) for aln_hit, aln_pos in aln_ctxts_per_strand[strand]])
                
            yield ref_pos, metrics

#    def __iter__(self):
#        return self.pulseKineticsIterator()
#    
#    def iterateWithMoleculeIDs(self):
#        ''' Same as the default iterator, but also reports vectors of MoleculeIDs as the third element in the tuple, organized the same way as metrics.
#        '''
#        return self.pulseKineticsIterator(reportMoleculeIDs = True)
#    
#    def pulseKineticsIterator(self, reportMoleculeIDs = False):
#        for ref_pos, aln_ctxts in self.contexts:
#            # note: aln_ctxts is an array of tuples (AlignmentHit, position)
#            if self.aln_context_filter:
#                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+' and self.aln_context_filter.filter(c)],
#                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-' and self.aln_context_filter.filter(c)]}
#            else:
#                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+'],
#                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-']}
#            metrics = {}; molecule_ids = {}
#            for metric in self.pulse_tables:
#                metrics[metric] = {}
#                for strand in '+', '-':
#                    metrics[metric][strand] = numpy.array([self._getMetricValue(metric, aln_hit, aln_pos) for aln_hit, aln_pos in aln_ctxts_per_strand[strand]])
#                    if reportMoleculeIDs:
#                        molecule_ids[metric][strand] = numpy.array([aln_hit.moleculeId for aln_hit, aln_pos in aln_ctxts_per_strand[strand]])
#            
#            if reportMoleculeIDs:
#                yield ref_pos, metrics, molecule_ids
#            else:
#                yield ref_pos, metrics

    def iterateMedianPulseKinetics(self, reportMoleculeIDs = False):
        ''' Returns (ref_pos, metric_values, counts, [molecule_ids])
        '''
        pass
    
    def _getMetricValue(self, metric, aln_hit, aln_pos):
        if metric in self.extended_metrics and metric not in aln_hit.pulseInfo:
            self._loadExtendedMetricIntoAlnHit(metric, aln_hit)

        value = aln_hit[metric][aln_pos]
        
        if self.pulse_metric_normalizer:
            # TODO: make me smarter
            value *= self.pulse_metric_normalizer.getPolRateForAlnID(aln_hit.alignmentId)
        
        return value
    
    def _loadExtendedMetricIntoAlnHit(self, metric, aln_hit):
        # TODO: move this to dmtk.io.cmph5.CmpH5AlnHit
        aln_hit.pulseInfo[metric] = numpy.empty(aln_hit.alignedLength)
        aln_hit.pulseInfo[metric].fill(numpy.nan)
        lo_offset, hi_offset = self.upstream_enzyme_footprint, self.downstream_enzyme_footprint
        if aln_hit.target_strand == '-': # CHECK ME
            lo_offset, hi_offset = hi_offset, lo_offset
        if metric == 'TransitTime':
            for rel_ref_pos in xrange(lo_offset, aln_hit.target_end-aln_hit.target_start-hi_offset):
                my_col = aln_hit._target2aln[rel_ref_pos]
                upstream_col = aln_hit._target2aln[rel_ref_pos-lo_offset]
                downstream_col = aln_hit._target2aln[rel_ref_pos+hi_offset]
                if upstream_col > downstream_col: upstream_col, downstream_col = downstream_col, upstream_col
                aln_hit.pulseInfo['TransitTime'][my_col] = nansum([nansum(aln_hit.pulseInfo['IPD'][upstream_col:downstream_col]),
                                                                   nansum(aln_hit.pulseInfo['PulseWidth'][upstream_col:downstream_col])])
        elif metric == 'LocalErrorRate': # CHECK ME
            for rel_ref_pos in xrange(lo_offset, aln_hit.target_end-aln_hit.target_start-hi_offset):
                my_col = aln_hit._target2aln[rel_ref_pos]
                upstream_col = aln_hit._target2aln[rel_ref_pos-lo_offset]
                downstream_col = aln_hit._target2aln[rel_ref_pos+hi_offset]
                if upstream_col > downstream_col: upstream_col, downstream_col = downstream_col, upstream_col
                nmatches = float(sum([aln_hit.alignedQuery[i] == aln_hit.alignedTarget[i] for i in range(upstream_col, downstream_col)]))
                
                aln_hit.pulseInfo['LocalErrorRate'][my_col] = 1 - (nmatches / (lo_offset + hi_offset + 1))
        elif metric == 'IABD': # CHECK ME
            offset = 1 if aln_hit.target_strand == '-' else -1
            for rel_ref_pos in xrange(aln_hit.target_end-aln_hit.target_start):
                my_col = aln_hit._target2aln[rel_ref_pos]
                upstream_col = aln_hit._target2aln[rel_ref_pos-offset]
                ipd_sum = nansum([nansum(aln_hit.pulseInfo['IPD'][upstream_col+1:my_col+1]),
                                  nansum(aln_hit.pulseInfo['PulseWidth'][upstream_col+1:my_col])])
                aln_hit.pulseInfo['IABD'][my_col] = ipd_sum
    
    def iterateByZMW(self):
        ''' Same as the default iterator, but the dict is keyed by metric, then strand, then zmw (used for non-positional modifications).
        '''
        raise NotImplementedError()
    
    # TODO: implement average_subreads
    # TODO: implement accuracy requirement rules ("anchor=1", -1..0, ...)
    def iterateWithContext(self, upstream_context_len=2, downstream_context_len=5, average_subreads=False,
                           report_alignment_ids=False, report_zmw_ids=False,
                           required_motif=None):
        ''' Returns a KineticContextSet (contains fields: reference position, strand, sequence context, metrics).
        Metrics is a dict keyed by metric name, containing matrices of metric values.
        Each row in a matrix corresponds to an alignment (subread), and each column corresponds
        to an offset from the position of interest. Optionally, returns a vector
        of alignment IDs that the observations (rows) came from as the fifth element in the tuple.
        required_motif, if set, must be a tuple ('...XX', 'Y', 'ZZ...') describing the sequence motif before, at, and after the position of interest.
        '''
        if average_subreads: report_zmw_ids = True
        
        def _hasFlanks(aln_ctxt):
            aln_hit, aln_pos = aln_ctxt
            flank1 = abs(aln_hit.alnToTargetPos[aln_pos] - aln_hit.alnToTargetPos[0])
            flank2 = abs(aln_hit.alnToTargetPos[aln_pos] - aln_hit.alnToTargetPos[-1])
            return min(flank1, flank2) >= max(upstream_context_len, downstream_context_len)*2
        
        def _getContextSequence(aln_hit, ref_pos, upstream_l, downstream_l, direction):
            ctxt_start_in_ref = ref_pos - (upstream_l * direction)
            ctxt_end_in_ref = ref_pos + (downstream_l * direction)
            aln_pos_lo = aln_hit.targetToAlnPos[ctxt_start_in_ref-aln_hit.target_start]
            aln_pos_hi = aln_hit.targetToAlnPos[ctxt_end_in_ref-aln_hit.target_start]
            return aln_hit.alignedTarget[aln_pos_lo:aln_pos_hi+1].replace('-', '')
        
        class MotifMismatch(Exception):
            pass
        
        for ref_pos, aln_ctxts in self.contexts:
            # note: aln_ctxts is an array of tuples (AlignmentHit, position)
            if self.aln_context_filter:
                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+' and self.aln_context_filter.filter(c) and _hasFlanks(c)],
                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-' and self.aln_context_filter.filter(c) and _hasFlanks(c)]}
            else:
                aln_ctxts_per_strand = {'+': [c for c in aln_ctxts if c[0].target_strand == '+' and _hasFlanks(c)],
                                        '-': [c for c in aln_ctxts if c[0].target_strand == '-' and _hasFlanks(c)]}

            for strand in '+', '-':
                try:
                    metrics = {}; alignment_ids = []; zmw_ids = []
                    subread_counts = None
                    direction = 1 if strand is '+' else -1
                    if len(aln_ctxts_per_strand[strand]) == 0: continue
                    aln_ctxt, extended_aln_ctxt, motif_match = None, None, False
                    
                    for metric in self.pulse_tables:
                        metrics[metric] = numpy.empty([len(aln_ctxts_per_strand[strand]), upstream_context_len+downstream_context_len+1])
                    for aln_idx in xrange(len(aln_ctxts_per_strand[strand])):
                        aln_hit, aln_pos = aln_ctxts_per_strand[strand][aln_idx]
                        # WARNING: ref_pos in rolled template reference refers to the position in the template, not in the rolled reference.
                        # To properly index into the alignment, replace ref_pos with one from the alignment's index.
                        ref_pos_from_aln = aln_hit.alnToTargetPos[aln_pos]
    
                        if required_motif and not motif_match:
                            my_motif = _getContextSequence(aln_hit, ref_pos_from_aln, len(required_motif[0]), len(required_motif[2]), direction)
                            
                            if my_motif.upper() != (required_motif[0]+required_motif[1]+required_motif[2]).upper():
                                raise MotifMismatch
                            else:
                                motif_match = True
                        
                        ctxt_start_in_ref = ref_pos_from_aln - (upstream_context_len * direction)
                        for p in xrange(upstream_context_len+downstream_context_len+1):
                            my_ref_pos = ctxt_start_in_ref + (p * direction)
                            my_aln_pos = aln_hit.targetToAlnPos[my_ref_pos-aln_hit.target_start]
                            for metric in self.pulse_tables:
                                metrics[metric][aln_idx, p] = aln_hit[metric][my_aln_pos]

                        if report_alignment_ids:
                            alignment_ids.append(aln_hit.query_id)
                        if report_zmw_ids:
                            zmw_ids.append(aln_hit.moleculeId)

                        if aln_ctxt is None:
                            aln_ctxt = _getContextSequence(aln_hit, ref_pos_from_aln, upstream_context_len, downstream_context_len, direction)
                            extended_aln_ctxt = _getContextSequence(aln_hit, ref_pos_from_aln, upstream_context_len*2, downstream_context_len*2, direction)
                    
                    # TODO: support early data reduction by specifying minimum # subreads here
                    if average_subreads:
                        unique_zmw_ids = list(set(zmw_ids))
                        zmw_ids = numpy.array(zmw_ids)
                        subread_averaged_metrics = {}
                        subread_counts = numpy.empty(len(unique_zmw_ids), dtype=int)
                        for metric in self.pulse_tables:
                            subread_averaged_metrics[metric] = numpy.empty([len(unique_zmw_ids), upstream_context_len+downstream_context_len+1])
                            for i in range(len(unique_zmw_ids)):
                                subread_rows = metrics[metric][zmw_ids==unique_zmw_ids[i], :]
                                subread_averaged_metrics[metric][i, :] = numpy.ma.masked_invalid(subread_rows).mean(axis=0)
                                subread_counts[i] = subread_rows.shape[0]
                        metrics = subread_averaged_metrics
                        zmw_ids = unique_zmw_ids
                        # FIXME: fix alignment_ids here
                        alignment_ids = None

                    kcs = KineticContextSet(ref_pos, strand, aln_ctxt, extended_aln_ctxt, metrics)
                    
                    if report_alignment_ids: kcs.alignment_ids = alignment_ids
                    if report_zmw_ids: kcs.zmw_ids = zmw_ids
                    if average_subreads: kcs.subread_counts = subread_counts
                    
                    yield kcs
                except MotifMismatch:
                    pass
