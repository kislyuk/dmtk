"""
PulseKinetics: Functions for loading and tabulating SMRT sequencing kinetic data,
and computing statistics for inferring DNA modifications

Author: Andrey Kislyuk
"""

import os, sys, re
import numpy, scipy, scipy.stats, logging
from dmtk.io import cmph5
from dmtk.io.cmph5 import CmpH5ColIterators
from dmtk.model.kinetics import PulseKineticsUtils, PulseStatsTable
from dmtk.model.kinetics.PulseKineticsIterators import AlignedPulseKineticsIterator

default_pulse_stats_table = None
default_upstream_context_len = 2
default_downstream_context_len = 5

def scoreKineticEvents(cmpH5, ref_info, ref_start=0, ref_end=None, template_length=None,
                       pulse_stats_table=default_pulse_stats_table, required_motif=None,
                       ctrl_cmpH5=None):
    """ Score statistically significant excursions in pulse kinetics, indicating possible DNA modifications or other events.
    (Will support a plug-in framework for probability models)
    
    pulse_stats_table can be either a PulseStatsTable object or a string containing a compressed pickle filename.
    required_motif, if set, must be a tuple ('...XX', 'Y', 'ZZ...') describing the sequence motif before, at, and after the position of interest.
    
    If ctrl_cmph5 is specified, it is used instead of pulse_stats_table.
    """
    i = AlignedPulseKineticsIterator(cmpH5, ref_info, ref_start=ref_start, ref_end=ref_end, template_length=template_length)
    
    if type(pulse_stats_table) == str:
        pulse_stats_table = PulseStatsTable.loadPulseStatsTable(pulse_stats_table)
    
    if ctrl_cmpH5:
        ctrl_i = AlignedPulseKineticsIterator(ctrl_cmpH5, ref_info, ref_start=ref_start, ref_end=ref_end, template_length=template_length)
        ctrl_i = ctrl_i.iterateWithContext(upstream_context_len=default_upstream_context_len,
                                           downstream_context_len=default_downstream_context_len,
                                           required_motif=required_motif,
                                           report_alignment_ids=True)
        
        next_ctrl_context_set = ctrl_i.next()
        
        my_classifier = PulseKineticsUtils.CaseControlKSWithProfile(pulse_stats_table)

        for context_set in i.iterateWithContext(upstream_context_len=default_upstream_context_len,
                                                downstream_context_len=default_downstream_context_len,
                                                required_motif=required_motif,
                                                #report_alignment_ids=False):
                                                report_alignment_ids=True):
            try:
                while next_ctrl_context_set.ref_pos < context_set.ref_pos or next_ctrl_context_set.strand != context_set.strand:
                    next_ctrl_context_set = ctrl_i.next()
            except StopIteration:
                break

            try:
                score = my_classifier.classify(context_set, next_ctrl_context_set)
                if score is None: continue
            except (PulseKineticsUtils.ClassifierException, PulseStatsTable.ContextMatchingException):
                continue
            
            # FIXME: alignment_ids misnomer
            read_cov = len(set([re.sub("/[^/]+$", "", query_id) for query_id in context_set.alignment_ids]))        
            subread_cov = context_set.metrics['IPD'].shape[0]
            
            yield context_set.ref_pos, context_set.strand, score, read_cov, subread_cov
    else:
        my_classifier = PulseKineticsUtils.WindowedRampedLLR(pulse_stats_table)
    
        for context_set in i.iterateWithContext(upstream_context_len=pulse_stats_table.upstream_context_len,
                                                downstream_context_len=pulse_stats_table.downstream_context_len,
                                                required_motif=required_motif,
                                                #report_alignment_ids=False):
                                                report_alignment_ids=True):
            try:
                score = my_classifier.classify(context_set)
                if score is None: continue
            except (PulseKineticsUtils.ClassifierException, PulseStatsTable.ContextMatchingException):
                continue
            
            # FIXME: alignment_ids misnomer
            read_cov = len(set([re.sub("/[^/]+$", "", query_id) for query_id in context_set.alignment_ids]))        
            subread_cov = context_set.metrics['IPD'].shape[0]
            
            yield context_set.ref_pos, context_set.strand, score, read_cov, subread_cov

