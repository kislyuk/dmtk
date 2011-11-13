"""
PulseKinetics: Functions for loading and tabulating SMRT sequencing kinetic data,
and computing statistics for inferring DNA modifications

Author: Andrey Kislyuk
"""

import logging, cPickle, bz2
import numpy, scipy
from dmtk.io import cmph5
from dmtk.io.cmph5 import CmpH5ColIterators

def loadPulseStatsTable(filename):
    ''' Loads a pulse stats table from a compressed pickle. Expects bz2 compression.
    '''
    logging.info("Unpickling pulse stats table in %s" % (filename))
        
    in_fh = bz2.BZ2File(filename)
    pulse_stats_table = cPickle.load(in_fh)
    logging.info("Loaded %s" % (in_fh.name))
    in_fh.close()
    return pulse_stats_table

class ContextMatchingException(Exception):
    pass

class PulseStatsTable:
    """ Extract pulse metrics (IPD, PulseWidth, TransitTime, pkmid), normalize them, and emit a table
    of pulse statistics (mean, median, dispersion of each metric) and coverage associated with each sequence context.
    Optionally, append to a previously created table.
    """
    i2base = 'ATGCN'
    base2i = dict(zip(i2base, range(len(i2base))))
    base2i.update(dict(zip(i2base.lower(), range(len(i2base)))))

    def __init__(self, upstream_context_len=2, downstream_context_len=5, pulse_metrics=['IPD', 'PulseWidth']):
        self.upstream_context_len = upstream_context_len
        self.downstream_context_len = downstream_context_len
        self.context_len = upstream_context_len + downstream_context_len + 1
        self.pulse_metrics = pulse_metrics
        self.means = {}
        self.coverages = {}
        for metric in self.pulse_metrics:
            self.coverages[metric] = numpy.zeros(len(self.i2base)**self.context_len, dtype=int)
            self.means[metric] = numpy.zeros(len(self.i2base)**self.context_len)
            #self.means[metric] = numpy.empty(len(self.i2base)**self.context_len)
            #self.means[metric].fill(numpy.nan)

    def _hasFlanks(self, aln_hit, aln_pos):
        flank1 = abs(aln_hit.alnToTargetPos[aln_pos] - aln_hit.alnToTargetPos[0])
        flank2 = abs(aln_hit.alnToTargetPos[aln_pos] - aln_hit.alnToTargetPos[-1])
        return min(flank1, flank2) >= max(self.upstream_context_len, self.downstream_context_len)

    # TODO: median and dispersion don't work with O(1) space updating
    # consider streaming median approximation
    # TODO: quantify loss of precision from direct update of mean, if it's high, change to summation with subsequent division pass
    def addContextsFromCmpH5(self, cmpH5):
        logging.debug("Adding contexts from %s to pulse stats table" % (cmpH5.name))
        aln_hits_processed = 0
        for ref_info in cmpH5.refInfoIterator():
            for aln_hit in cmpH5.refGroupFromFullName(ref_info.fullName).alnHitIterator():
                assert(aln_hit.query_strand == '+')
                aln_hit.alignedQuery = aln_hit.alignedQuery.upper()
                aln_hit.alignedTarget = aln_hit.alignedTarget.upper()

                for aln_pos in range(aln_hit.alignedLength):
                    if aln_hit.alignedTarget[aln_pos] == '-': continue
                    
                    if not self._hasFlanks(aln_hit, aln_pos): continue
                    
                    ref_pos = aln_hit.alnToTargetPos[aln_pos]

                    ctxt_begin_ref_pos = ref_pos - (self.upstream_context_len if aln_hit.target_strand is '+' else self.downstream_context_len)
                    ctxt_end_ref_pos = ref_pos + (self.downstream_context_len if aln_hit.target_strand is '+' else self.upstream_context_len)
                    
                    aln_pos_lo = aln_hit.targetToAlnPos[ctxt_begin_ref_pos-aln_hit.target_start]
                    aln_pos_hi = aln_hit.targetToAlnPos[ctxt_end_ref_pos-aln_hit.target_start]
                    if aln_hit.target_strand is '-': aln_pos_lo, aln_pos_hi = aln_pos_hi, aln_pos_lo
                    seq_context = aln_hit.alignedTarget[aln_pos_lo:aln_pos_hi+1].replace('-', '')

                    assert(len(seq_context) <= self.context_len)
                    if len(seq_context) < self.context_len: continue

                    # TODO: computing full marginals is hard, so this interim solution sweeps wildcards from either side only
                    # FIXME: this does not properly trim data

                    matching_seq_contexts = [seq_context]
                    for i in range(1, self.upstream_context_len+1):
                        if seq_context[i+1] == 'N': continue
                        matching_seq_contexts.append('N'*i + seq_context[i:])
                    
                    for i in range(1, self.downstream_context_len+1):
                        if seq_context[-i] == 'N': continue
                        matching_seq_contexts.append(seq_context[:-i] + 'N'*i)
                    
                    for metric in self.pulse_metrics:
                        for c in matching_seq_contexts:
                            idx = self.context2idx(c)
                            self.coverages[metric][idx] += 1
                            v = aln_hit[metric][aln_pos]
                            if numpy.isnan(v): continue
                            
                            # FIXME: this does not properly trim data
                            if metric is 'IPD' and v > 20: continue

                            cov = self.coverages[metric][idx]
                            #if numpy.isnan(self.means[metric][idx]): self.means[metric][idx] = 0.0
                            self.means[metric][idx] = (cov * self.means[metric][idx] + v) / (cov+1)
                aln_hits_processed += 1
                if aln_hits_processed % 1000 == 0:
                    progress = 100 * float(aln_hits_processed) / cmpH5.numAlnHits
                    logging.debug("%.2f%% done" % (progress))
            
        logging.debug("Added contexts from %s to pulse stats table" % (cmpH5.name))

    def addContexts(self, contexts):
        c_lo, c_hi = self.context_extent
        
        for context in contexts:
            assert(len(context.template_upstream_nt_context) >= -c_lo)
            assert(len(context.template_downstream_nt_context) >= c_hi)
            seq_context = context.template_upstream_nt_context[c_lo:] + context.template_nt + context.template_downstream_nt_context[:c_hi]

            # TODO: computing full marginals is hard, so this interim solution sweeps wildcards from either side only
            matching_seq_contexts = [seq_context]
            for i in range(1, -c_lo+1):
                if seq_context[i+1] == 'N': continue
                matching_seq_contexts.append('N'*i + seq_context[i:])
            
            for i in range(1, c_hi+1):
                if seq_context[-i] == 'N': continue
                matching_seq_contexts.append(seq_context[:-i] + 'N'*i)
            
            for metric in self.pulse_metrics:
                v = context.metrics[metric]
                if numpy.isnan(v): continue
                
                for c in matching_seq_contexts:
                    idx = self.context2idx(c)

                    cov = self.coverages[metric][idx]
                    self.means[metric][idx] = (cov * self.means[metric][idx] + v) / (cov+1)
                    self.coverages[metric][idx] += 1

    def getMetricForContext(self, metric, context):
        i = self.context2idx(context)
        return self.means[metric][i], self.coverages[metric][i]

    # FIXME: the current sweeping approach is biased. Implement proper smoothing
    def getMetricForClosestContext(self, metric, seq_context, min_coverage=10):
        ''' Finds the metric value for the closest context satisfying the minimum coverage
        '''
        
        i = self.context2idx(seq_context)
        
        if self.coverages[metric][i] > min_coverage:
            return self.means[metric][i], self.coverages[metric][i], seq_context
        else:
            for i in range(1, self.upstream_context_len+1):
                matching_context = 'N'*i + seq_context[i:]
                j = self.context2idx(matching_context)
                if self.coverages[metric][j] > min_coverage:
                    return self.means[metric][j], self.coverages[metric][j], matching_context
            for i in range(1, self.downstream_context_len+1):
                matching_context = seq_context[:-i] + 'N'*i
                j = self.context2idx(matching_context)
                if self.coverages[metric][j] > min_coverage:
                    return self.means[metric][j], self.coverages[metric][j], matching_context
            # TODO: this will never be invoked until sweeping from both sides is done further up
            upstream_masked_seq_context = 'N'*self.upstream_context_len + seq_context[self.upstream_context_len:]
            for i in range(1, self.downstream_context_len+1):
                matching_context = upstream_masked_seq_context[:-i] + 'N'*i
                j = self.context2idx(matching_context)
                if self.coverages[metric][j] > min_coverage:
                    return self.means[metric][j], self.coverages[metric][j], matching_context
        raise ContextMatchingException
            
    def __str__(self):
        s = self.__class__.__name__ + ":\n"
        for i in range(len(self.coverages.values()[0])):
            s += self.idx2context(i)+': '
            for metric in self.pulse_metrics:
                s += "%s: c=%d, u=%.4f; " % (metric, self.coverages[metric][i], self.means[metric][i])
            s += "\n"
        return s

    def context2idx(self, c):
        assert(len(c)==self.context_len)
        i = 0
        L = len(self.i2base)
        for x in range(0, len(c)):
            i += self.base2i[c[-1-x]] * (L ** x)
        return i
    
    def idx2context(self, i):
        c = ''
        for x in range(self.context_len-1, -1, -1):
            stride = len(self.i2base) ** x
            c += self.i2base[i / stride]
            i %= stride
        return c
    
    def save(self, filename):
        ''' Saves this table to a compressed pickle using bz2 compression.
        Suitable for loading with PulseStatsTable.loadPulseStatsTable()
        '''
        out_fh = bz2.BZ2File(filename, 'w')
        cPickle.dump(self, out_fh, cPickle.HIGHEST_PROTOCOL)
        out_fh.close()
        logging.debug("Wrote table to file %s" % (filename))
