"""
PulseKineticsUtils: Functions for loading and tabulating SMRT sequencing kinetic data,
and computing statistics for inferring DNA modifications

Author: Andrey Kislyuk
"""

import os, sys, logging, re, bz2, csv
import numpy, scipy, scipy.stats, random
from numpy import mean, median, std, var, sum, sqrt, sort, log, log10, nansum, nonzero, product, floor
from scipy.stats import ks_2samp, ttest_ind, ttest_rel
from numpy.random import randint
from multiprocessing import Pool, cpu_count

from dmtk.io import cmph5
from dmtk.io.cmph5 import CmpH5ColIterators
from dmtk.model.kinetics.PulseKineticsIterators import AlignedPulseKineticsIterator
from dmtk.model.kinetics.Normalization import AlignedPulseKineticsNormalizer


def llr2samp(sample1, sample2, min_samples=10):
    """ Compute maximum likelihood estimates of exponential distribution fit
    of two datasets separately and together. Return
    log P(sample1|fit1) + log P(sample2|fit2) - log P(sample1, sample2|common_fit)
    (large values suggest data came from a different distribution)
    If either dataset contains less than min_samples points, return None
    """
    if not isinstance(sample1, numpy.ndarray): sample1 = numpy.array(sample1)
    if not isinstance(sample2, numpy.ndarray): sample2 = numpy.array(sample2)
    
    if len(sample1) < min_samples or len(sample2) < min_samples:
        return None
    
    exp_mean1 = mean(sample1)
    exp_mean2 = mean(sample2)
    mixed_exp_mean = (exp_mean1 * len(sample1) + exp_mean2 * len(sample2))/(len(sample1)+len(sample2))
    
    loglik1 = sum((1/exp_mean1)*numpy.exp(-sample1/exp_mean1))
    loglik2 = sum((1/exp_mean2)*numpy.exp(-sample2/exp_mean2))
    
    mixed_loglik = sum((1/mixed_exp_mean)*numpy.exp(-sample1/mixed_exp_mean)) \
        + sum((1/mixed_exp_mean)*numpy.exp(-sample2/mixed_exp_mean))
    
    llr = loglik1 + loglik2 - mixed_loglik
    return llr

def slow_llr2samp(sample1, sample2, min_samples=10):
    # old version using map
    if len(sample1) < min_samples or len(sample2) < min_samples:
        return None
    
    exp_mean1 = mean(sample1)
    exp_pdf1 = lambda(x): (1/exp_mean1)*numpy.exp(-x/exp_mean1)
    exp_mean2 = mean(sample2)
    exp_pdf2 = lambda(x): (1/exp_mean2)*numpy.exp(-x/exp_mean2)
    mixed_exp_mean = (sum(sample1)+sum(sample2))/(len(sample1)+len(sample2))
    mixed_exp_pdf = lambda(x): (1/mixed_exp_mean)*numpy.exp(-x/mixed_exp_mean)
    
    loglik1 = sum(map(exp_pdf1, sample1))
    loglik2 = sum(map(exp_pdf2, sample2))
    mixed_loglik = sum(map(mixed_exp_pdf, sample1)) + sum(map(mixed_exp_pdf, sample2))
    llr = loglik1 + loglik2 - mixed_loglik
    return llr


class KineticContextStatistics:
    """ Extracts statistics such as mean and standard error for kinetic metrics,
    grouped per reference position.
    Expects a dict with pulse metric names as keys,
    outputs of PulseKineticsUtils.AlignedPulseKineticsCollator.getMetricValues(cmpH5_fn, ref_info) as values
    """
    def __init__(self, pulse_metric_collators):
        self.pulse_metrics = pulse_metric_collators
        
        self._memoized_metric_data = None
        self._memoized_metric_metadata = None
    
    def _getMetricDataForPos(self, ref_pos, metric, trim, coverage_ceiling, strand=None):
        if (ref_pos, metric, trim, coverage_ceiling, strand) == self._memoized_metric_metadata:
            return self._memoized_metric_data
        
        if strand is None:
            fwd_strand_data = self.pulse_metrics[metric][ref_pos]['+']
            rev_strand_data = self.pulse_metrics[metric][ref_pos]['-']
            if len(fwd_strand_data) == 0:
                data = rev_strand_data
            elif len(rev_strand_data) == 0:
                data = fwd_strand_data
            else:
                data = numpy.concatenate((fwd_strand_data, rev_strand_data))
        else:
            data = self.pulse_metrics[metric][ref_pos][strand]
        if numpy.rank(data) > 1: data = data.flatten()
        if trim > 0.0: data = self._trim_both(data, trim)
        if coverage_ceiling is not None and len(data) > int(coverage_ceiling):
            data = numpy.array(random.sample(data, int(coverage_ceiling))) # without replacement

        nanmask = numpy.isfinite(data)
        if sum(nanmask) > 0: data = data[nanmask]
        
        self._memoized_metric_metadata = (ref_pos, metric, trim, coverage_ceiling, strand)
        self._memoized_metric_data = data
        
        return data

    def _statForMetric(self, stat_fn, metric, trim, ref_position, coverage_ceiling, strand):
        if isinstance(ref_position, int):
            data = self._getMetricDataForPos(ref_position, metric, trim, coverage_ceiling, strand)
            return stat_fn(data)
        else:
            return self._iterableStatForMetric(stat_fn, metric, trim, ref_position, coverage_ceiling, strand)
    
    def _iterableStatForMetric(self, stat_fn, metric, trim, ref_positions, coverage_ceiling, strand):
        if ref_positions is None: ref_positions = sorted(self.metrics.keys())
        for ref_pos in ref_positions:
            data = self._getMetricDataForPos(ref_pos, metric, trim, coverage_ceiling, strand)
            yield (ref_pos, stat_fn(data))

    def minMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): min(x), metric, trim, ref_positions, coverage_ceiling, strand)
    
    def maxMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): max(x), metric, trim, ref_positions, coverage_ceiling, strand)

    def meanMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): mean(x), metric, trim, ref_positions, coverage_ceiling, strand)
    
    def medianMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): median(x), metric, trim, ref_positions, coverage_ceiling, strand)

    def stdevMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): std(x), metric, trim, ref_positions, coverage_ceiling, strand)

    def bootstrapSEofMeanMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): self._bootstrap_se_of_mean(x), metric, trim, ref_positions, coverage_ceiling, strand)

    def bootstrapCIofMeanMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        return self._statForMetric(lambda(x): self._bootstrap_ci_of_mean(x), metric, trim, ref_positions, coverage_ceiling, strand)
    
    def expOverdispersionOfMetric(self, metric, trim=0.0, ref_positions=None, coverage_ceiling=None, strand=None):
        ''' Assuming exponentiality, compute overdispersion (variance over mean).
        Significant deviation from 1 is evidence of departure from exponentiality '''
        return self._statForMetric(lambda(x): var(x)/mean(x), metric, trim, ref_positions, coverage_ceiling, strand)
    
    def _trim_both(self, my_list, trim_per_tail):
        """ Trim a fraction of data from both sides of the distribution of the
        data in the list, as sorted by the default sort.
        
        @param my_list: Data to trim
        @param trim_per_tail: Fraction of data to trim per tail
        """
        if my_list is None or len(my_list) == 0 or trim_per_tail == 0: return my_list
        trimmed_list = scipy.stats.trimboth(numpy.sort(my_list), trim_per_tail)
        return trimmed_list.flatten()
        
    def _bootstrap_ci_of_mean(self, x):
        if not isinstance(x, numpy.ndarray): x = numpy.array(x)
        if (len(x) < 1): return (numpy.nan, numpy.nan)
        idx = randint(len(x), size=(100, len(x)))
        bx = x[idx] # resamples x with replacement
        means = numpy.mean(bx, axis=1)
        #conf_interval = (means[25], means[975])
        conf_interval = (min(means), max(means))
        return conf_interval

    def _bootstrap_se_of_mean(self, x):
        if not isinstance(x, numpy.ndarray): x = numpy.array(x)
        if (len(x) < 1): return numpy.nan
        bins = 16
        #idx = randint(len(x), size=(bins, int(len(x)/bins)))
        idx = randint(len(x), size=(bins, len(x)))
        bx = x[idx] # resamples x with replacement
        means = numpy.mean(bx, axis=1)
        se = numpy.std(means) / numpy.sqrt(bins)
        return se


class ClassifierException(Exception):
    pass

class AlignedPulseKineticsCollator:
    ''' Collates a single pulse metric per reference position and strand. Uses a multithreaded pool of readers.
    Caution: loads all values of the specified pulse metric for the speicfied range into memory.
    cmpH5_filenames, ref_infos are lists of equal length.
    Example: to retrieve IPDs for a methylated/native sample and a conrol/WGA sample:
    cmpH5_filenames=["native.cmp.h5", "control.cmp.h5"], ref_infos=[ref_info1, ref_info2]
    There should be no readers open in this process on the given cmp.h5s.
    '''
    def __init__(self, cmpH5_filenames, ref_infos, ref_start=0, ref_end=None, pulse_metric='IPD', template_length=None, num_threads=None,
                 aln_context_filters=None, normalize_kinetics=True, target_coverage_ceiling=None):
        if num_threads is None: num_threads = cpu_count()
        if ref_end is None: ref_end = max(map(lambda(x): x.length, ref_infos))
        
        stride = max(int((ref_end-ref_start)/400)*100, 200)
        stride = min(stride, 100000)
        
        jobs = []
        for coords_start in range(ref_start, ref_end+1, stride):
            for i in range(len(cmpH5_filenames)):
                jobs.append([cmpH5_filenames[i], ref_infos[i], dict(pulse_metric=pulse_metric,
                                                                    ref_start=coords_start,
                                                                    ref_end=min(coords_start+stride-1, ref_end),
                                                                    template_length=template_length,
                                                                    aln_context_filter=aln_context_filters[i] if aln_context_filters else None,
                                                                    normalize_kinetics=normalize_kinetics,
                                                                    target_coverage_ceiling=target_coverage_ceiling)])
        if num_threads > 1:
            p = Pool(num_threads, maxtasksperchild=1)
            job_results = p.map_async(_AlignedPulseKineticsCollatorWorker, jobs).get(99999999)
        else:
            job_results = map(_AlignedPulseKineticsCollatorWorker, jobs)
    
        self.collated_metric_values = {}
        for cmpH5_filename, ref_info, metric_values in job_results:
            self.collated_metric_values.setdefault(cmpH5_filename, {})
            self.collated_metric_values[cmpH5_filename].setdefault(ref_info.fullName, {})
            for ref_pos, metric_values_this_ref_pos in metric_values.iteritems():
                # assert(ref_pos not in self.collated_metric_values[cmpH5_filename][ref_info.fullName])
                self.collated_metric_values[cmpH5_filename][ref_info.fullName][ref_pos] = {'+': metric_values_this_ref_pos['+'],
                                                                                           '-': metric_values_this_ref_pos['-']}
        if num_threads > 1:
            try:
                p.close()
            except AttributeError:
                pass
    
    def getMetricValues(self, cmpH5_filename, ref_info, ref_pos=None, strand=None):
        if ref_pos is None:
            return self.collated_metric_values[cmpH5_filename][ref_info.fullName]
        elif strand is None:
            return self.collated_metric_values[cmpH5_filename][ref_info.fullName][ref_pos]
        else:
            return self.collated_metric_values[cmpH5_filename][ref_info.fullName][ref_pos][strand]

def _AlignedPulseKineticsCollatorWorker(args):
    cmpH5_filename, ref_info, opts = args
    pulse_metric = opts['pulse_metric']
    metric_values = {}
    cmpH5 = cmph5.factory.create(cmpH5_filename)
    pkn = AlignedPulseKineticsNormalizer(cmpH5) if opts['normalize_kinetics'] else None
    for ref_pos, metrics in AlignedPulseKineticsIterator(cmpH5, ref_info, pulse_tables=[pulse_metric],
                                                         ref_start=opts['ref_start'], ref_end=opts['ref_end'],
                                                         template_length=opts['template_length'], aln_context_filter=opts['aln_context_filter'],
                                                         pulse_metric_normalizer=pkn, target_coverage_ceiling=opts['target_coverage_ceiling']):
        metric_values[ref_pos] = {}
        metric_values[ref_pos]['+'] = metrics[pulse_metric]['+']
        metric_values[ref_pos]['-'] = metrics[pulse_metric]['-']

    cmpH5.close()
    return cmpH5_filename, ref_info, metric_values


class StandardAlnContextFilter:
    def __init__(self, min_alignment_length = None,
                 min_alignment_zscore = None,
                 anchor = None,
                 upstream_anchor = None,
                 downstream_anchor = None,
                 restrict_to_movies = None,
                 restrict_by_strand = None):
        self.min_alignment_length = min_alignment_length
        self.min_alignment_zscore = min_alignment_zscore
        self.anchor = anchor
        self.upstream_anchor = upstream_anchor
        self.downstream_anchor = downstream_anchor
        self.stats = {}
        self.restrict_to_movies = restrict_to_movies
        self.restrict_by_strand = restrict_by_strand
        
    def filter(self, aln_context):
        aln_hit, aln_pos = aln_context
        if self.min_alignment_length is not None:
            if aln_hit.query_end - aln_hit.query_start + 1 < self.min_alignment_length:
                #self._stats['rejected_aln_hits_lo_length'] += 1
                return False
        
        if self.min_alignment_zscore is not None:
            if aln_hit.zScore < self.min_alignment_zscore:
                #self._stats['rejected_aln_hits_lo_zscore'] += 1
                return False
        
        if self.restrict_to_movies:
            movie_name, hole_number = aln_hit.query_id.split('/')[0:2]
            if movie_name not in self.restrict_to_movies:
                return False
        
        if self.restrict_by_strand:
            if aln_hit.target_strand != self.restrict_by_strand:
                return False
        
        if self.anchor == 0:
            if aln_hit.alignedQuery[aln_pos] != aln_hit.alignedTarget[aln_pos]:
                return False
        elif self.anchor > 0:
            if aln_hit.alignedQuery[aln_pos-self.anchor:aln_pos+self.anchor+1] != aln_hit.alignedTarget[aln_pos-self.anchor:aln_pos+self.anchor+1]:
                return False
            if len(aln_hit.alignedQuery[aln_pos-self.anchor:aln_pos+self.anchor+1]) != self.anchor*2 + 1:
                return False

        return True

class WhitelistAlnContextFilter(StandardAlnContextFilter):
    def __init__(self, zmw_whitelist=None, **kwargs):
        self.zmw_whitelist = zmw_whitelist
        self.last_query_id = None
        self.last_aln_context_accepted = False

        StandardAlnContextFilter.__init__(self, **kwargs)
        
    def filter(self, aln_context):
        aln_hit, aln_pos = aln_context
        
        if self.last_query_id != aln_hit.query_id:
            movie_name, hole_number = aln_hit.query_id.split('/')[0:2]
            if movie_name+'/'+hole_number in self.zmw_whitelist:
                self.last_aln_context_accepted = True
            else:
                self.last_aln_context_accepted = False
            self.last_query_id = aln_hit.query_id
        
        return StandardAlnContextFilter.filter(self, aln_context) if self.last_aln_context_accepted else False



# TODO: move me to dmtk.io.GffIO
class AggregatingFeatureWriter:
    ''' Aggregates adjacent GFF features into a single feature. Wraps around dmtk.io.GffIO.GffWriter '''
    def __init__(self, gff_writer):
        self.gff_queues = {'+': {}, '-': {}}
        self.gff_writer = gff_writer

    def addFeature(self, gff_feature):
        # assumes ordered input
        # if feature is adjacent to current queue's last element, append
        # else write aggregate feature, restart queue
        self.gff_queues[gff_feature.strand].setdefault(gff_feature.seqid, [])
        
        cur_queue = self.gff_queues[gff_feature.strand][gff_feature.seqid]
        if len(cur_queue) > 0 and cur_queue[-1].start != gff_feature.start - 1:
            self._write_queue(cur_queue)
            self.gff_queues[gff_feature.strand][gff_feature.seqid] = []
        
        self.gff_queues[gff_feature.strand][gff_feature.seqid].append(gff_feature)

    def _write_queue(self, queue):
        scores = numpy.array(map(lambda(x): float(x['score']), queue))
        top_scoring_feature_i = nonzero(scores == max(scores))[0][0]
        top_scoring_feature = queue[top_scoring_feature_i]
        
        top_scoring_feature.start = queue[0].start
        top_scoring_feature.end = queue[-1].end

        self.gff_writer.writeRecord(top_scoring_feature)

    def __del__(self):
        for strand, queues_this_strand in self.gff_queues.iteritems():
            map(self._write_queue, queues_this_strand.itervalues())


class AlignedPulseKineticsClassifier:
    ''' Abstract class for a classifier that operates on KineticContextSets.
    '''
    def __init__(self, pulse_stats_table):
        raise NotImplementedError()
    
    def classify(self, kinetic_context_set, control_kinetic_context_set=None):
        ''' Takes as input a kinetic context set (set of observations at an
        aligned reference position) and a kinetic context table (a data structure
        for expected observation values). Produces as output a single value
        representing a classification score.
        '''
        raise NotImplementedError()


TRIM_ABOVE_MEDIAN_MULTIPLE=64

def exp_loglik(sample, exp_mean, trim=0.05):
    assert(type(sample) == numpy.ndarray and sample.ndim == 1)
    # this is too slow - filter upstream and/or use median
#    sample = sample[numpy.isfinite(sample)]
#    trimmed_sample = sample[numpy.logical_and(sample > 0, sample < median(sample)*TRIM_ABOVE_MEDIAN_MULTIPLE)]
#    return log(1/exp_mean)*len(trimmed_sample) - numpy.sum(trimmed_sample)/exp_mean
    return log(1/exp_mean)*len(sample) - numpy.sum(sample)/exp_mean

def norm_exp_loglik(sample, exp_mean, trim=0.05):
    assert(type(sample) == numpy.ndarray and sample.ndim == 1)
    # this is too slow - filter upstream and/or use median
    #sample = sample[numpy.isfinite(sample)]
    #trimmed_sample = sample[numpy.logical_and(sample > 0, sample < median(sample)*TRIM_ABOVE_MEDIAN_MULTIPLE)]
    #return log(1/exp_mean) - numpy.sum(trimmed_sample)/(exp_mean*len(trimmed_sample))
    return log(1/exp_mean) - numpy.sum(sample)/(exp_mean*len(sample))


class scoringStats:
    """ Struct for reporting counters. """
    def __init__(self, **entries): self.__dict__.update(entries)
    def __str__(self): return self.__class__.__name__+': '+self.__dict__.__str__()


class WindowedRampedLLR(AlignedPulseKineticsClassifier):
    def __init__(self, pulse_stats_table):
        self.pulse_stats_table = pulse_stats_table
        
        ramp_lo, ramp_hi = 0.2, 1.0
        self.ramp = map(lambda(x): ramp_lo + (ramp_hi-ramp_lo)*x/(pulse_stats_table.upstream_context_len),
                        range(1, pulse_stats_table.upstream_context_len+1))
        self.ramp.append(1.0)
        self.ramp.extend(map(lambda(x): ramp_hi - (ramp_hi-ramp_lo)*x/(pulse_stats_table.downstream_context_len),
                             range(pulse_stats_table.downstream_context_len)))

        # TODO: support multiple, user profiles
        self.profile = numpy.ones(pulse_stats_table.context_len) * 4

        # For mC only: expect a drop immediately following position of interest
        #self.profile[pulse_stats_table.upstream_context_len+1] = 0.75
        
        self.context_len = pulse_stats_table.context_len
        self.scoring_stats = scoringStats(low_coverage_in_table=0,
                                          context_substitutions=0)
        self.min_lookup_table_coverage = 10
        self.min_read_coverage = 3
        self.min_subread_coverage = 3
        
        self.ipd_clip_threshold = 30.0
        self.ipd_discard_threshold = 60.0

    
    def classify0(self, kinetic_context_set):
        base_lls = []
        elev_lls = []
        bs_lls = []
        # TODO: expand to other metrics
        #for metric in kinetic_context_set.metrics.keys():
        
        metric = 'IPD'
        if kinetic_context_set.metrics[metric].shape[0] < self.min_read_coverage:
            return None
        # TODO: Evaluate each row independently, computing per-molecule measurements, for mixture estimation, etc
        
        extended_aln_ctxt = kinetic_context_set.extended_aln_ctxt
        for context_pos in range(self.context_len):
            subcontext = extended_aln_ctxt[context_pos:context_pos+self.context_len]
            if len(subcontext) < self.context_len:
                raise ClassifierException
            #expected_mean, ctxt_coverage = pulse_stats_table.getMetricForContext(metric, subcontext)
            expected_mean, ctxt_coverage, derived_context = self.pulse_stats_table.getMetricForClosestContext(metric, subcontext, self.min_lookup_table_coverage)
            if ctxt_coverage < self.min_lookup_table_coverage:
                self.scoring_stats.low_coverage_in_table += 1
                continue
            if derived_context != subcontext:
                #logging.debug("Substituted context %s for %s" % (derived_context, subcontext))
                self.scoring_stats.context_substitutions += 1

#                lls['base'].append(norm_exp_loglik(kinetic_context_set.metrics[metric][:, context_pos], expected_mean) * self.ramp[context_pos])
#                lls['elev'].append(norm_exp_loglik(kinetic_context_set.metrics[metric][:, context_pos], expected_mean*self.profile[context_pos]) * self.ramp[context_pos])
#                lls['bs'].append(norm_exp_loglik(kinetic_context_set.metrics[metric][:, context_pos], expected_mean*1000) * self.ramp[context_pos])

            my_metric_values = kinetic_context_set.metrics[metric][:, context_pos]
            # if var(my_metric_values)/mean(my_metric_values) > self.max_overdisp:
            
            
            
            # TODO: for performance, nans must be filtered out before this stage
            my_metric_values = my_metric_values[numpy.isfinite(my_metric_values)]
            my_metric_values = my_metric_values[my_metric_values < self.ipd_discard_threshold]
            my_metric_values[my_metric_values > self.ipd_clip_threshold] = self.ipd_clip_threshold
            my_sum = numpy.sum(my_metric_values)
#            base_lls.append((log(1/expected_mean)*len(my_metric_values) - my_sum/expected_mean) * self.ramp[context_pos])
#            elev_lls.append((log(1/(expected_mean*self.profile[context_pos]))*len(my_metric_values) - my_sum/(expected_mean*self.profile[context_pos])) * self.ramp[context_pos])

            base_lls.append((log(1/expected_mean) - my_sum/(expected_mean*len(my_metric_values))) * self.ramp[context_pos])
            elev_lls.append((log(1/(expected_mean*self.profile[context_pos])) - my_sum/(expected_mean*len(my_metric_values)*self.profile[context_pos])) * self.ramp[context_pos])
            
            #bs_factor = 100
            #bs_lls.append(log(1/(expected_mean*bs_factor*self.profile[context_pos]))*len(my_metric_values) - my_sum/(expected_mean*bs_factor*self.profile[context_pos]))
            
            
            #base_lls.append(exp_loglik(my_metric_values, expected_mean) * self.ramp[context_pos])
            #elev_lls.append(exp_loglik(my_metric_values, expected_mean*self.profile[context_pos]) * self.ramp[context_pos])
            
            
        #print >> sys.stderr, "Scored", kinetic_context_set
        #logging.debug("score = %f - %f = %f" % (numpy.nansum(elev_lls), numpy.nansum(base_lls), numpy.nansum(elev_lls) - numpy.nansum(base_lls)))
        #logging.debug("BS LL = %f" % numpy.nansum(bs_lls))
        #if numpy.nansum(bs_lls) > max(numpy.nansum(elev_lls), numpy.nansum(base_lls)):
            #logging.debug("triggered")
            #return None
        return numpy.nansum(elev_lls) - numpy.nansum(base_lls)

    
    def classify1(self, kinetic_context_set):
        ipdrs = []
        # TODO: expand to other metrics
        
        metric = 'IPD'
        if kinetic_context_set.metrics[metric].shape[0] < self.min_read_coverage:
            return None
        # TODO: Evaluate each row independently, computing per-molecule measurements, for mixture estimation, etc
        
        extended_aln_ctxt = kinetic_context_set.extended_aln_ctxt
        for context_pos in range(self.context_len):
            subcontext = extended_aln_ctxt[context_pos:context_pos+self.context_len]
            if len(subcontext) < self.context_len:
                raise ClassifierException
            #expected_mean, ctxt_coverage = pulse_stats_table.getMetricForContext(metric, subcontext)
            expected_mean, ctxt_coverage, derived_context = self.pulse_stats_table.getMetricForClosestContext(metric, subcontext, self.min_lookup_table_coverage)
            if ctxt_coverage < self.min_lookup_table_coverage:
                self.scoring_stats.low_coverage_in_table += 1
                continue
            if derived_context != subcontext:
                #logging.debug("Substituted context %s for %s" % (derived_context, subcontext))
                self.scoring_stats.context_substitutions += 1

            my_metric_values = kinetic_context_set.metrics[metric][:, context_pos]
            my_metric_values = my_metric_values[numpy.isfinite(my_metric_values)]
            my_median = median(my_metric_values)
            my_ipdr = my_median/expected_mean
            my_od = var(my_metric_values)/mean(my_metric_values)
            if my_od > 1:
                my_ipdr /= my_od
            #print len(my_metric_values), my_ipdr, my_od
            ipdrs.append(my_ipdr * self.ramp[context_pos])

        #print kinetic_context_set.ref_pos, kinetic_context_set.strand, ipdrs
        return numpy.nansum(ipdrs)
    
    def classify2(self, kinetic_context_set):
        metric = 'IPD'
        if kinetic_context_set.metrics[metric].shape[0] < self.min_read_coverage:
            return None
        
        return numpy.nansum(kinetic_context_set.metrics[metric])/kinetic_context_set.metrics[metric].shape[0]


    def classify3(self, kinetic_context_set):
        norm_ipds = []
        # TODO: expand to other metrics
        
        metric = 'IPD'
        if kinetic_context_set.metrics[metric].shape[0] < self.min_read_coverage:
            return None
        # TODO: Evaluate each row independently, computing per-molecule measurements, for mixture estimation, etc
        
        extended_aln_ctxt = kinetic_context_set.extended_aln_ctxt
        for context_pos in range(self.context_len):
            subcontext = extended_aln_ctxt[context_pos:context_pos+self.context_len]
            if len(subcontext) < self.context_len:
                raise ClassifierException
            #expected_mean, ctxt_coverage = pulse_stats_table.getMetricForContext(metric, subcontext)
            expected_mean, ctxt_coverage, derived_context = self.pulse_stats_table.getMetricForClosestContext(metric, subcontext, self.min_lookup_table_coverage)
            if ctxt_coverage < self.min_lookup_table_coverage:
                self.scoring_stats.low_coverage_in_table += 1
                continue
            if derived_context != subcontext:
                #logging.debug("Substituted context %s for %s" % (derived_context, subcontext))
                self.scoring_stats.context_substitutions += 1

            my_metric_values = kinetic_context_set.metrics[metric][:, context_pos]
            my_metric_values = my_metric_values[numpy.isfinite(my_metric_values)]
            my_metric_values = my_metric_values[my_metric_values < self.ipd_discard_threshold]
            
            if len(my_metric_values) < 1: return None
            
            my_metric_values[my_metric_values > self.ipd_clip_threshold] = self.ipd_clip_threshold
            my_metric_values /= expected_mean
            norm_ipds.append(numpy.nansum(my_metric_values) * self.ramp[context_pos] / len(my_metric_values))

        #print kinetic_context_set.ref_pos, kinetic_context_set.strand, ipdrs
        #return numpy.nansum(norm_ipds)
        return numpy.nansum(norm_ipds)

    def classify(self, kinetic_context_set):
        all_norm_ipds = []
        # TODO: expand to other metrics
        
        metric = 'IPD'
        if kinetic_context_set.metrics[metric].shape[0] < self.min_read_coverage:
            return None
        
        # FIXME: alignment_ids misnomer
        read_ids = [re.sub("/[^/]+$", "", query_id) for query_id in kinetic_context_set.alignment_ids]
        unique_read_ids = set(read_ids)
        read_ids = numpy.array(read_ids)
#        for movie_id, records_this_movie in itertools.groupby(my_records, lambda(x): x["MovieID"]):
        
        #read_cov = len(set([re.sub("/[^/]+$", "", query_id) for query_id in context_set.alignment_ids]))
        
        extended_aln_ctxt = kinetic_context_set.extended_aln_ctxt
        for read_id in unique_read_ids:
            norm_ipds = []
            for context_pos in range(self.context_len):
                subcontext = extended_aln_ctxt[context_pos:context_pos+self.context_len]
                if len(subcontext) < self.context_len:
                    raise ClassifierException
                #expected_mean, ctxt_coverage = pulse_stats_table.getMetricForContext(metric, subcontext)
                expected_mean, ctxt_coverage, derived_context = self.pulse_stats_table.getMetricForClosestContext(metric, subcontext, self.min_lookup_table_coverage)
                if ctxt_coverage < self.min_lookup_table_coverage:
                    self.scoring_stats.low_coverage_in_table += 1
                    continue
                if derived_context != subcontext:
                    #logging.debug("Substituted context %s for %s" % (derived_context, subcontext))
                    self.scoring_stats.context_substitutions += 1
    
                #my_metric_values = kinetic_context_set.metrics[metric][:, context_pos]
                my_metric_values = kinetic_context_set.metrics[metric][read_ids==read_id, context_pos]
                my_metric_values = my_metric_values[numpy.isfinite(my_metric_values)]
                my_metric_values = my_metric_values[my_metric_values < self.ipd_discard_threshold]
                
                if len(my_metric_values) < 1: return None
                
                my_metric_values[my_metric_values > self.ipd_clip_threshold] = self.ipd_clip_threshold
                my_metric_values /= expected_mean
                norm_ipds.append(numpy.nansum(my_metric_values) * self.ramp[context_pos] / len(my_metric_values))
            all_norm_ipds.append(numpy.nansum(norm_ipds))

        #print kinetic_context_set.ref_pos, kinetic_context_set.strand, ipdrs

#        multiread_attenuator = 0.25
#        return numpy.exp(numpy.nansum(numpy.log(all_norm_ipds)+numpy.log(multiread_attenuator)))
    
        return numpy.nansum(all_norm_ipds)


class CaseControlKSWithProfile(AlignedPulseKineticsClassifier):
    ''' Case-control classification using KS 2-sample test and context-specific mean ratios profile
    Assumes full enrichment of modified molecules in case sample (and 0 in control)
    '''
    def __init__(self, priors):
        self.min_coverage = 10
        
        '''
        table_fn="/home/UNIXHOME/akislyuk/projects/extractpulsedata/neb/m5C_pRRS_Sau3aim_51_compare/m5C_pRRS_Sau3aim_51_new/IPD_per_ref_pos.csv.bz2"

        with bz2.BZ2File(table_fn, 'r') as input_csv:
            csv_reader = csv.reader(input_csv)
            header = csv_reader.next()
            for line in csv_reader:
                ref_pos_data = dict(zip(header, line))
                ref_pos = int(ref_pos_data['RefPosn'])

                for strand in '+', '-':
                    pass

        '''
        
        self.expected_ratios_table = None
        
    
    def classify(self, kinetic_context_set, control_kinetic_context_set=None):
        metric='IPD'
        trim_per_tail = 0.05
        
        my_pvalues=numpy.empty(kinetic_context_set.metrics[metric].shape[1])
        my_pvalues.fill(numpy.nan)
        
        my_logratios=numpy.empty(kinetic_context_set.metrics[metric].shape[1])
        my_logratios.fill(numpy.nan)
        
        for context_pos in range(kinetic_context_set.metrics[metric].shape[1]):
            my_metric_values = kinetic_context_set.metrics[metric][:, context_pos]
            my_metric_values = scipy.stats.trimboth(numpy.sort(my_metric_values), trim_per_tail)
            my_log_metric_values = log(my_metric_values[numpy.isfinite(my_metric_values)] + 0.001)
            
            my_ctrl_metric_values = control_kinetic_context_set.metrics[metric][:, context_pos]
            my_ctrl_metric_values = scipy.stats.trimboth(numpy.sort(my_ctrl_metric_values), trim_per_tail)
            my_log_ctrl_metric_values = log(my_ctrl_metric_values[numpy.isfinite(my_ctrl_metric_values)] + 0.001)
            
            if len(my_log_metric_values) < 10 or len(my_log_ctrl_metric_values) < 10:
                #print "passed", kinetic_context_set.ref_pos, len(my_log_metric_values), len(my_log_ctrl_metric_values)
                return None
            
#            print "Case:", my_log_metric_values
#            print "Ctrl:", my_log_ctrl_metric_values
            tvalue, pvalue = ks_2samp(my_log_metric_values, my_log_ctrl_metric_values)
#            print kinetic_context_set.ref_pos, kinetic_context_set.strand, context_pos, pvalue
            my_pvalues[context_pos] = pvalue
            my_logratios[context_pos] = mean(my_log_metric_values)/mean(my_log_ctrl_metric_values)
        
        my_expected_logratios=numpy.empty(kinetic_context_set.metrics[metric].shape[1])
        my_expected_logratios.fill(0)
        
        profile_dist = sum((my_logratios - my_expected_logratios)**2)
        
        return sum(-log10(my_pvalues + 1e-300))
