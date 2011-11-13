"""
PulseKineticsNormalizer: Functions for normalizing SMRT sequencing kinetic data,
e.g. compensation for per-movie or per-subread variation

Author: Andrey Kislyuk
"""

import os, sys, itertools
import numpy, scipy, scipy.stats, random
from numpy import mean, median, std, var, sum, nansum, sqrt, sort, log

import logging

pol_rate_cache = {}
subread_phase_cache = {}

def getPolRate(aln_hit):
    ''' Compute, cache and return per-subread polymerization rate.
    '''
    # use start time if available: (aln_hit['StartTime'][-1] - aln_hit['StartTime'][0])
    if aln_hit.target_id not in pol_rate_cache:
        pol_rate_cache[aln_hit.target_id] = (aln_hit.target_end - aln_hit.target_start) / (sum(aln_hit['IPD']) + sum(aln_hit['PulseWidth']))
        assert(pol_rate_cache[aln_hit.target_id] > 0)
    return pol_rate_cache[aln_hit.target_id]


class AlignedPulseKineticsNormalizer:
    ''' Normalizes kinetic parameters by advance rate, local polymerization rate,
    or global metric average in the corresponding movie or zmw to account for intractable biases.
    '''
    def __init__(self, cmpH5, normalizable_metrics=['PulseWidth', 'IPD']):
        self.cmpH5 = cmpH5
        self.normalizable_metrics = normalizable_metrics
        self.aln_phases = None
        self.aln_pol_rates = None
        self.movie_pol_rates = None
        self.zmw_pol_rates = None
        self.aln_subread_ids = None

        self.HARD_MAX_IPD = 64.0
        self.MIN_ALIGNED_BASES_PER_MOVIE = 1e6
    
    def _cacheSubreadKinetics(self):
        ''' Compute subread phasing estimates.
        Phase 0 corresponds to slow/odd/forward strand half-laps around the smrtbell, phase 1 is fast/even/reverse strand half-laps.
        None means phase could not be determined.
        '''
        logging.debug("Caching subread kinetics in %s" % (self.cmpH5.name))
        self.aln_phases, self.aln_pol_rates, self.movie_pol_rates, self.zmw_pol_rates = {}, {}, {}, {}
        self.aln_subread_ids = {}
        my_records = self.cmpH5.asRecArray(["MovieID", "HoleNumber", "rStart", "AlnID", "AlnGroupID", "tStart", "tEnd", "Offset_begin", "Offset_end", "RCRefStrand"])
        my_records.sort(order=["MovieID", "HoleNumber", "rStart"])
        aln_group_cache = {}
        for movie_id, records_this_movie in itertools.groupby(my_records, lambda(x): x["MovieID"]):
            time_in_alns_this_movie = {0: 0.0, 1: 0.0}
            bases_in_alns_this_movie = {0: 0, 1: 0}
            self.zmw_pol_rates.setdefault(movie_id, {})
            for hole_number, aln_records_this_hole in itertools.groupby(records_this_movie, lambda(x): x["HoleNumber"]):
                subread_id = 0; phase = 0
                time_in_alns_this_zmw = {0: [], 1: []}
                bases_in_alns_this_zmw = {0: [], 1: []}
                aln_ids_this_zmw = {0: [], 1: []}
                for aln in aln_records_this_hole:
                    my_aln_group_id, my_aln_id = aln["AlnGroupID"], aln["AlnID"]
                    if my_aln_group_id not in aln_group_cache:
                        aln_group_cache[my_aln_group_id] = self.cmpH5.alnGroupById(my_aln_group_id)
                    my_aln_group = aln_group_cache[my_aln_group_id]
                    
                    my_ipds = my_aln_group['IPD'][aln["Offset_begin"]:aln["Offset_end"]]
                    my_pws = my_aln_group['PulseWidth'][aln["Offset_begin"]:aln["Offset_end"]]
                    
                    my_ipds[my_ipds > self.HARD_MAX_IPD] = numpy.nan
                    
                    time_in_aln = nansum(my_ipds) + nansum(my_pws)
                    my_num_bases = aln["tEnd"] - aln["tStart"]
                    
                    time_in_alns_this_zmw[phase].append(time_in_aln)
                    bases_in_alns_this_zmw[phase].append(my_num_bases)
                    aln_ids_this_zmw[phase].append(my_aln_id)
                    
                    self.aln_pol_rates[my_aln_id] = my_num_bases / time_in_aln
                    self.aln_subread_ids[my_aln_id] = subread_id

                    subread_id += 1; phase = subread_id % 2
                rate0 = sum(bases_in_alns_this_zmw[0]) / sum(time_in_alns_this_zmw[0])
                rate1 = sum(bases_in_alns_this_zmw[1]) / sum(time_in_alns_this_zmw[1])
                if rate0 > rate1:
                    bases_in_alns_this_zmw[0], bases_in_alns_this_zmw[1] = bases_in_alns_this_zmw[1], bases_in_alns_this_zmw[0]
                    time_in_alns_this_zmw[0], time_in_alns_this_zmw[1] = time_in_alns_this_zmw[1], time_in_alns_this_zmw[0]
                    aln_ids_this_zmw[0], aln_ids_this_zmw[1] = aln_ids_this_zmw[1], aln_ids_this_zmw[0]
                
                self.zmw_pol_rates[movie_id][hole_number] = {}
                for phase in 0, 1:
                    self.zmw_pol_rates[movie_id][hole_number][phase] = sum(bases_in_alns_this_zmw[phase]) / sum(time_in_alns_this_zmw[phase])
                    for aln_id in aln_ids_this_zmw[phase]:
                        self.aln_phases[aln_id] = phase
                    time_in_alns_this_movie[phase] += sum(time_in_alns_this_zmw[phase])
                    bases_in_alns_this_movie[phase] += sum(bases_in_alns_this_zmw[phase])
                
            self.movie_pol_rates[movie_id] = {}
            for phase in 0, 1:
                self.movie_pol_rates[movie_id][phase] = bases_in_alns_this_movie[phase] / time_in_alns_this_movie[phase]

    def getPhaseForAlnID(self, aln_id):
        if self.aln_phases is None:
            self._cacheSubreadKinetics()
        return self.aln_phases[aln_id]

    def getPolRateForAlnID(self, aln_id):
        if self.aln_pol_rates is None:
            self._cacheSubreadKinetics()
        return self.aln_pol_rates[aln_id]

    def getSubreadIDForAlnID(self, aln_id):
        ''' Here, subread ID refers to the index of the subread within the read.
        '''
        if self.aln_subread_ids is None:
            self._cacheSubreadKinetics()
        return self.aln_subread_ids[aln_id]
    
    def _cacheMoviePolRate(self, movie_name, cmpH5):
        total_time, total_num_bases = 0.0, 0
        total_time_per_zmw, total_num_bases_per_zmw = {}, {}
        #metric_values = dict(zip(self.normalizable_metrics, ([] for i in self.normalizable_metrics)))
        for ref_group in cmpH5.refGroupIterator():
            for aln_group in ref_group.alnGroupIterator():
                # Assumes each aln_group contains alignments from only one movie
                aln_rows = aln_group._myAlignmentIndexRows()
                aln_group_movie_name = cmpH5["/MovieInfo"].asDict("ID", "Name")[ int(cmpH5["/AlnInfo"].asRecArray()[aln_rows[0]]["MovieID"]) ]
                if movie_name != aln_group_movie_name: continue
                for aln_row in aln_rows:
                    alnInfoRow = cmpH5["/AlnInfo"].asRecArray( )[ aln_row ]
                    
                    aln_movie_name = cmpH5["/MovieInfo"].asDict("ID","Name",cache=True)[ int(alnInfoRow["MovieID"]) ]
                    assert(aln_movie_name == movie_name)
                    
                    my_zmw = alnInfoRow["HoleNumber"]
                    total_time_per_zmw.setdefault(my_zmw, 0.0)
                    total_num_bases_per_zmw.setdefault(my_zmw, 0)
                    my_ipds = aln_group['IPD'][alnInfoRow["Offset_begin"]:alnInfoRow["Offset_end"]]
                    my_pws = aln_group['PulseWidth'][alnInfoRow["Offset_begin"]:alnInfoRow["Offset_end"]]
                    
                    my_ipds[my_ipds > self.HARD_MAX_IPD] = self.HARD_MAX_IPD
                    
                    time_in_aln = nansum(my_ipds) + nansum(my_pws)
                    total_time += time_in_aln
                    total_time_per_zmw[my_zmw] += time_in_aln
                    
                    total_num_bases += (alnInfoRow["tEnd"] - alnInfoRow["tStart"])
                    total_num_bases_per_zmw[my_zmw] += (alnInfoRow["tEnd"] - alnInfoRow["tStart"])
                    
                    # TODO: fit for global speedup (linear f(startTime))
                    #time_in_aln = aln_group["StartTime"][alnInfoRow["Offset_end"]] - aln_group["StartTime"][alnInfoRow["Offset_begin"]]
        if total_num_bases < self.MIN_ALIGNED_BASES_PER_MOVIE:
            raise ValueError("Insufficient data to normalize by movie")
        self.movie_pol_rates[movie_name] = total_num_bases / total_time
        self.zmw_pol_rates[movie_name] = dict([(zmw, total_num_bases_per_zmw[zmw]/total_time_per_zmw[zmw]) for zmw in total_time_per_zmw.keys()])
    
    def normalizeAlnHitByMoviePolRate(self, aln_hit, cmpH5):
        # TODO: this won't work with split cmp.h5s
        #  - may use an array of cmp.h5s and open once for reading, read in all tables, then close
        # Lazily compute global movie/zmw averages per movie
        # per-movie cache doesn't need to be GCd, per-zmw maybe does
        movie_name = aln_hit.query_id.split('/')[0]
        if movie_name not in self.movie_pol_rates:
            self._cacheMoviePolRate(movie_name, cmpH5)
        for metric in self.normalizable_metrics:
            aln_hit.pulseInfo[metric+'_normalized'] = aln_hit.pulseInfo[metric] * self.movie_pol_rates[movie_name]

    def normalizeAlnHitByZMWPolRate(self, aln_hit, cmpH5):
        movie_name, zmw = aln_hit.query_id.split('/')[0:2]
        if movie_name not in self.movie_pol_rates:
            self._cacheMoviePolRate(movie_name, cmpH5)
        for metric in self.normalizable_metrics:
            aln_hit.pulseInfo[metric+'_normalized'] = aln_hit.pulseInfo[metric] * self.zmw_pol_rates[movie_name][zmw]
