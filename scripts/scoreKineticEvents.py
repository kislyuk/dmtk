#!/usr/bin/env python

import sys, os, re, logging, cPickle, bz2, copy, shutil
import numpy, scipy
from numpy import ceil, nonzero
from optparse import OptionParser

import time
from random import random

from multiprocessing import Pool, cpu_count

from dmtk.io import cmph5, GffIO
from dmtk.model.kinetics import PulseKinetics, PulseStatsTable

logging.basicConfig(level=logging.DEBUG)

def _gff_feature(seqName, start, end, strand, type = '.'):
    r = GffIO.Gff3Record(seqName = seqName)
    r.start = start; r.end = end; r.strand = strand; r.type = type
    return r

def ske_worker(worker_args):
    #logging.debug("worker %d accepted job: %s" % (os.getpid(), worker_args.__str__()))
    cmph5_fn, ref_info, ske_worker_pulse_stats_table, opts = worker_args
    
    cmpH5 = cmph5.factory.create(cmph5_fn)
    ctrl_cmpH5 = cmph5.factory.create(opts.ctrl) if opts.ctrl else None
    
    events_gff = GffIO.GffWriter(open(opts.gff_output, 'w'))
    
    agg_writer = aggregatingFeatureWriter(events_gff)
    
    logging.debug("Scoring events in %s:%s" % (cmpH5.name, ref_info.fullName))
    total_events = 0
    
    required_motif = None
    if opts.motif_at: required_motif = (opts.motif_upstream, opts.motif_at, opts.motif_downstream)

    try:
        for event_info in PulseKinetics.scoreKineticEvents(cmpH5, ref_info,
                                                           ref_start = opts.ref_start,
                                                           ref_end = opts.ref_end,
                                                           template_length = opts.template_length,
                                                           pulse_stats_table = ske_worker_pulse_stats_table,
                                                           ctrl_cmpH5 = ctrl_cmpH5,
                                                           required_motif = required_motif):
            ref_pos, strand, score, read_cov, subread_cov = event_info
            if opts.min_score is not None and score < opts.min_score: continue
            r = _gff_feature(seqName = ref_info.fullName, start = ref_pos+1, end = ref_pos+1, strand = strand, type = '.')
            r.put('score', "%.4f" % score)
            r.put('read_cov', "%d" % read_cov)
            r.put('subread_cov', "%d" % subread_cov)
            #r.put('cov', coverage_per_strand[reference][strand][pos])
            
            #events_gff.writeRecord(r)
            agg_writer.addFeature(r)
            total_events += 1
        
        logging.debug("Wrote %d events to %s" % (total_events, opts.gff_output))
    except IndexError:
        logging.debug("Ignoring IndexError")
        pass

    del agg_writer # write out all the queues
    events_gff.close()
    cmpH5.close()


class aggregatingFeatureWriter:
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


class scoreKineticEvents:
    """ Score statistically significant excursions in pulse kinetics, indicating possible DNA modifications or other events.
    """

    def __init__(self):
        self.__parseArgs()

    def __parseArgs(self):
        """Handle command line args"""
        usage = "Usage: %prog [options] input.cmp.h5 (--ctrl=ctrl_input.cmp.h5|--pulse_stats_table=table.pickle.bz2) --gff_output=out.gff"
        parser = OptionParser(usage=usage, description=self.__class__.__doc__)

        parser.add_option("--pulse_stats_table", help="Input pulse stats table file")
        parser.add_option("-c", "--ctrl", help="cmp.h5 file with control data aligned to the same reference")
        parser.add_option("--motif_upstream", help="Reference sequence motif to require before any position of interest")
        parser.add_option("--motif_at", help="Reference sequence motif to require at any position of interest")
        parser.add_option("--motif_downstream", help="Reference sequence motif to require after any position of interest")
        parser.add_option("-t", "--template_length", help="LIMS Template Length", type="int")
        parser.add_option("--reference_name", help="Only scan pulse data for this reference")
        parser.add_option("--ref_start", help="Start of window of interest in reference coordinates", type="int", default=0)
        parser.add_option("--ref_end", help="End of window of interest in reference coordinates", type="int")
        parser.add_option("-o", "--gff_output", help="Name of GFF file to generate, listing top events")
        parser.add_option("--algorithm", help="Name of detection algorithm to use")
        parser.add_option("--min_score", help="Minimum score reported by the algorithm to print to GFF file", type="float", default=80.0)
        parser.add_option("--chunk_size", help="Maximum number of bases to process at once", type="int", default=int(1e7))
        parser.add_option("--threads", help="Maximum number of threads to run", type="int", default=cpu_count())
        (self.opts, args) = parser.parse_args()

        if len(args) not in range(1, 2):
            parser.error("Expected 1 arguments")

        self.input_cmp = args[0]

        if not os.path.isfile(self.input_cmp):
            raise IOError("File %s not found" % self.input_cmp)
        
        if self.opts.gff_output is None:
            self.opts.gff_output = re.sub(".cmp.h5$", "_kinetic_events.gff", self.input_cmp)
        
        if self.opts.motif_upstream or self.opts.motif_at or self.opts.motif_downstream:
            if self.opts.motif_upstream is None or self.opts.motif_at is None or self.opts.motif_downstream is None:
                parser.error("Options motif_upstream, motif_at, motif_downstream must be supplied together")

    def run(self):
        cmpH5 = cmph5.factory.create(self.input_cmp)

        ske_worker_pulse_stats_table = PulseStatsTable.loadPulseStatsTable(self.opts.pulse_stats_table) if self.opts.pulse_stats_table else None 
        
        jobs = []; gff_shards = []
        for ref_info in cmpH5.refInfoIterator():
            if self.opts.reference_name and ref_info.fullName != self.opts.reference_name: continue
            if self.opts.ref_start or self.opts.ref_end:
                jobs.append([self.input_cmp, ref_info, ske_worker_pulse_stats_table, self.opts])
            else:
                # TODO: retain overlaps exceeding context length, then delete duplicates, to avoid loss of context at chunk boundaries
                for ref_start in range(0, ref_info.length, self.opts.chunk_size):
                    my_opts = copy.deepcopy(self.opts)
                    my_opts.ref_start = ref_start
                    my_opts.ref_end = ref_start+self.opts.chunk_size-1
                    my_opts.gff_output = re.sub(".gff$", "_%s_%d_%d.gff" % (ref_info.fullName, my_opts.ref_start, my_opts.ref_end), self.opts.gff_output)
                    jobs.append([self.input_cmp, ref_info, ske_worker_pulse_stats_table, my_opts])
                    gff_shards.append(my_opts.gff_output)

        logging.debug("Job manifest for %s: %d jobs" % (cmpH5.name, len(jobs)))
        cmpH5.close()


        if self.opts.threads > 1:
            p = Pool(self.opts.threads)
            p.map_async(ske_worker, jobs).get(99999999)
        else:
            map(ske_worker, jobs)

        
        full_gff = open(self.opts.gff_output, 'w')
        for gff_filename in gff_shards:
            shutil.copyfileobj(open(gff_filename, 'r'), full_gff)
        full_gff.close()
        
        map(os.unlink, gff_shards)
        
        logging.debug("Wrote job results to %s" % (self.opts.gff_output))


if __name__ == "__main__":
    r = scoreKineticEvents()
    r.run()
    #import cProfile; cProfile.run('r.run()', 'skeprof')
