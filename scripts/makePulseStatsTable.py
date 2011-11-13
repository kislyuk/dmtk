#!/usr/bin/env python

import sys, os, re, logging, cPickle, bz2
import numpy, scipy
from optparse import OptionParser

from dmtk.model.kinetics import PulseStatsTable
from dmtk.io import cmph5

logging.basicConfig(level=logging.DEBUG)


class makePulseStatsTable:
    """ Extract pulse metrics (IPD, PulseWidth, TransitTime, pkmid), normalize them, and emit a table
    of pulse statistics (mean, median, dispersion of each metric) and coverage associated with each sequence context.
    Optionally, append to a previously created table.
    """

    def __init__(self):
        self.__parseArgs()

    def __parseArgs(self):
        """Handle command line args"""
        usage = "Usage: %prog [options] input.cmp.h5 output_table.pickle.bz2"
        parser = OptionParser(usage=usage, description=self.__class__.__doc__)

        parser.add_option("-m", "--metrics", help="Pulse Metrics (comma-separated). Default: IPD, PulseWidth")
        parser.add_option("-u", "--upstream_context", help="Upstream context nucleotides to extract and report. Default: 0", type="int")
        parser.add_option("-d", "--downstream_context", help="Downstream context nucleotides to extract and report. Default: 0", type="int")
        parser.add_option("-t", "--previous_table", help="Pulse stats table to append to. Default: None")

        parser.set_defaults(metrics = 'IPD,PulseWidth',
                            upstream_context = 2,
                            downstream_context = 7,
                            previous_table = None)
        (self.opts, args) = parser.parse_args()

        if len(args) not in range(2, 3):
            parser.error("Expected 2 arguments")

        self.input, self.output = args
        
        if not os.path.isfile(self.input):
            raise IOError("File %s not found" % self.input)

        self.opts.metrics = self.opts.metrics.split(',')

        for metric in self.opts.metrics:
            if metric not in cmph5.pulseTables:
                raise StandardError('Unknown pulse metric '+metric+' specified. Known metrics are: '+", ".join(cmph5.pulseTables))
#        self.augmented_pulse_metrics = self.opts.metrics[:]
#        self.augmented_pulse_metrics.append('TransitTime')

    def run(self):
        if self.opts.previous_table:
            in_fh = open(self.opts.previous_table, 'rb')
            t = cPickle.load(in_fh)
            in_fh.close()
        else:
            t = PulseStatsTable.PulseStatsTable(upstream_context_len=self.opts.upstream_context,
                                                downstream_context_len=self.opts.downstream_context,
                                                pulse_metrics=self.opts.metrics)
        cmpH5 = cmph5.factory.create(self.input)
        t.addContextsFromCmpH5(cmpH5)
        
        t.save(self.output)


if __name__ == "__main__":
    r = makePulseStatsTable()
    r.run()
