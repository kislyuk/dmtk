#!/usr/bin/env python

import sys
import numpy
from optparse import OptionParser
from dmtk.io import cmph5
from dmtk.io.cmph5 import Basemap
from dmtk.io.FastaIO import FastaEntry, writeFastaEntry

usage="Usage: %prog input.cmp.h5 output.fasta"

parser = OptionParser(usage=usage)
(opts, args) = parser.parse_args()

if len(args) not in range(2, 3):
	parser.print_help()
	parser.error("Expected 2 arguments, got %d" % len(args))

(input_cmph5_fn, output_fasta_fn) = args

with open(output_fasta_fn, 'w') as fasta_fh:
	ch = cmph5.factory.create(input_cmph5_fn)
	for ref_info in ch.refInfoIterator():
		ref_group = ch.refGroupFromFullName(ref_info.fullName)
		if not ref_group.hasConsensus:
			sys.exit("Reference group for %s has no consensus calls" % ref_info.fullName)
		
		consensus_array = numpy.empty_like(ref_group.consensusCalls._dataset)
		ref_group.consensusCalls._dataset.read_direct(consensus_array)
		consensus = "".join(Basemap[consensus_array])
		writeFastaEntry(fasta_fh, FastaEntry(ref_info.fullName, consensus))
