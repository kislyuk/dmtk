#!/usr/bin/env python

import sys, os, csv, tempfile, zipfile, xml.dom.minidom, re, subprocess, shutil, logging, bz2
import numpy, random, h5py
from scipy.stats import ks_2samp, ttest_ind, ttest_rel
from numpy import mean, median, std, sum, sqrt, sort, log, floor
from matplotlib import pyplot
from optparse import OptionParser
from glob import glob
from multiprocessing import Pool, cpu_count, current_process
# import shlex 

logging.basicConfig(level=logging.DEBUG)

from dmtk.model.kinetics.PulseKineticsUtils import AlignedPulseKineticsCollator, KineticContextStatistics, StandardAlnContextFilter, llr2samp

from dmtk.io import cmph5, GffIO, FastaIO, ConsensusIO
from dmtk.io.cmph5 import CmpH5ColIterators
from dmtk.model.SeqUtils import revcomp
from dmtk.io.ReferenceEntry import ReferenceEntry

#from dmtk.smrtpipe.SmrtPipeContext import SmrtPipeContext
#class spc_opts: debug = False; D = []
#TEMPDIR = SmrtPipeContext(options=spc_opts).config["TMP"]
#tempfile.tempdir = TEMPDIR

SMRTPIPE_JOB_GLOB = "/mnt/secondary/Smrtpipe/martin/8081/static/analysisJob/%d/%06d/data/aligned_reads.cmp.h5"

class ComparePulseData:
    """Command line tool used to pull out pulse metrics relevant to DNA modification detection studies."""

    def __init__(self):
        self.tempfiles = []
        self.__parseArgs()

    def __parseArgs(self):
        """Handle command line args"""
        usage = "Usage: %prog [options] modified[.cmp.h5|SMRTpipe JobID] control[.cmp.h5|SMRTpipe JobID]"
        parser = OptionParser(usage=usage, description=self.__class__.__doc__)

        parser.add_option("-t", "--template_length", help="Template Length", type="int")
        parser.add_option("-m", "--metrics", help="Pulse Metrics (comma-separated)", default='IPD,PulseWidth')
        parser.add_option("-f", "--modified_positions_fwd", help="Modified Positions (forward strand, comma-separated)")
        parser.add_option("-b", "--modified_positions_rev", help="Modified Positions (reverse strand, comma-separated)")
        parser.add_option("-l", "--min_length", help="Minimum aligned read length to use", type="int", default=50)
        parser.add_option("-z", "--min_z_score", help="Minimum z-score required to admit read", type="float", default=None)
        parser.add_option("-a", "--anchor", help="Number of alignment positions to either side we require to match perfectly", type="int", default=0)
        parser.add_option("-o", "--output_dir", help="Directory to write output files to", default=self.__class__.__name__+".out")
        parser.add_option("-c", "--coverage_ceiling", help="Artificially limit coverage to this number", type="int")
        parser.add_option("-s", "--strand", help="Restrict analysis to reads which aligned to this reference strand")
        parser.add_option("-y", "--ylim", help="Set y-axis limit on ratio plots. Default: None", type="float")
        parser.add_option("--modified_movies", help="Restrict modified dataset analysis to reads from these movies (comma-separated)")
        parser.add_option("--control_movies", help="Restrict control dataset analysis to reads from these movies (comma-separated)")
        parser.add_option("--ref_name", help="Process alignments to this reference")
        parser.add_option("--modified_template_name", help="Admit only alignments to this template as modified data")
        parser.add_option("--control_template_name", help="Admit only alignments to this template as control data")
        parser.add_option("--reference_repository", help="Reference repository to get sequence from", default="/mnt/secondary/Smrtpipe/repository")
        parser.add_option("--window_start", help="Start of window of interest in reference coordinates. Produce results only for this window", type="int", default=0)
        parser.add_option("--window_end", help="End of window of interest in reference coordinates. Produce results only for this window", type="int")
        parser.add_option("--write_zipfile", help="Compress output into a zip file. Default: False", action='store_true', default=False)
        parser.add_option("--motifs", help="Motifs of interest (comma-separated)")
        parser.add_option("--pos_in_motif", help="Position to highlight in motif(s)", type="int", default=0)
        parser.add_option("--exp_trim", help="Exponential distribution trim", type="float", default=0.05)
        parser.add_option("--make_report", help="Produce report with plots. Default: False", action='store_true', default=False)
        parser.add_option("--min_coverage", help="Minimum coverage to admit data (pass-through option for reports only)")
        parser.add_option("--target_coverage_ceiling", type="int", default=20000)
        parser.add_option("--save_plot_csv", help="Save CSV files for plots. Default: False", action='store_true', default=False)
        parser.add_option("--normalize_kinetics", help="Use kinetics normalization", action='store_true', default=False)
        parser.add_option("--plot_ipdr_histograms", help="Plot IPD ratio histograms for positions of interest", action='store_true', default=False)
        parser.add_option("--num_threads", help="Number of CPU threads to use", type="int", default=max(int(cpu_count()/4), 1))
        parser.add_option("--tempdir", default=tempfile.gettempdir())
        parser.add_option("--debug", action='store_true', default=False)
        (self.opts, args) = parser.parse_args()
        
        if self.opts.debug: logging.basicConfig(level=logging.DEBUG)
        
        if len(args) not in range(2, 3):
            parser.print_help()
            parser.error("Expected 2 arguments, got %d" % len(args))

        self.input = {'Modified': args[0], 'Control': args[1]}
        self.smrtpipe_job_names = {}

        if self.opts.ref_name:
            self.opts.modified_template_name = self.opts.control_template_name = self.opts.ref_name
        self.restrict_ref_id = {'Modified': self.opts.modified_template_name, 'Control': self.opts.control_template_name}
        
        self.opts.modified_positions_fwd = [ int(p) for p in self.opts.modified_positions_fwd.split(",") ] if self.opts.modified_positions_fwd != None else []
        self.opts.modified_positions_rev = [ int(p) for p in self.opts.modified_positions_rev.split(",") ] if self.opts.modified_positions_rev != None else []
        
        self.opts.metrics = self.opts.metrics.split(',')
        
        self.opts.motifs = self.opts.motifs.split(',') if self.opts.motifs else []

        self.opts.movies = {'Modified': self.opts.modified_movies.split(',') if self.opts.modified_movies else None,
                            'Control': self.opts.control_movies.split(',') if self.opts.control_movies else None}

#        for metric in self.opts.metrics:
#            if metric not in dmtk.io.CmpH5IO.pulseTables:
#                raise StandardError('Unknown pulse metric '+metric+' specified. Known metrics are: '+", ".join(dmtk.io.CmpH5IO.pulseTables))
        self.cmp_filenames = {}; self.ref_infos = {}
        self.consensus_calls = None
        self.reference_sequence = None

    def run(self):
        for condition in 'Modified', 'Control':
            cmph5_file, smrtpipe_job_name, reference_location = self._findCmpH5(self.input[condition])
            
            cmp_fh = h5py.File(cmph5_file, 'r')
            version = str(cmp_fh.attrs["Version"]) if "Version" in cmp_fh.attrs else None
            cmp_fh.close()
            if version == "PB_v0.1":
                cmph5_file = self._convertCmpTo12(cmph5_file)
            
            self.cmp_filenames[condition] = cmph5_file
            self.smrtpipe_job_names[condition] = smrtpipe_job_name

            cmp_fh = cmph5.factory.create(self.cmp_filenames[condition])
            for ref_info in cmp_fh.refInfoIterator():
                if self.restrict_ref_id[condition] and ref_info.fullName != self.restrict_ref_id[condition]:
                    continue
                self.ref_infos[condition] = ref_info; break
            
            if condition == 'Modified':
                self.reference_sequence = self._getReferenceSequence(reference_location)
                ref_name = self.ref_infos[condition].fullName
                ref_group = cmp_fh.refGroupFromFullName(ref_name)
                variants_gff_file = os.path.join(os.path.dirname(self.cmp_filenames[condition]), "variants.gff")
                if not os.path.exists(variants_gff_file):
                    variants_gff_file += ".gz"
                
                if ref_group.hasConsensus and os.path.exists(variants_gff_file):
                    variants_fh = GffIO.GffReader(variants_gff_file)
                    mapper = ConsensusIO.ConsensusReferenceMapper(cmp_fh, variants_fh, ref_names=[ref_name], gap_value=-1)
                    self.consensus_calls = mapper.consensusInRefSpace(self.ref_infos[condition])
                else:
                    self.consensus_calls = 'N'

            cmp_fh.close()
            
            if self.restrict_ref_id[condition] and condition not in self.ref_infos:
                raise StandardError('No data found matching reference name %s for condition %s' % (self.restrict_ref_id[condition], condition))

        if self.opts.window_end is None:
            self.opts.window_end = self.ref_infos['Modified'].length
            if self.opts.template_length:
                self.opts.window_end = self.opts.template_length
        
        if self.opts.window_end < self.opts.window_start:
            raise StandardError('Window start exceeds window end') 

        if self.opts.motifs:
            my_sequence = self.reference_sequence if self.reference_sequence != None else self.consensus_calls
            if my_sequence != None:
                for motif in self.opts.motifs:
                    fwd_motif_pos = [ m.start() + self.opts.pos_in_motif for m in re.finditer(motif, my_sequence, re.IGNORECASE) ]
                    self.opts.modified_positions_fwd.extend(fwd_motif_pos)
                    rev_motif_pos = [ len(my_sequence) - 1 - m.start() - self.opts.pos_in_motif for m in re.finditer(motif, revcomp(my_sequence), re.IGNORECASE) ]
                    self.opts.modified_positions_rev.extend(rev_motif_pos)
            else:
                sys.exit("Motif specified but no reference or consensus sequence found to find motifs in")
        
        if self.opts.write_zipfile:
            output_file = self.opts.output_dir
            self.opts.output_dir = tempfile.mkdtemp()
        else:
            if not os.path.isdir(self.opts.output_dir): os.mkdir(self.opts.output_dir)
        
        stride = 40000
        jobs = []
        for window_start in range(self.opts.window_start, self.opts.window_end+1, stride):
            jobs.append([window_start, min(window_start+stride-1, self.opts.window_end),
                         self.cmp_filenames, self.ref_infos, self.reference_sequence, self.consensus_calls, self.opts])
        
        if self.opts.num_threads > 1:
            p = Pool(self.opts.num_threads, maxtasksperchild=1)
            job_results = p.map_async(_ComparePulseDataWorker, jobs).get(99999999)
            try:
                p.close()
            except AttributeError:
                pass
        else:
            job_results = map(_ComparePulseDataWorker, jobs)
        
        worker_table_files = dict(zip(self.opts.metrics, ([] for i in self.opts.metrics)))
        for fileset in job_results:
            for metric, file in fileset.iteritems():
                worker_table_files[metric].append(file)

        table_files = {}
        for metric in self.opts.metrics:
#            table_files[metric] = tempfile.mkstemp(prefix=metric+"_", suffix='.csv.bz2')[1]
            table_files[metric] = os.path.join(self.opts.output_dir, metric+"_per_ref_pos.csv.bz2")
            with bz2.BZ2File(table_files[metric], 'w') as concat_csv_fh:
                header = bz2.BZ2File(worker_table_files[metric][0], 'r').readline()
                concat_csv_fh.write(header)
                for worker_csv in worker_table_files[metric]:
                    with bz2.BZ2File(worker_csv, 'r') as worker_csv_fh:
                        header = worker_csv_fh.readline()
                        for line in worker_csv_fh:
                            concat_csv_fh.write(line)
        
        if self.opts.make_report:
            child_process_handles = []
            for metric in self.opts.metrics:
                invoke_string = "makePulseKineticsReport.py '%s'" % table_files[metric]
                invoke_string += " --metric=%s" % metric
                invoke_string += " --title='%s (%s) vs. %s (%s)'" % (self.input['Modified'], self.smrtpipe_job_names['Modified'], self.input['Control'], self.smrtpipe_job_names['Control'])
                invoke_string += " --ref_name='%s'" % self.ref_infos['Modified'].fullName
                invoke_string += " --output='%s'" % self.opts.output_dir
                invoke_string += " --tempdir='%s'" % self.opts.tempdir
                if self.opts.min_coverage:
                    invoke_string += " --min_coverage=%s" % self.opts.min_coverage
                if self.opts.template_length: # proxy for LIMS template, which should be visualized with inverted strand sense
                    invoke_string += " --invert_strand_sense"
                if self.opts.save_plot_csv:
                    invoke_string += " --save_plot_csv"
                if self.opts.ylim is not None:
                    invoke_string += " --ymax=%f" % self.opts.ylim
                if self.opts.debug:
                    invoke_string += " --debug"
                logging.debug("Running %s" % invoke_string)
                ph = subprocess.Popen(invoke_string, shell=True)
                child_process_handles.append(ph)
            
            for ph in child_process_handles:
                if ph:
                    returncode = ph.wait()
                    if returncode: sys.exit(returncode)
        
        if self.opts.write_zipfile:
            z = zipfile.ZipFile(output_file, 'w')
            for f in os.listdir(self.opts.output_dir): z.write(os.path.join(self.opts.output_dir, f), f)
            z.close()

    def _findCmpH5(self, job_id_or_cmp_fn):
        full_cmp_fn = None; smrtpipe_job_title = ''; reference_location = None
        
        # If input is a string of digits, interpret as job id
        if re.search("^\d+$", job_id_or_cmp_fn):
            SMRTPIPE_JOB_GLOBS = ("/mnt/secondary/Smrtpipe/martin/8081/static/analysisJob/%d/%06d/data/aligned_reads.cmp.h5",
                                  "/mnt/secondary/Smrtanalysis/userdata/jobs/%03d/%06d/data/aligned_reads.cmp.h5")
            candidates = []
            for SMRTPIPE_JOB_GLOB in SMRTPIPE_JOB_GLOBS:
                candidates.extend(glob(SMRTPIPE_JOB_GLOB % (int(job_id_or_cmp_fn)/1000, int(job_id_or_cmp_fn))))
            if len(candidates) != 1:
                raise IOError, "Unable to locate unique cmp.h5 file for SMRTpipe Job %06d" % int(job_id_or_cmp_fn)
            job_id_or_cmp_fn = candidates[0]
        else:
            if not os.path.exists(job_id_or_cmp_fn) or job_id_or_cmp_fn[-7:] != ".cmp.h5":
                raise IOError, "The input specified (%s) does not appear to be a valid SMRTpipe Job ID or .cmp.h5 file." % str(job_id_or_cmp_fn)

        full_cmp_fn = os.path.abspath(job_id_or_cmp_fn)

        (toc_xml_fn, n_substitutions) = re.subn('data/.+\.cmp\.h5$', 'toc.xml', job_id_or_cmp_fn)
        if n_substitutions == 1 and os.path.exists(toc_xml_fn):
            my_metadata_xml = xml.dom.minidom.parse(toc_xml_fn)
            for xml_elt in my_metadata_xml.getElementsByTagName('attribute'):
                if xml_elt.getAttribute('name') == 'Job Name':
                    smrtpipe_job_title = xml_elt.childNodes[0].data.replace('\"', '')
            for xml_elt in my_metadata_xml.getElementsByTagName('apParam'):
                if xml_elt.getAttribute('name') == 'global.reference':
                    reference_location = xml_elt.childNodes[0].data
                    if not os.path.isabs(reference_location):
                        reference_location = os.path.join(os.path.dirname(full_cmp_fn), "..", reference_location)
        else:
            cmp_h = h5py.File(full_cmp_fn, 'r')
            if 'Repository' in cmp_h.attrs:
                reference_location = cmp_h.attrs['Repository'].replace('pbi://secondary/references', self.opts.reference_repository)
            else:
                logging.info("Could not get reference sequence location")
            cmp_h.close()
            
        return full_cmp_fn, smrtpipe_job_title, reference_location

    def _convertCmpTo12(self, input_filename, overwrite=False):
        logging.info('Converting cmp to 1.2 spec')
        if not (os.access(input_filename, os.W_OK) and overwrite):
            temp_cmph5 = tempfile.mkstemp(suffix='.cmp.h5')[1]
            logging.info("Creating working copy of file %s in %s" % (input_filename, temp_cmph5))
            shutil.copyfile(input_filename, temp_cmph5)
            input_filename = temp_cmph5
            self.tempfiles.append(temp_cmph5)
        
        invoke_string = "cmpH5Convert.py '%s' >/dev/null 2>&1" % input_filename
        logging.info("Running "+invoke_string)
        retcode = subprocess.call(invoke_string, shell=True)
        if retcode != 0: raise OSError("Error while converting file %s to 1.2 spec" % self.input_cmph5_file)
        return input_filename

    def _getReferenceSequence(self, reference_location, condition='Modified', sequence_name=None):
        reference_entry = ReferenceEntry(referenceDir = reference_location)
        if sequence_name == None: sequence_name = self.ref_infos[condition].fullName
        contig = reference_entry.getContig(sequence_name)
        if contig == None: return None
        for entry in FastaIO.SimpleFastaReader(contig.path):
            if entry.name == sequence_name or entry.name.split("|")[1] == sequence_name:
                return entry.sequence
        logging.info("Could not get reference sequence (no match to sequence name %s)" % sequence_name)
        return None

    def __del__(self):
#        for stat, value in self._stats.iteritems():
#            logging.info("%s(%s): %s = %d" % (self.__class__.__name__, self.input_cmph5_file, stat, value))
        for tf in self.tempfiles: os.unlink(tf)


def _ComparePulseDataWorker(args):
    window_start, window_end, cmp_filenames, ref_infos, reference_sequence, consensus_calls, opts = args
    logging.debug("_ComparePulseDataWorker started on %s:%d..%d in %s" % (ref_infos['Modified'].fullName, window_start, window_end, cmp_filenames['Modified']))
    
    current_process().daemon = False
    
    my_filters = {}
    for condition in 'Modified', 'Control':
        my_filters[condition] = StandardAlnContextFilter(min_alignment_length=opts.min_length,
                                                       min_alignment_zscore=opts.min_z_score,
                                                       anchor=opts.anchor,
                                                       restrict_to_movies = opts.movies[condition],
                                                       restrict_by_strand = opts.strand)
    
    pulse_kinetics_values = {'Modified': {}, 'Control': {}}
    kc_stats = {}
    
    for metric in opts.metrics:
        collator = AlignedPulseKineticsCollator([cmp_filenames['Modified'], cmp_filenames['Control']],
                                                [ref_infos['Modified'], ref_infos['Control']],
                                                ref_start=window_start,
                                                ref_end=window_end,
                                                pulse_metric=metric,
                                                template_length=opts.template_length,
                                                aln_context_filters=[my_filters['Modified'], my_filters['Control']],
                                                num_threads=opts.num_threads,
                                                normalize_kinetics=opts.normalize_kinetics,
                                                target_coverage_ceiling=opts.target_coverage_ceiling)
        for condition in 'Modified', 'Control':
            pulse_kinetics_values[condition][metric] = collator.getMetricValues(cmp_filenames[condition], ref_infos[condition])

    for condition in 'Modified', 'Control':
        kc_stats[condition] = KineticContextStatistics(pulse_kinetics_values[condition])
    
    worker_csv_out_filenames = {}
    for metric in opts.metrics:
        worker_csv_out_filenames[metric] = tempfile.mkstemp(prefix=metric+"_worker"+str(os.getpid())+"_", suffix='.csv.bz2')[1]
        csv_fh = bz2.BZ2File(worker_csv_out_filenames[metric], 'w')
        csv_writer = csv.writer(csv_fh)
        header = ['RefPosn', 'NucleotideFwd', 'NucleotideRev', 'RefNucleotideFwd', 'RefNucleotideRev',
                  'ModifiedMeanFwd', 'ModifiedMedianFwd', 'ModifiedStdevFwd', 'ModifiedSpreadLoFwd', 'ModifiedSpreadHiFwd', 'ModifiedCoverageFwd',
                  'ModifiedMeanRev', 'ModifiedMedianRev', 'ModifiedStdevRev', 'ModifiedSpreadLoRev', 'ModifiedSpreadHiRev', 'ModifiedCoverageRev',
                  'ControlMeanFwd', 'ControlMedianFwd', 'ControlStdevFwd', 'ControlSpreadLoFwd', 'ControlSpreadHiFwd', 'ControlCoverageFwd',
                  'ControlMeanRev', 'ControlMedianRev', 'ControlStdevRev', 'ControlSpreadLoRev', 'ControlSpreadHiRev', 'ControlCoverageRev',
                  'LLRFwd', 'LLRRev', 'PvalueFwd', 'PvalueRev', 'OfInterestFwd', 'OfInterestRev']
        
        csv_writer.writerow(header)
        row_buffer = []
        
        def _preprocess_row(row):
            for i in range(len(row)):
                if isinstance(row[i], float) and header[i] != 'PvalueFwd' and header[i] != 'PvalueRev':
                    row[i] = "%.8f" % row[i]
            return row
        
        use_se = True
        
        ctrl_posns = pulse_kinetics_values['Control'][metric].keys()
        for ref_pos in sorted(pulse_kinetics_values['Modified'][metric].keys()):
            if ref_pos not in ctrl_posns:
                # TODO: collect stats on matching vs. no-match positions
                #logging.info("Reference position %d in modified dataset but not in control dataset" % ref_pos)
                continue

            try:
                my_nt = consensus_calls[ref_pos].upper()
                my_nt_rc = revcomp(consensus_calls[ref_pos]).upper()
            except IndexError:
                my_nt, my_nt_rc = 'N', 'N'
            
            my_ref_nt, my_ref_nt_rc = 'N', 'N'
            if reference_sequence:
                my_ref_nt = reference_sequence[ref_pos].upper()
                my_ref_nt_rc = revcomp(reference_sequence[ref_pos]).upper()
            
            my_report_row = [str(ref_pos), my_nt, my_nt_rc, my_ref_nt, my_ref_nt_rc]
            
            my_stats = {'Modified': {'+': {}, '-': {}}, 'Control': {'+': {}, '-': {}}}
            for condition in 'Modified', 'Control':
                for strand in '+', '-':
                    my_stats[condition][strand]['mean'] = kc_stats[condition].meanMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand=strand)
                    my_stats[condition][strand]['median'] = kc_stats[condition].medianMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand=strand)
                    my_stats[condition][strand]['stdev'] = kc_stats[condition].stdevMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand=strand)
                    #(m2_ci_lo_fwd, m2_ci_hi_fwd) = kc_stats['Control'].bootstrapCIofMeanMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand='+')
                    if use_se:
#                        se = kc_stats[condition].bootstrapSEofMeanMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand=strand)
                        se = my_stats[condition][strand]['stdev'] / sqrt(len(kc_stats[condition]._getMetricDataForPos(ref_pos, metric, trim=0.0, coverage_ceiling=opts.coverage_ceiling, strand=strand)))
                        my_stats[condition][strand]['spread_lo'] = my_stats[condition][strand]['mean'] - se
                        my_stats[condition][strand]['spread_hi'] = my_stats[condition][strand]['mean'] + se
                    else:
                        my_stats[condition][strand]['spread_lo'] = my_stats[condition][strand]['mean'] - my_stats[condition][strand]['stdev']
                        my_stats[condition][strand]['spread_hi'] = my_stats[condition][strand]['mean'] + my_stats[condition][strand]['stdev']
                    my_stats[condition][strand]['coverage'] = len(kc_stats[condition]._getMetricDataForPos(ref_pos, metric, trim=0.0, coverage_ceiling=opts.coverage_ceiling, strand=strand))
                    my_report_row.extend([my_stats[condition][strand]['mean'],
                                          my_stats[condition][strand]['median'],
                                          my_stats[condition][strand]['stdev'],
                                          my_stats[condition][strand]['spread_lo'],
                                          my_stats[condition][strand]['spread_hi'],
                                          my_stats[condition][strand]['coverage']])
                    
            my_values = {'Modified': {'+': {}, '-': {}}, 'Control': {'+': {}, '-': {}}}
            my_logvalues = {'Modified': {'+': {}, '-': {}}, 'Control': {'+': {}, '-': {}}}
            for condition in 'Modified', 'Control':
                for strand in '+', '-':
                    my_values[condition][strand] = kc_stats[condition]._getMetricDataForPos(ref_pos, metric, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand=strand)
                    my_logvalues[condition][strand] = log(my_values[condition][strand] + 0.001)
                        
            llr_fwd = llr2samp(my_values['Modified']['+'], my_values['Control']['+'], min_samples=10)
            llr_rev = llr2samp(my_values['Modified']['-'], my_values['Control']['-'], min_samples=10)
            
#            tvalue_fwd, pvalue_fwd = ttest_ind(my_logvalues['Modified']['+'], my_logvalues['Control']['+'])
#            tvalue_rev, pvalue_rev = ttest_ind(my_logvalues['Modified']['-'], my_logvalues['Control']['-'])

            tvalue_fwd, pvalue_fwd, tvalue_rev, pvalue_rev = numpy.nan, numpy.nan, numpy.nan, numpy.nan
            if len(my_logvalues['Modified']['+']) > 0 and len(my_logvalues['Control']['+']) > 0:
                tvalue_fwd, pvalue_fwd = ks_2samp(my_logvalues['Modified']['+'], my_logvalues['Control']['+'])
                
#                bg_tvalues, bg_pvalues = [], []
#                for i in range(10):
#                    s = random.sample(my_logvalues['Control']['+'], len(my_logvalues['Control']['+']))
#                    rD, rp = ks_2samp(s[:len(s)/2], s[len(s)/2:])
#                    bg_tvalues.append(rD)
#                    bg_pvalues.append(rp)
#                
#                print pvalue_fwd, mean(bg_pvalues)
                
                
            if len(my_logvalues['Modified']['-']) > 0 and len(my_logvalues['Control']['-']) > 0:
                tvalue_rev, pvalue_rev = ks_2samp(my_logvalues['Modified']['-'], my_logvalues['Control']['-'])
            of_interest_fwd = 1 if ref_pos in opts.modified_positions_fwd else 0
            of_interest_rev = 1 if ref_pos in opts.modified_positions_rev else 0
                        
            my_report_row.extend([llr_fwd, llr_rev, pvalue_fwd, pvalue_rev, of_interest_fwd, of_interest_rev])
            
            row_buffer.append(_preprocess_row(my_report_row)) 
            if len(row_buffer) > 1000:
                csv_writer.writerows(row_buffer)
                row_buffer = []
            
            if opts.plot_ipdr_histograms and metric == 'IPD' and of_interest_fwd:
                ctrl_mean = kc_stats['Control'].meanMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand='+')
                # subread averaging and min subread req/multiple line plots for diff. subread counts goes here
                values = my_values['Modified']['+'] / ctrl_mean
                plot_ipdr_histogram(ref_infos['Modified'].fullName, ref_pos, '+', metric, values, opts)
            if opts.plot_ipdr_histograms and metric == 'IPD' and of_interest_rev:
                ctrl_mean = kc_stats['Control'].meanMetric(metric, ref_positions=ref_pos, trim=opts.exp_trim, coverage_ceiling=opts.coverage_ceiling, strand='-')
                values = my_values['Modified']['-'] / ctrl_mean
                plot_ipdr_histogram(ref_infos['Modified'].fullName, ref_pos, '-', metric, values, opts)
        
        csv_writer.writerows(row_buffer)
        csv_fh.close()
    return worker_csv_out_filenames

def plot_ipdr_histogram(ref_name, ref_pos, strand, metric, values, opts):
    """ Assumes values vector is sorted
    """
    pyplot.figure()
    pyplot.hist(values, numpy.arange(0, 20, 0.1))
    #pyplot.plot(numpy.cumsum(values) / numpy.sum(values))
    # pyplot.axis('equal'); pyplot.xlim(0, 1); pyplot.ylim(0, 1)
    pyplot.title("%s Ratio, %s:%d (%s)" % (metric, ref_name, ref_pos, strand))
    pyplot.xlabel(metric); pyplot.ylabel("Count")
    pyplot.savefig(os.path.join(opts.output_dir, "poi_histogram_%s_%d_%s.pdf" % (ref_name, ref_pos, strand)), format='pdf')

'''
# TODO: fix power report
#            for condition in 'Modified', 'Control':
#                for context in kin_contexts[condition].iterateAllContexts(ref_start=window_start, ref_end=window_end):
#                    if window_start is not None and context.ref_pos < window_start: continue
#                    if window_end is not None and context.ref_pos > window_end: continue
#                    of_interest = 1 if context.ref_pos in self.modified_positions else 0
#                    report = [condition, str(context.movie_name)+'/'+str(context.hole_number), str(context.template_pos), str(of_interest), str(context.num_passes)]
#                    #report = [condition, context.read_name, str(context.template_pos), str(of_interest), str(context.num_passes)]
#                    #report = [condition, context.read_name, str(context.ref_pos), str(of_interest), str(context.num_passes), context.read_aln_string, context.template_aln_string]
#                    # report.append(str(len(context.dropped_pulses)))
#                    for metric in 'IPD', 'PulseWidth':
#                        if metric not in self.augmented_pulse_metrics: continue
#                        report.append("%.4f" % (context.metrics[metric]/context.pol_rate))
#                    for metric in self.augmented_pulse_metrics:
#                        report.append("%.4f" % (context.metrics[metric]))
#                    #print ",".join(report)
#                    reports.setdefault(context.ref_pos, [])
#                    reports[context.ref_pos].append(report)
#
#            for ref_pos in sorted(reports.keys()):
#                this_pos_reports = reports[ref_pos]
#                if self.opts.coverage_ceiling is not None and len(this_pos_reports) > self.opts.coverage_ceiling:
#                    this_pos_reports = random.sample(this_pos_reports, self.opts.coverage_ceiling) # without replacement
#                for report in this_pos_reports:
#                    csv_report_fh.write(",".join(report) + "\n")

    if len(self.modified_positions) > 0:
        llrs = {}; bg_llrs = []; mp_llrs = {}; llr_rs = {}
        for mp in self.modified_positions: mp_llrs[mp] = []
        for ref_pos in sorted(pulse_kinetics_values['Modified'][metric].keys()):
            if window_start is not None and ref_pos < window_start: continue
            if window_end is not None and ref_pos > window_end: continue
            s1 = kc_stats['Modified']._getMetricDataForPos(ref_pos, 'IPD', trim=0.0, coverage_ceiling=opts.coverage_ceiling)
            try:
                s2 = kc_stats['Control']._getMetricDataForPos(ref_pos, 'IPD', trim=0.0, coverage_ceiling=opts.coverage_ceiling)
            except KeyError:
                continue
            llr = llr2samp(s1, s2)
            if llr is None: continue
            llrs[ref_pos] = llr
            in_mp_range = False
            for mp in self.modified_positions:
                if ref_pos in range(mp-3, mp+9):
                    mp_llrs[mp].append(llrs[ref_pos])
                    in_mp_range = True
            if not in_mp_range: bg_llrs.append(llrs[ref_pos])
        
        if len(bg_llrs) > 0:
            for mp in self.modified_positions:
                if len(mp_llrs[mp]) > 0: llr_rs[mp] = mean(mp_llrs[mp]) / mean(bg_llrs)
            
            power_report_fh = open(os.path.join(opts.output_dir, 'detection_power_metrics.txt'), 'w')
            print >> power_report_fh, "Goodness: %.2f (version 1)" % mean(llr_rs.values())
            
            print >> power_report_fh, "%d positions, mean LLRr = %.2f, min=%.2f, max=%.2f, median=%.2f" % (len(llr_rs), mean(llr_rs.values()), min(llr_rs.values()), max(llr_rs.values()), median(llr_rs.values()))
            print >> power_report_fh, "%d positions with LLRr <= 1" % numpy.nansum(numpy.array(llr_rs.values()) <= 1)        
            print >> power_report_fh, "Mean ramped LLR:"
            print >> power_report_fh, "Mean background LLR: %.2f" % mean(bg_llrs)
            #coverages = kc_stats['Modified'].coverage_map.values(); coverages.extend(kc_stats['Control'].coverage_map.values())
            #print >> power_report_fh, "Coverage: mean=%.1f, median=%.1f, min=%.1f, max=%.1f" % (mean(coverages), median(coverages), min(coverages), max(coverages))
    
            print >> power_report_fh, "Position, Mean windowed LLR, LLRr"
            for mp in self.modified_positions:
                if len(mp_llrs[mp]) > 0:
                    print >> power_report_fh, ','.join([str(mp), "%.4f" % mean(mp_llrs[mp]), "%.4f" % llr_rs[mp]])
                else:
                    print >> power_report_fh, ','.join([str(numpy.nan), str(numpy.nan)])
            power_report_fh.close()
'''


if __name__ == "__main__":
    epd = ComparePulseData()
    epd.run()
#    import cProfile; cProfile.run('epd.run()', 'epdprof')
