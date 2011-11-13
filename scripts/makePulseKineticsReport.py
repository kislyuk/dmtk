#!/usr/bin/env python

import sys, os, tempfile, logging, subprocess, csv, gzip, bz2, shutil
import numpy
from numpy import float64, log10, product, isfinite, isinf, sqrt
from collections import deque

from string import Template
from optparse import OptionParser

from dmtk.plot.GraphItem import *
from dmtk.plot.RScript import *

from dmtk.io import GffIO, BedIO
from dmtk.model.kinetics import PulseKineticsUtils
from dmtk.io.BedIO import BedWriter, BedRecord

RSCRIPT_DEFAULT = "pulseKineticsCaseControlPlot.R"
R_BIN = 'R'

class MakePulseKineticsReport:
    """Command line tool used to visualize pulse metrics relevant to DNA modification detection studies."""
    def __init__(self):
        self.__parseOptions()
        
    def __parseOptions(self):
        usage = 'Usage: %prog [--help] [options] pulseKinetics.csv.bz2'
        parser = OptionParser(usage=usage, description=__doc__)
        
        parser.add_option("-o", "--output", help="Output directory for report files", default=".")
        parser.add_option("--script", help="Location of R script", default=RSCRIPT_DEFAULT)
        parser.add_option("--metric", help="Name of pulse metric being reported on")
        parser.add_option("--title", help="Title of plot")
        parser.add_option("--ref_name", help="Reference name")
        parser.add_option("--xmin", help="X minimum", default='NULL')
        parser.add_option("--xmax", help="X maximum", default='NULL')
        parser.add_option("--ymin", help="Y minimum", default='NULL')
        parser.add_option("--ymax", help="Y maximum", default='NULL')
        parser.add_option("--ymin_covplot", help="Y minimum (coverage plots)", default='NULL')
        parser.add_option("--ymax_covplot", help="Y maximum (coverage plots)", default='NULL')
        parser.add_option("--make_norm_plots", help="Make plots with normalized baselines", default=True)
        parser.add_option("--invert_strand_sense", help="Invert strand sense in plots", default=False, action='store_true')
        parser.add_option("--save_plot_csv", help="Save CSV files for plots", default=False, action='store_true')
        parser.add_option("--find_rois", help="Find regions of interest and save to GFF file", default=True, action='store_true')
        parser.add_option("--min_coverage", help="Minimum coverage to admit data", type='int', default=30)
        parser.add_option("--tempdir", default=tempfile.gettempdir())
        parser.add_option("--debug", action='store_true', default=True)
        self.opts, self.args = parser.parse_args()
        
        if self.opts.debug: logging.basicConfig(level=logging.DEBUG)
        
        if len(self.args) != 1:
            parser.error('Expected 1 arguments')

        self.csvbz_input_filename = self.args[0]

    def _generateParameters(self, report_type, input_csv_fn):
        self.params = {}
        self.params['inputData'] = input_csv_fn
        self.params['outputFile'] = str(self.opts.metric)+"_"+report_type+".pdf"
        self.params['outputFile'] = os.path.join(self.opts.output, self.params['outputFile'])
        self.params['ylabel'] = str(self.opts.metric)+" "+report_type
        self.params['title'] = str(self.opts.metric)+" "+report_type+": "+str(self.opts.title)
        self.params['ref_name'] = self.opts.ref_name
        self.params['xmin'] = self.opts.xmin
        self.params['xmax'] = self.opts.xmax
        if re.search("Coverage$", report_type):
            self.params['ymin'] = self.opts.ymin_covplot
            self.params['ymax'] = self.opts.ymax_covplot
        else:
            self.params['ymin'] = self.opts.ymin
            self.params['ymax'] = self.opts.ymax
        self.params['baseline'] = "0.0"
        if re.search("Ratio$", report_type):
            self.params['baseline'] = "1.0"
        self.params['top_strand_label'] = "Forward strand template>"
        self.params['bottom_strand_label'] = "<Reverse strand template"
        if self.opts.invert_strand_sense:
            self.params['top_strand_label'], self.params['bottom_strand_label'] = self.params['bottom_strand_label'], self.params['top_strand_label']

        return self.params
    
    def _generateScript(self, report_type, input_csv_fn):
        script = RScript(self.opts.script)
        script.construct(self._generateParameters(report_type, input_csv_fn))
        return script

    def _findROIs(self):
        ''' Use a simple rule to find regions of interest and output them into a gff file
        This will eventually need refactoring into the dmtk.model.kinetics framework;
        the logic here will reduce to "aggregate overlapping regions with significant p-values" (i.e. AggregatingFeatureWriter)
        '''
        poi = {'+': set(), '-': set()}
        gff_outfile = os.path.join(self.opts.output, str(self.opts.metric)+"_ROIs.gff")

        events_gff = GffIO.GffWriter(open(gff_outfile, 'w'))
        agg_writer = PulseKineticsUtils.AggregatingFeatureWriter(events_gff)

        # FIXME: HACK - see notes above
        min_coverage = 10
        min_coverage_extent = 8 # to either side
        min_ratio = 4

        min_coverage = 10
        min_coverage_extent = 4 # to either side
        min_ratio = 5
        
        coverage_queue = {'+': deque([], maxlen=min_coverage_extent*2 + 1),
                          '-': deque([], maxlen=min_coverage_extent*2 + 1)}
        ratio_queue = {'+': deque([], maxlen=min_coverage_extent + 1),
                       '-': deque([], maxlen=min_coverage_extent + 1)}
        num_events = 0
        
        with bz2.BZ2File(self.csvbz_input_filename, 'r') as input_csv:
            csv_reader = csv.reader(input_csv)
            header = csv_reader.next()
            prev_ref_pos = -1
            for line in csv_reader:
                ref_pos_data = dict(zip(header, line))
                ref_pos = int(ref_pos_data['RefPosn'])
                if ref_pos != prev_ref_pos + 1: # require continuous record in window - assumes ordered gff
                    for strand in '+', '-':
                        coverage_queue[strand].clear()
                        ratio_queue[strand].clear()
                    prev_ref_pos = ref_pos
                    continue
                prev_ref_pos = ref_pos
                
                coverage_queue['+'].append(min(int(ref_pos_data['ModifiedCoverageFwd']),
                                               int(ref_pos_data['ControlCoverageFwd'])))
                coverage_queue['-'].append(min(int(ref_pos_data['ModifiedCoverageRev']),
                                               int(ref_pos_data['ControlCoverageRev'])))
#                ratio_queue['+'].append(float64(ref_pos_data['ModifiedMeanFwd'])/float64(ref_pos_data['ControlMeanFwd']))
#                ratio_queue['-'].append(float64(ref_pos_data['ModifiedMeanRev'])/float64(ref_pos_data['ControlMeanRev']))
                ratio_queue['+'].append(float64(ref_pos_data['ModifiedMedianFwd'])/float64(ref_pos_data['ControlMedianFwd']))
                ratio_queue['-'].append(float64(ref_pos_data['ModifiedMedianRev'])/float64(ref_pos_data['ControlMedianRev']))
                
                for strand in '+', '-':
                    if len(coverage_queue[strand]) < coverage_queue[strand].maxlen: continue
                    if ratio_queue[strand][0] > min_ratio and min(coverage_queue[strand]) > min_coverage:
                        my_poi = ref_pos - min_coverage_extent
                        poi[strand].add(my_poi)
                        
                        r = GffIO.Gff3Record(seqName = self.opts.ref_name)
                        r.start = my_poi+1; r.end = my_poi+1; r.strand = strand; r.type = '.'
                        r.put('score', "%.4f" % ratio_queue[strand][0])
                        r.put('read_cov', "%d" % min(coverage_queue[strand]))
                        r.put('subread_cov', "%d" % min(coverage_queue[strand]))

                        agg_writer.addFeature(r)
                        num_events += 1
        del agg_writer # write out all the queues
        events_gff.close()
        logging.debug("Wrote %d events to %s" % (num_events, gff_outfile))

    def _findBaselineBias(self):
        mean_metrics = {'Modified': {'+': 1, '-': 1}, 'Control': {'+': 1, '-': 1}}
        median_metrics = {'Modified': {'+': 1, '-': 1}, 'Control': {'+': 1, '-': 1}}
        if self.opts.make_norm_plots:
            mean_metric_values = {'Modified': {'+': {}, '-': {}}, 'Control': {'+': {}, '-': {}}}
            median_metric_values = {'Modified': {'+': {}, '-': {}}, 'Control': {'+': {}, '-': {}}}
            poi = {'+': set(), '-': set()}
            with bz2.BZ2File(self.csvbz_input_filename, 'r') as input_csv:
                csv_reader = csv.reader(input_csv)
                header = csv_reader.next()
                for line in csv_reader:
                    ref_pos_data = dict(zip(header, line))
                    ref_pos = int(ref_pos_data['RefPosn'])
                    for condition in 'Modified', 'Control':
                        for strand in '+', '-':
                            mean_metric_values[condition][strand][ref_pos] = float(ref_pos_data[condition+'Mean'+('Fwd' if strand == '+' else 'Rev')])
                            median_metric_values[condition][strand][ref_pos] = float(ref_pos_data[condition+'Median'+('Fwd' if strand == '+' else 'Rev')])
                    if ref_pos_data['OfInterestFwd'] == '1':
                        for p in range(ref_pos-2, ref_pos+9):
                            poi['+'].add(p)
                    if ref_pos_data['OfInterestRev'] == '1':
                        for p in range(ref_pos-9, ref_pos+2):
                            poi['-'].add(p)
            
            for condition in 'Modified', 'Control':
                for strand in '+', '-':
                    my_means = numpy.array([m for p, m in mean_metric_values[condition][strand].iteritems() if p not in poi[strand]])
                    nanmask = numpy.isfinite(my_means)
                    if sum(nanmask) > 0: my_means = my_means[nanmask]
                    if len(my_means) == 0: continue
                    mean_metrics[condition][strand] = numpy.mean(my_means)
            for condition in 'Modified', 'Control':
                for strand in '+', '-':
                    my_medians = numpy.array([m for p, m in median_metric_values[condition][strand].iteritems() if p not in poi[strand]])
                    nanmask = numpy.isfinite(my_medians)
                    if sum(nanmask) > 0: my_medians = my_medians[nanmask]
                    if len(my_medians) == 0: continue
                    median_metrics[condition][strand] = numpy.mean(my_medians)
        
        return mean_metrics, median_metrics

    def _getWindowPvalues(self):
        gff_outfile = os.path.join(self.opts.output, str(self.opts.metric)+"_ROIs_by_pvalue.gff")
        events_gff = GffIO.GffWriter(open(gff_outfile, 'w'))
        agg_writer = PulseKineticsUtils.AggregatingFeatureWriter(events_gff)
        
        window_pvalues = {'+': {}, '-': {}}
        bg_scores = {}
        
        min_coverage = 10
        min_coverage_extent = 3 # to either side - TODO: asymmetric queue to reflect footprint
        #min_score  # to report to gff - FIXME: requires shuffle calibration, coverage calibration
        
        coverage_queue = {'+': deque([], maxlen=min_coverage_extent*2 + 1),
                          '-': deque([], maxlen=min_coverage_extent*2 + 1)}
        pvalue_queue = {'+': deque([], maxlen=min_coverage_extent*2 + 1),
                       '-': deque([], maxlen=min_coverage_extent*2 + 1)}
        num_events = 0
        
        
        with bz2.BZ2File(self.csvbz_input_filename, 'r') as input_csv:
            csv_reader = csv.reader(input_csv)
            header = csv_reader.next()
            prev_ref_pos = -1
            for line in csv_reader:
                ref_pos_data = dict(zip(header, line))
                ref_pos = int(ref_pos_data['RefPosn'])
                if ref_pos != prev_ref_pos + 1: # require continuous record in window - assumes ordered gff
                    for strand in '+', '-':
                        coverage_queue[strand].clear()
                        pvalue_queue[strand].clear()
                    prev_ref_pos = ref_pos
                    continue
                prev_ref_pos = ref_pos
                
                coverage_queue['+'].append(min(int(ref_pos_data['ModifiedCoverageFwd']),
                                               int(ref_pos_data['ControlCoverageFwd'])))
                coverage_queue['-'].append(min(int(ref_pos_data['ModifiedCoverageRev']),
                                               int(ref_pos_data['ControlCoverageRev'])))
                pvalue_queue['+'].append(float64(ref_pos_data['PvalueFwd']))
                pvalue_queue['-'].append(float64(ref_pos_data['PvalueRev']))
                
                for strand in '+', '-':
                    if len(coverage_queue[strand]) < coverage_queue[strand].maxlen: continue
                    if min(coverage_queue[strand]) > min_coverage:
                        my_score = -log10(product(pvalue_queue[strand]) + 1e-300)
                        bg_scores.setdefault(min(coverage_queue[strand]), [])
                        bg_scores[min(coverage_queue[strand])].append(my_score)
        for strand in '+', '-':
            coverage_queue[strand].clear()
            pvalue_queue[strand].clear()
        
        for cov in bg_scores.keys():
            bg_scores[cov] = numpy.median(bg_scores[cov])
        
        
        with bz2.BZ2File(self.csvbz_input_filename, 'r') as input_csv:
            csv_reader = csv.reader(input_csv)
            header = csv_reader.next()
            prev_ref_pos = -1
            for line in csv_reader:
                ref_pos_data = dict(zip(header, line))
                ref_pos = int(ref_pos_data['RefPosn'])
                if ref_pos != prev_ref_pos + 1: # require continuous record in window - assumes ordered gff
                    for strand in '+', '-':
                        coverage_queue[strand].clear()
                        pvalue_queue[strand].clear()
                    prev_ref_pos = ref_pos
                    continue
                prev_ref_pos = ref_pos
                
                coverage_queue['+'].append(min(int(ref_pos_data['ModifiedCoverageFwd']),
                                               int(ref_pos_data['ControlCoverageFwd'])))
                coverage_queue['-'].append(min(int(ref_pos_data['ModifiedCoverageRev']),
                                               int(ref_pos_data['ControlCoverageRev'])))
                pvalue_queue['+'].append(float64(ref_pos_data['PvalueFwd']))
                pvalue_queue['-'].append(float64(ref_pos_data['PvalueRev']))
                
                for strand in '+', '-':
                    if len(coverage_queue[strand]) < coverage_queue[strand].maxlen: continue
                    if min(coverage_queue[strand]) > min_coverage:
                        my_pvalue = product(pvalue_queue[strand]) + 1e-300
                        my_score = -log10(my_pvalue)
                        my_norm_score = my_score / bg_scores[min(coverage_queue[strand])]
                        window_pvalues[strand][ref_pos - min_coverage_extent] = my_pvalue
                        
                        if my_norm_score > 2: # TODO: this needs pvalue calibration by shuffling
                            my_poi = ref_pos - min_coverage_extent
                            
                            
                            r = GffIO.Gff3Record(seqName = self.opts.ref_name)
                            r.start = my_poi+1; r.end = my_poi+1; r.strand = strand; r.type = '.'
                            r.put('score', "%.4f" % my_score)
                            r.put('norm_score', "%.4f" % my_norm_score)
                            r.put('pvalue', "%f" % my_pvalue)
                            r.put('read_cov', "%d" % min(coverage_queue[strand]))
                            r.put('subread_cov', "%d" % min(coverage_queue[strand]))
    
                            agg_writer.addFeature(r)
                            num_events += 1

        del agg_writer # write out all the queues
        events_gff.close()
        logging.debug("Wrote %d events to %s" % (num_events, gff_outfile))
        
        return window_pvalues

    def run(self):
        mean_metrics, median_metrics = self._findBaselineBias()
        window_pvalues = self._getWindowPvalues()
        if self.opts.find_rois:
            self._findROIs()
        
        report_filenames = {}; report_fhs = {}; report_csv_writers = {}
        bed_report_filenames = {}; bed_report_fhs = {}; bed_writers = {}
        
        header = ['RefPosn', 'myVarFwdStrand', 'spreadLoFwdStrand', 'spreadHiFwdStrand',
                  'myVarRevStrand', 'spreadLoRevStrand', 'spreadHiRevStrand',
                  'ModifiedFwdStrand', 'ModifiedRevStrand',
                  'nucleotide', 'nucleotideRevStrand', 'refNucleotide', 'refNucleotideRevStrand',
                  'effectiveCoverageFwdStrand', 'effectiveCoverageRevStrand']
        
        #report_types = 'ratio', 'diff', 'llr', 'norm_ratio', 'log_norm_ratio', 'norm_diff', 'norm_llr', 'm1', 'm2', 'od1', 'od2', 'od_ratio':
        report_types = ('Ratio', 'NormalizedRatio', 'NormalizedMedianRatio', 'diff', 'LLR', 'nlog10Pvalue', 'ModifiedMean', 'ControlMean', 'EffectiveModifiedCoverage', 'EffectiveControlCoverage', 'EffectiveJointCoverage')
        bed_report_types = ('NormalizedRatio', 'NormalizedMedianRatio', 'nlog10Pvalue')

        for report_type in report_types:
            report_filenames[report_type] = tempfile.mkstemp(prefix="makePulseKineticsReport_"+str(self.opts.metric)+"_"+report_type, suffix='.csv')[1]
            report_fhs[report_type] = open(report_filenames[report_type], 'w')
            report_csv_writers[report_type] = csv.writer(report_fhs[report_type])
            report_csv_writers[report_type].writerow(header)

        for report_type in bed_report_types:
            bed_report_filenames[report_type] = os.path.join(self.opts.output, str(self.opts.metric)+"_"+report_type+".bed.gz")
            bed_report_fhs[report_type] = gzip.open(bed_report_filenames[report_type], 'w')
            bed_writers[report_type] = BedWriter(bed_report_fhs[report_type])
            track_description = str(self.opts.metric)+" "+report_type+": "+str(self.opts.title)
            bed_writers[report_type].writeHeader(report_type, track_description, 1)

        def _preprocess_row(row):
            for i in range(len(row)):
                if row[i] == None or str(row[i]) == 'nan' or str(row[i]) == 'inf':
                    row[i] = 'NA'
                else:
                    row[i] = str(row[i])
            return row
        
        def _cast_row(row):
            for i in range(len(row)):
                try:
                    row[i] = int(row[i])
                except ValueError:
                    try:
                        row[i] = float64(row[i])
                    except ValueError:
                        pass
            return row
        
        with bz2.BZ2File(self.csvbz_input_filename, 'r') as input_csv:
            csv_reader = csv.reader(input_csv)
            header = csv_reader.next()
            for line in csv_reader:
                line = _cast_row(line)
                ref_pos_data = dict(zip(header, line))
                fwd_diff = ref_pos_data['ModifiedMeanFwd'] - ref_pos_data['ControlMeanFwd']
                rev_diff = ref_pos_data['ModifiedMeanRev'] - ref_pos_data['ControlMeanRev']

                fwd_ratio = ref_pos_data['ModifiedMeanFwd'] / ref_pos_data['ControlMeanFwd']
                rev_ratio = ref_pos_data['ModifiedMeanRev'] / ref_pos_data['ControlMeanRev']

                fwd_se1 = (ref_pos_data['ModifiedSpreadHiFwd'] - ref_pos_data['ModifiedMeanFwd']) * sqrt(ref_pos_data['ModifiedCoverageFwd'])
                fwd_se2 = (ref_pos_data['ControlSpreadHiFwd'] - ref_pos_data['ControlMeanFwd']) * sqrt(ref_pos_data['ControlCoverageFwd'])
                rev_se1 = (ref_pos_data['ModifiedSpreadHiRev'] - ref_pos_data['ModifiedMeanRev']) * sqrt(ref_pos_data['ModifiedCoverageRev'])
                rev_se2 = (ref_pos_data['ControlSpreadHiRev'] - ref_pos_data['ControlMeanRev']) * sqrt(ref_pos_data['ControlCoverageRev'])

                fwd_tvar = fwd_se2**2 * ref_pos_data['ModifiedMeanFwd']**2 / ref_pos_data['ControlMeanFwd']**4 + fwd_se1**2 / ref_pos_data['ControlMeanFwd']**2
                fwd_tvar /= min(ref_pos_data['ModifiedCoverageFwd'], ref_pos_data['ControlCoverageFwd'])
                fwd_ratio_se = sqrt(fwd_tvar)

                rev_tvar = rev_se2**2 * ref_pos_data['ModifiedMeanRev']**2 / ref_pos_data['ControlMeanRev']**4 + rev_se1**2 / ref_pos_data['ControlMeanRev']**2
                rev_tvar /= min(ref_pos_data['ModifiedCoverageRev'], ref_pos_data['ControlCoverageRev'])
                rev_ratio_se = sqrt(rev_tvar)
                
#                 t_value = 12.71 # 2-sided, 1-dof 95% CI
#                 t_value = 1.0 # disable SEM -> CI conversion
#                 fwd_se1 = ref_pos_data['ModifiedSpreadHiFwd'] - ref_pos_data['ModifiedMeanFwd']
#                 fwd_se2 = ref_pos_data['ControlSpreadHiFwd'] - ref_pos_data['ControlMeanFwd']
#                 fwd_norm_se1 = fwd_se1/ref_pos_data['ModifiedMeanFwd']
#                 fwd_norm_se2 = fwd_se2/ref_pos_data['ControlMeanFwd']                
#                 fwd_ratio_se = fwd_ratio * sqrt(fwd_norm_se1**2 + fwd_norm_se2**2) * t_value
                
#                 rev_se1 = ref_pos_data['ModifiedSpreadHiRev'] - ref_pos_data['ModifiedMeanRev']
#                 rev_se2 = ref_pos_data['ControlSpreadHiRev'] - ref_pos_data['ControlMeanRev']
#                 rev_norm_se1 = rev_se1/ref_pos_data['ModifiedMeanRev']
#                 rev_norm_se2 = rev_se2/ref_pos_data['ControlMeanRev']
#                 rev_ratio_se = rev_ratio * sqrt(rev_norm_se1**2 + rev_norm_se2**2) * t_value

                norm_fwd_ratio = (ref_pos_data['ModifiedMeanFwd']/mean_metrics['Modified']['+']) / (ref_pos_data['ControlMeanFwd']/mean_metrics['Control']['+'])
                norm_rev_ratio = (ref_pos_data['ModifiedMeanRev']/mean_metrics['Modified']['-']) / (ref_pos_data['ControlMeanRev']/mean_metrics['Control']['-'])

                norm_fwd_median_ratio = (ref_pos_data['ModifiedMedianFwd']/median_metrics['Modified']['+']) / (ref_pos_data['ControlMedianFwd']/median_metrics['Control']['+'])
                norm_rev_median_ratio = (ref_pos_data['ModifiedMedianRev']/median_metrics['Modified']['-']) / (ref_pos_data['ControlMedianRev']/median_metrics['Control']['-'])
                
                try:
                    pvalue_fwd = window_pvalues['+'][ref_pos_data['RefPosn']]
                except KeyError:
                    pvalue_fwd = 1.0
                try:
                    pvalue_rev = window_pvalues['-'][ref_pos_data['RefPosn']]
                except KeyError:
                    pvalue_rev = 1.0

#                if min(ref_pos_data['ModifiedCoverageFwd'], ref_pos_data['ControlCoverageFwd']) < self.opts.min_coverage:
#                    fwd_ratio, norm_fwd_ratio, norm_fwd_median_ratio, fwd_diff, ref_pos_data['LLRFwd']= None, None, None, None, None
#                if min(ref_pos_data['ModifiedCoverageRev'], ref_pos_data['ControlCoverageRev']) < self.opts.min_coverage:
#                    rev_ratio, norm_rev_ratio, norm_rev_median_ratio, rev_diff, ref_pos_data['LLRRev'] = None, None, None, None, None

                my_report_lines = {}
                for report_type in report_types:
                    my_report_lines[report_type] = [ref_pos_data['RefPosn']]

                my_report_lines['Ratio'].extend([fwd_ratio, fwd_ratio-fwd_ratio_se, fwd_ratio+fwd_ratio_se, rev_ratio, rev_ratio-rev_ratio_se, rev_ratio+rev_ratio_se])
                my_report_lines['NormalizedRatio'].extend([norm_fwd_ratio, norm_fwd_ratio-fwd_ratio_se, norm_fwd_ratio+fwd_ratio_se, norm_rev_ratio, norm_rev_ratio-rev_ratio_se, norm_rev_ratio+rev_ratio_se])
                my_report_lines['NormalizedMedianRatio'].extend([norm_fwd_median_ratio, norm_fwd_median_ratio, norm_fwd_median_ratio, norm_rev_median_ratio, norm_rev_median_ratio, norm_rev_median_ratio])
                my_report_lines['diff'].extend([fwd_diff, fwd_diff, fwd_diff, rev_diff, rev_diff, rev_diff])
                my_report_lines['LLR'].extend([ref_pos_data['LLRFwd'], ref_pos_data['LLRFwd'], ref_pos_data['LLRFwd'],
                                               ref_pos_data['LLRRev'], ref_pos_data['LLRRev'], ref_pos_data['LLRRev']])
#                my_report_lines['nlog10Pvalue'].extend([-log10(ref_pos_data['PvalueFwd']), -log10(ref_pos_data['PvalueFwd']), -log10(ref_pos_data['PvalueFwd']),
#                                               -log10(ref_pos_data['PvalueRev']), -log10(ref_pos_data['PvalueRev']), -log10(ref_pos_data['PvalueRev'])])
                my_report_lines['nlog10Pvalue'].extend([-log10(pvalue_fwd), -log10(pvalue_fwd), -log10(pvalue_fwd),
                                               -log10(pvalue_rev), -log10(pvalue_rev), -log10(pvalue_rev)])
                my_report_lines['ModifiedMean'].extend([ref_pos_data['ModifiedMeanFwd'], ref_pos_data['ModifiedSpreadLoFwd'], ref_pos_data['ModifiedSpreadHiFwd'],
                                              ref_pos_data['ModifiedMeanRev'], ref_pos_data['ModifiedSpreadLoRev'], ref_pos_data['ModifiedSpreadHiRev']])
                my_report_lines['ControlMean'].extend([ref_pos_data['ControlMeanFwd'], ref_pos_data['ControlSpreadLoFwd'], ref_pos_data['ControlSpreadHiFwd'],
                                              ref_pos_data['ControlMeanRev'], ref_pos_data['ControlSpreadLoRev'], ref_pos_data['ControlSpreadHiRev']])
                my_report_lines['EffectiveModifiedCoverage'].extend([ref_pos_data['ModifiedCoverageFwd'], None, None,
                                                     ref_pos_data['ModifiedCoverageRev'], None, None])
                my_report_lines['EffectiveControlCoverage'].extend([ref_pos_data['ControlCoverageFwd'], None, None,
                                                     ref_pos_data['ControlCoverageRev'], None, None])
                my_report_lines['EffectiveJointCoverage'].extend([min(ref_pos_data['ModifiedCoverageFwd'], ref_pos_data['ControlCoverageFwd']), None, None,
                                                          min(ref_pos_data['ModifiedCoverageRev'], ref_pos_data['ControlCoverageRev']), None, None])
                
                for report_type in report_types:
                    my_report_lines[report_type].extend([ref_pos_data['OfInterestFwd'], ref_pos_data['OfInterestRev'],
                                                         ref_pos_data['NucleotideFwd'], ref_pos_data['NucleotideRev'],
                                                         ref_pos_data['RefNucleotideFwd'], ref_pos_data['RefNucleotideRev'],
                                                         min(ref_pos_data['ModifiedCoverageFwd'], ref_pos_data['ControlCoverageFwd']),
                                                         min(ref_pos_data['ModifiedCoverageRev'], ref_pos_data['ControlCoverageRev'])])
                    
                    if self.opts.invert_strand_sense:
                        line = my_report_lines[report_type][:]
                        my_report_lines[report_type] = [line[0]] + line[4:7] + line[1:4] + [line[8], line[7], line[10], line[9], line[12], line[11], line[14], line[13]]
                    
                    report_csv_writers[report_type].writerow(_preprocess_row(my_report_lines[report_type]))

                common_bed_kwargs = dict(chrom = self.opts.ref_name, chromStart = ref_pos_data['RefPosn'], chromEnd = ref_pos_data['RefPosn'])
                if min(ref_pos_data['ModifiedCoverageFwd'], ref_pos_data['ControlCoverageFwd']) > self.opts.min_coverage:
                    common_bed_kwargs['strand'] = '+'
                    bed_writers['NormalizedRatio'].writeRecord(BedRecord(name = 'NormalizedRatio', score = norm_fwd_ratio, **common_bed_kwargs))
                    bed_writers['NormalizedMedianRatio'].writeRecord(BedRecord(name = 'NormalizedMedianRatio', score = norm_fwd_median_ratio, **common_bed_kwargs))
                    bed_writers['nlog10Pvalue'].writeRecord(BedRecord(name = 'nlog10Pvalue', score = -log10(pvalue_fwd), **common_bed_kwargs))

                if min(ref_pos_data['ModifiedCoverageRev'], ref_pos_data['ControlCoverageRev']) > self.opts.min_coverage:
                    common_bed_kwargs['strand'] = '-'
                    bed_writers['NormalizedRatio'].writeRecord(BedRecord(name = 'NormalizedRatio', score = norm_rev_ratio, **common_bed_kwargs))
                    bed_writers['NormalizedMedianRatio'].writeRecord(BedRecord(name = 'NormalizedMedianRatio', score = norm_rev_median_ratio, **common_bed_kwargs))
                    bed_writers['nlog10Pvalue'].writeRecord(BedRecord(name = 'nlog10Pvalue', score = -log10(pvalue_rev), **common_bed_kwargs))
        
        child_process_handles = []
        for report_type in report_types:
            report_fhs[report_type].close()
            
            if self.opts.save_plot_csv:
                shutil.copy(report_filenames[report_type], self.opts.output)
        
        for report_type in bed_report_types:
            bed_report_fhs[report_type].close()

        for report_type in report_types:
            rScript = self._generateScript(report_type, report_filenames[report_type])
            scriptFile = rScript.save()
            rCmd = '%s CMD BATCH %s /dev/null' % (R_BIN, scriptFile)
            
            logging.debug("Running '%s' (%s %s)" % (rCmd, self.opts.metric, report_type))
            ph = subprocess.Popen(rCmd, shell=True)
            child_process_handles.append(ph)
            
#            try:
#                subprocess.check_call(rCmd, shell=True)
#            except subprocess.CalledProcessError:
#                sys.stderr.write("R Command Failed: (%s), Unable to generate Pulse Kinetics Report.\n" % rCmd)
#                return

        for ph in child_process_handles:
            if ph: ph.wait()

        
        '''
        report = GraphReportItem()
        report.layout = 'onecolumn'
        report.title = 'Pulse Kinetics Visualizations'

        graphItem = GraphItem()
        graphItem.title = ''       
        graphItem.addImage(os.path.basename(self.params['accuracyComparison']))
        report.addGraph(graphItem)
        
        graphItem = GraphItem()
        graphItem.title = ''       
        graphItem.addImage(os.path.basename(self.params['etaFit']))
        report.addGraph(graphItem)
        
        report.saveToFile(sys.stdout)
        '''

if __name__=='__main__':
    app = MakePulseKineticsReport()
    app.run()
