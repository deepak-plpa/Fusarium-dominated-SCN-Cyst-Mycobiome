#!/usr/bin/env python
#
# Trim all ab1 files and merge them all into a single file
#
# Taejoon Byun <taejoon@umn.edu>
# last : Feb 28 2017

# @Deepak:  if the python interpreter in your machine complains about 
#           missing module(s), please run the following command:
#           `pip --user install <module>`
import abifpy as abi
import sys, copy, glob, os

class QualityRule(object):
    ''' An abstract class that defines the interface among quality rules.''' 

    trimmed = None      # the trimmed trace to check
    original = None     # the original trace to check against

    def __init__(self, original, trimmed):
        self.original = original
        self.trimmed = trimmed

    def is_satisfied(self):
        ''' :returns: True if the corresponding rule is is_satisfied. ''' 
        raise NotImplementedError


class RuleNoConsecutiveNs(QualityRule): 
    ''' Rule 1: No more than ``count`` consecutive Ns ''' 

    count = 0   # number of consecutive Ns to allow

    def __init__(self, original, trimmed, count=2):
        QualityRule.__init__(self, original, trimmed)
        self.count = count

    def is_satisfied(self):
        cnt_n = 0
        for val in self.trimmed.seq:
            if val is 'N':
                cnt_n += 1
            elif cnt_n is not 0:
                # not 'N'
                cnt_n = 0
            if cnt_n >= 2:
                print '\t[error] More than %d consecutive Ns at: %s' \
                        % (self.count, self.trimmed.name)
                return False
        return True

class RuleTraceLength(QualityRule): 
    ''' Rule 2: Trimmed sequence shall be as long as 70% of the length of the 
        original trace. '''
    def __init__(self, original, trimmed):
        QualityRule.__init__(self, original, trimmed)

    def is_satisfied(self):
        if len(self.trimmed.seq) < 0.7 * len(self.original.seq):
            # length is lesser than 70% of the original, not useful
            print '\t[error] Too short: ' + self.trimmed.name
            return False
        return True

class RuleAverageQuality(QualityRule): 
    ''' Rule 3: Average qual value shall be greater than or equal to 30. '''
    def __init__(self, original, trimmed):
        QualityRule.__init__(self, original, trimmed)

    def _get_average(self, lst):
        ''' :returns: the average value of all the int/floats in the ``lst``.'''
        return reduce(lambda x, y: float(x + y), lst) / len(lst)

    def is_satisfied(self):
        if self._get_average(self.trimmed.qual_val) < 30:
            # intensity too low
            print '\t[error] Intensity too low: ' + self.trimmed.name
            return False
        return True
        

class Ab1Trimmer(object):

    _rules = []

    def __init__(self, rules):
        ''' Initialize the trimmer.
        :param rules: the name of the rules (subclasses of ``QualityRule``) to 
        check'''
        self._rules = rules

    def _get_average(self, lst):
        ''' :returns: the average value of all the int/floats in the ``lst``.'''
        return reduce(lambda x, y: float(x + y), lst) / len(lst)

    def _is_qualified(self, original, trimmed):
        ''' Cross check the original and trimmed trace against all the rules.
        :param original:    the original trace
        :param trim:        the trimmed trace
        :returns:           True if qualified, False if not qualified.
        '''
        rules = [Rule(original, trimmed) for Rule in self._rules]
        for rule in rules:
            if not rule.is_satisfied():
                return False
        return True

    def _get_head_ind_to_trim(self, seq, qual_vals):
        ''' A simple heuristic for trimming a sequence.
        :returns: the index of the sequence to start trimming from.'''
        # @Deepak:  Please modify this function to implement how you trim the 
        #           sequence.
        cnt_consecutive_good = 0
        head_ind = 0    # return value
        # trim the head
        assert len(qual_vals) == len(seq)
        for i in range(0, len(qual_vals) / 3):
            # only look at the first 1/3 of the whole sequence
            if qual_vals[i] < max(qual_vals) / 2:
                # quality value is lesser than a half of the max qual value
                head_ind = i + 1
                cnt_consecutive_good = 0
            elif seq[i] == 'N':
                # 'N' is not good
                head_ind = i + 1
        return head_ind

    def _trim(self, trace):
        ''' Trim a trace
        :param trace: a trace to trim
        :returns: (nullable) the trimmed trace if the trace satisfies all the 
        rules, None if it doesn't. '''
        print 'Trimming ' + trace.name
        trimmed = copy.copy(trace)
        maximum = max(trace.qual_val)
        average = self._get_average(trimmed.qual_val)
        print '\tlength: %d, max %d, avg: %d' \
                % (len(trace.seq), maximum, average)
        # trim
        qual_val = copy.deepcopy(trace.qual_val)
        seq = copy.deepcopy(trace.seq)
        start = self._get_head_ind_to_trim(seq, qual_val)
        qual_val.reverse()
        seq = seq[::-1]     # reverse
        end = self._get_head_ind_to_trim(seq, qual_val)
        end = len(trace.seq) - 1 - end
        trimmed.seq = trimmed.seq[start:end]
        print '\ttrimming from %d to %d' % (start, end)
        if not self._is_qualified(trace, trimmed):
            print '\t[error] Discarding an unqualified trace: ' + trimmed.name
            return None
        else:
            return trimmed

    def _save_fasta(self, trimmed, fname):
        print '\tsaving a trimmed trace to ' + fname
        trimmed.export(fname, fmt='fasta')

    def _trim_all_ab1s(self, ab1_files):
        print 'Trimming %d ab1 files' % len(ab1_files)
        original_list = [abi.Trace(fname) for fname in ab1_files]
        for trace in original_list:
            trimmed = self._trim(trace)
            if trimmed is not None:
                yield trimmed

    def trim_single_file(self, fname):
        ''' Trim a single file and save it to fasta format. 
        Do nothing if the trimmed trace fails any validation rule.'''
        original = abi.Trace(fname)
        trimmed = self._trim(original)
        if trimmed is not None:
            out_fname = fname.split('.')[-2] + '_trimmed.fas'
            self._save_fasta(trimmed, out_fname)

    def trim_and_merge(self, dirname):
        ''' Trim a bunch of ``ab1`` files in a directory and merge them into a
        single ``.fasta`` file.
        :param dirname: The name of a directory where ``.ab`` files are stored.
        '''
        def get_fas_fname(outdir, ab1fname):
            return outdir + ab1fname.split('/')[-1].split('.')[-2] + '.fas'

        ab1_files = glob.glob(dirname + '/*.ab1')
        outdir = dirname + '/trimmed'
        try:
            os.mkdir(outdir)
        except OSError as e:
            pass
        ## trim and save all qualified traces ##
        for trimmed in self._trim_all_ab1s(ab1_files):
            self._save_fasta(trimmed, '%s/%s.fas' % (outdir, trimmed.name))

        ## merge all ##
        merged_fname = outdir + '/merged.fas'
        fasta_files = glob.glob(outdir + '/*.fas')
        print 'Merging %d/%d files to %s' \
                % (len(fasta_files), len(ab1_files), merged_fname)
        merged = []
        for fasta in fasta_files:
            with open(fasta, 'r') as f:
                merged += f.readlines()
        with open(merged_fname, 'w') as f:
            f.writelines(merged)
        

def trim_single_file(rules):
    def fail():
        print 'Invalid argument'
        print 'Usage: $ %s <ab1_file>' % sys.argv[0]
        sys.exit()

    if len(sys.argv) != 2:
        fail()
    fname = sys.argv[1]
    if not os.path.isfile(fname):
        print 'Invalid file: %s' % fname
        fail()
    trimmer = Ab1Trimmer(rules)
    trimmer.trim_single_file(fname)

def trim_a_bunch(rules):
    def fail():
        print 'Invalid argument'
        print 'Usage: $ %s <ab1_directory>' % sys.argv[0]
        sys.exit()
    if len(sys.argv) != 2:
        fail()
    dirname = sys.argv[1]
    if not os.path.isdir(dirname):
        print 'Invalid directory: %s' % dirname
        fail()
    trimmer = Ab1Trimmer(rules)
    trimmer.trim_and_merge(dirname)


def main():
    # mode of operation
    SINGLE_FILE = False     # set this flag to True to run in single-file mode
    DIRECTORY = not SINGLE_FILE

    rules = [RuleNoConsecutiveNs, RuleTraceLength, RuleAverageQuality]
    if SINGLE_FILE:
        trim_single_file(rules)
        return
    if DIRECTORY:
        trim_a_bunch(rules)
        return


if __name__ == '__main__':
    main()

