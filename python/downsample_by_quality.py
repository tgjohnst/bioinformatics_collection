# -*- coding: utf-8 -*-
"""
@author: TJ
Downsamples a fastq file to a certain number of records per barcode.
Assumes the fastq has been preprocessed (e.g. PRESTO) to include barcode field in header.
By default, selects those n reads with the highest mean base quality
"""

import os
import sys
import argparse
import math
import heapq
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from presto.Annotation import parseAnnotation
import numpy as np


def parse_cmd(args):
    """ Parses user-supplied command line arguments """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=str, required=True)
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('-n', '--num_per_group', type=int, default=50)
    parser.add_argument('-b', '--barcode_field', type=str, required=True)
    parser.add_argument('-r', '--random_selection', action='store_true')
    optS = parser.parse_args(args)
    return optS


def mean_base_quality(record):
    ''' gets the mean phred base quality for a SeqRecord '''
    return np.mean(record.letter_annotations['phred_quality'])


class ComparableRecord:
    ''' 
    Stores a SeqRecord ID string but allows comparison based on mean base quality 
    This allows us to maintain fixed sized heaps, since heapq does not support
    custom comparators for functions other than nsmallest/nlargest and SeqRecord
    do not have built in comparison
    '''
    def __init__(self, key, mbq):
        self.key = key
        self.mbq = mbq
    
    def __lt__(self, other):
        return self.mbq < other.mbq


def shuffle_fastq(record_dict):
    '''
    Shuffles a fastq into random order and return a list of the shuffled keys
    '''
    # shuffle the keys of the dict with a set random.seed for determinism
    shuffled_keys = list(record_dict.keys())
    shuffled_keys.sort()
    random.seed('Cav is fast')
    random.shuffle(shuffled_keys)
    return shuffled_keys


def write_read_with_count(record, count, handle):
    record.description = record.description + '|ORIG_COUNT=' + str(count)
    handle.write(record.format('fastq'))


def get_totals_per_barcode(record_dict, barcode):
    '''
    Returns a dictionary with the number of reads sharing each barcode
    '''
    print('Counting barcodes...', file=sys.stderr)
    counts = defaultdict(int)
    for key in record_dict.keys():
        annot = parseAnnotation(record_dict[key].description)
        counts[annot[barcode]] += 1
    return counts


def first_n_per_group(record_dict, shuffled_keys, outfile, totals, barcode, num):
    '''
    Takes the first n reads per group from infile and writes to outfile
    '''
    print('Randomly selecting ' + str(num) + ' reads per barcode', file=sys.stderr)
    with open(str(outfile), 'w') as outh:
        seen = defaultdict(int)
        for key in shuffled_keys:
            annot = parseAnnotation(record_dict[key].description)
            seen[annot[barcode]] += 1
            if seen[annot[barcode]] <= num:
                write_read_with_count(record_dict[key], totals[annot[barcode]], outh)


def top_n_per_group(record_dict, outfile, totals, barcode, num):
    '''
    Takes the top n reads per group from infile based on the mean base quality.
    Write them out to outfile.
    '''
    print('Selecting top ' + str(num) + ' reads per barcode based on quality...', file=sys.stderr)
    top_reads = defaultdict(list)
    with open(str(outfile), 'w') as outh:
        best_quals = defaultdict(list)
        for key in record_dict.keys():
            annot = parseAnnotation(record_dict[key].description)
            bc = annot[barcode]
            # don't write out reads to file yet if conscount > subsample
            if totals[bc] > num:
                cr = ComparableRecord(key, mean_base_quality(record_dict[key]))
                # if we havent' seen num reads yet, just append
                if len(top_reads[bc]) < num:
                    top_reads[bc].append(cr)
                # if we have, keep fixed size heap and pushpop!
                else:
                    heapq.heappushpop(top_reads[bc], cr)
            else:
                write_read_with_count(record_dict[key], totals[bc], outh)
        # Finally, write out our collected top reads per barcode to the outfile
        for bc, reads in top_reads.items():
            for cr in reads:
                write_read_with_count(record_dict[cr.key], totals[bc], outh)


def main(argv=None):
    args = parse_cmd(argv[1:])
    # index the file without reading into memory
    record_dict = SeqIO.index(str(args.infile), 'fastq')
    totals = get_totals_per_barcode(record_dict, args.barcode_field)
    # Randomly select n reads per barcode
    if args.random_selection:
        shuffled = shuffle_fastq(record_dict)
        first_n_per_group(record_dict, shuffled, args.outfile, totals,
                          args.barcode_field, args.num_per_group)
    # default behavior - select top n reads by quality
    else:
        top_n_per_group(record_dict, args.outfile, totals,
                        args.barcode_field, args.num_per_group)


if __name__ == '__main__':
    main(sys.argv)

