'''
Functions to handle and process primers.
'''

import os
import sys
import re

def import_bed_file(bed):
    '''
    Import a BED file containing regions and return a matrix containing those
    regions.

    Arguments:
        * bed: full path to a BED file to process

    Return Value:
        Returns a list of primers from the BED file
    '''
    primers = []
    with open(bed, 'r') as bed_f:
        for line in bed_f:
            line = line.rstrip()
            tmp_data = line.split(sep='\t')
            primers.append(tmp_data)
    return primers


def create_primer_pairs(primers, left='_LEFT', right='_RIGHT'):
    '''
    From a list of primers, create a dictionary containing the left and right
    regions in a single entry.

    Arguments:
        * primers:  a list containing primers to process
        * left:     string pattern for the left primer in the id
        * right:    string pattern for the right primer in the id

    Return Value:
        Returns a dictionary with the primer ID as key and the following:
            * left_start:       start position of the left primer
            * left_end:         end position of the left primer
            * right_start:      start position of the right primer
            * right_end:        end position of the right primer
    '''
    primer_pairs = {}
    tmp_primer = {}
    pattern = ''.join(["(", left, '|', right, ")"])
    for primer in primers:
        primer_id = re.sub(pattern, '', primer[3])
        if primer_id not in primer_pairs.keys():
            primer_pairs[primer_id] = {}
        if re.search(left, primer[3]):
            primer_pairs[primer_id]['left_start'] = int(primer[1])
            primer_pairs[primer_id]['left_end'] = int(primer[2])
        elif re.search(right, primer[3]):
            primer_pairs[primer_id]['right_start'] = int(primer[1])
            primer_pairs[primer_id]['right_end'] = int(primer[2])
        else:
            print('skipping')
            print(primer)
    return primer_pairs


def create_amplicon_range(primer_pairs, pattern):
    '''
    Convert the primer pairs to a BED formatted amplicon range list.

    Arguments:
        * primer_pairs:     output from the create_primer_pairs function
        * pattern:          pattern in the amplicon key to get the amplicon
                            number

    Return Values:
        Returns a list of amplicon ranges in BED format
    '''
    amplicon_range = []
    for amplicon in primer_pairs:
        if all (k in primer_pairs[amplicon] for k in ('left_start', 'right_end')):
            amplicon_num = re.sub(pattern, '', amplicon)
            amplicon_range.append(create_bed_item(ref='MN908947.3',
                                                start=primer_pairs[amplicon]['left_end'],
                                                end=primer_pairs[amplicon]['right_start'],
                                                name=amplicon_num,
                                                score=amplicon))
    return amplicon_range


def create_bed_item(ref, start, end, name, score='100', strand='+'):
    '''
    Create a BED formatted list.

    Arguments:
        * ref: reference for genome/chromosome
        * start: start position for amplicon
        * end: end position for amplicon
        * name: name of amplicon
        * score: score of amplicon
        * strand: strand of the amplicon

    Return Values:
        Returns a list
    '''
    return [str(ref),
            str(start),
            str(end),
            str(name),
            str(score),
            str(strand)]


def create_unique_amplicons(amplicons, offset=30):
    '''
    Create a list of unique regions from the amplicons.

    Arguments:
        * amplicons:    a list of amplicons in BED format from
                        create_amplicon_range
        * offset:       nucleotide count offset (default: 30) 
    '''
    unique_amplicons = []
    for index in range(0, len(amplicons)):
        if index == 0:
            tmp_amplicon = [
                amplicons[index][0],
                amplicons[index][1],
                str(int(amplicons[index+1][1]) - int(offset)),
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        elif index == len(amplicons)-1:
            tmp_amplicon = [
                amplicons[index][0],
                str(int(amplicons[index-1][2]) + int(offset)),
                amplicons[index][2],
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        else:
            tmp_amplicon = [
                amplicons[index][0],
                str(int(amplicons[index-1][2]) + int(offset)),
                str(int(amplicons[index+1][1]) - int(offset)),
                amplicons[index][3],
                amplicons[index][4],
                amplicons[index][5]
            ]
        unique_amplicons.append(tmp_amplicon)
    return unique_amplicons
