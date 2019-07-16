#!/usr/bin/env python

import pysam
import pandas as pd
import sys
import logging

logger = logging.getLogger(__name__)  # module logger

def cis_check_two(bam_file, output_file, chromosome, pos1, pos2):

    hout1 = open(output_file + ".tmp1", 'w')
    hout2 = open(output_file + ".tmp2", 'w')
    
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # check
    if int(pos1) == int(pos2):
        print >> sys.stderr, "positions are same."
        sys.exit(1)

    # pos1,2
    left_read = []
    tabixErrorFlag1 = 0
    try:
        records_one = bamfile.pileup(str(chromosome), int(pos1)-1, int(pos1)+1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag1 = 1

    if tabixErrorFlag1 == 0:
        for pileupcolumn in records_one:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and int(pileupcolumn.pos) == int(pos1)-1:
                    left_read.append(str(pileupread.alignment.query_name))
                
    right_read = []
    tabixErrorFlag2 = 0
    try:
        records_two = bamfile.pileup(str(chromosome), int(pos2)-1, int(pos2)+1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag2 = 1

    if tabixErrorFlag2 == 0:
        for pileupcolumn in records_two:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and (int(pileupcolumn.pos) == int(pos2)-1 and str(pileupread.alignment.query_name) in left_read):
                        print >> hout2, str(pileupread.alignment.query_name) + '\t' + pileupread.alignment.query_sequence[pileupread.query_position]
                        right_read.append(str(pileupread.alignment.query_name))

    if tabixErrorFlag1 == 0:
        for pileupcolumn in bamfile.pileup(str(chromosome), int(pos1)-1, int(pos1)+1):
            for pileupread in pileupcolumn.pileups:
                if int(pileupcolumn.pos) == int(pos1)-1:
                    if not pileupread.is_del and not pileupread.is_refskip and str(pileupread.alignment.query_name) in right_read:
                        print >> hout1, str(pileupread.alignment.query_name) + '\t' + pileupread.alignment.query_sequence[pileupread.query_position]

    bamfile.close()
    hout1.close()
    hout2.close()



def cis_check_three(bam_file, output_file, chromosome, pos1, pos2, pos3):

    hout1 = open(output_file + ".tmp1", 'w')
    hout2 = open(output_file + ".tmp2", 'w')
    hout3 = open(output_file + ".tmp3", 'w')

    bamfile = pysam.AlignmentFile(bam_file, "rb")

    # check
    if (int(pos1) == int(pos2)) or (int(pos1) == int(pos3)) or (int(pos2) == int(pos3)):
        print >> sys.stderr, "positions are same."
        sys.exit(1)

    # pos1,2,3
    left_read = []
    tabixErrorFlag1 = 0
    try:
        records_one = bamfile.pileup(str(chromosome), int(pos1)-1, int(pos1)+1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag1 = 1

    if tabixErrorFlag1 == 0:
        for pileupcolumn in records_one:
            for pileupread in pileupcolumn.pileups:
                if int(pileupcolumn.pos) == int(pos1)-1:
                    left_read.append(str(pileupread.alignment.query_name))


    right_read = []
    tabixErrorFlag2 = 0
    try:
        records_two = bamfile.pileup(str(chromosome), int(pos2)-1, int(pos2)+1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag2 = 1

    if tabixErrorFlag2 == 0:
        for pileupcolumn in records_two:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and (int(pileupcolumn.pos) == int(pos2)-1 and str(pileupread.alignment.query_name) in left_read):
                    right_read.append(str(pileupread.alignment.query_name))


    third_read = []
    tabixErrorFlag3 = 0
    try:
        records_three= bamfile.pileup(str(chromosome), int(pos3)-1, int(pos3)+1)
    except Exception as inst:
        print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
        tabixErrorFlag3 = 1

    if tabixErrorFlag3 == 0:
        for pileupcolumn in records_three:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and (int(pileupcolumn.pos) == int(pos3)-1 and str(pileupread.alignment.query_name) in right_read):
                        third_read.append(str(pileupread.alignment.query_name))
                        print >> hout3, str(pileupread.alignment.query_name) + '\t' + pileupread.alignment.query_sequence[pileupread.query_position]


    if tabixErrorFlag1 == 0:
        for pileupcolumn in bamfile.pileup(str(chromosome), int(pos1)-1, int(pos1)+1):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and int(pileupcolumn.pos) == int(pos1)-1:
                    if str(pileupread.alignment.query_name) in third_read:
                        print >> hout1, str(pileupread.alignment.query_name) + '\t' + pileupread.alignment.query_sequence[pileupread.query_position]

    if tabixErrorFlag2 == 0:
        for pileupcolumn in bamfile.pileup(str(chromosome), int(pos2)-1, int(pos2)+1):
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip and int(pileupcolumn.pos) == int(pos2)-1:
                    if str(pileupread.alignment.query_name) in third_read:
                        print >> hout2, str(pileupread.alignment.query_name) + '\t' + pileupread.alignment.query_sequence[pileupread.query_position]

    bamfile.close()
    hout1.close()
    hout2.close()
    hout3.close()


def combine_two_results(input_file1, input_file2, output_file):

    hout = open(output_file, 'w')
    print >> hout, '\t'.join(["read_ID", "pos1", "pos2", "merge"])

    dict2 = {}
    with open(input_file2, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            dict2[F[0]] = F[1]
    hin.close()

    with open(input_file1, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print >> hout, '\t'.join(F) + '\t' + dict2[F[0]] + '\t' + ":".join([F[1], dict2[F[0]]])
    hin.close()


def combine_three_results(input_file1, input_file2, input_file3, output_file):

    hout = open(output_file, 'w')
    print >> hout, '\t'.join(["read_ID", "pos1", "pos2", "pos3", "merge"])

    dict2 = {}
    with open(input_file2, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            dict2[F[0]] = F[1]
    hin.close()

    dict3 = {}
    with open(input_file3, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            dict3[F[0]] = F[1]
    hin.close()

    with open(input_file1, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print >> hout, '\t'.join(F) + '\t' + dict2[F[0]] + '\t' + dict3[F[0]] + '\t' + ":".join([F[1], dict2[F[0]], dict3[F[0]]])
    hin.close()