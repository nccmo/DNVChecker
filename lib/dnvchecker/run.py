#!/usr/bin/env python

import pandas as pd
import subprocess, sys
import argparse
import listing, cis_check
import numexpr
import logging
import sys

def main(args):

    # logging
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    stdout_stream = logging.StreamHandler(sys.stdout)
    stdout_stream.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    stdout_stream.setFormatter(formatter)
    root.addHandler(stdout_stream)
    root.propagate = True
        
    listing.modify_df(args.MAF_file, args.output_file + ".tmp.txt", args.gene, args.sample)
    listing.DNV_filter(args.output_file + ".tmp.txt", args.output_file)

    listing.DNV_one_three_modification(args.output_file + ".DNV_one_three_to_remove.txt", args.output_file + ".DNV_one_three_to_remove.tmp.txt", args.ref_genome_path)
    listing.maf2maf(args.output_file + ".DNV_one_three_to_remove.tmp.txt", args.output_file + ".DNV_one_three_to_add.txt", 
                    args.ref_genome_path, args.vcf2maf_path, args.vep_path, args.vep_data_path, args.filter_vcf)

    listing.TNV_modification(args.output_file + ".TNV_to_remove.txt", args.output_file + ".TNV_to_remove.tmp.txt")
    listing.maf2maf(args.output_file + ".TNV_to_remove.tmp.txt", args.output_file + ".TNV_to_add.txt", 
                    args.ref_genome_path, args.vcf2maf_path, args.vep_path, args.vep_data_path, args.filter_vcf)

    listing.DNV_modification(args.output_file + ".DNV_to_remove.txt", args.output_file + ".DNV_to_remove.tmp.txt")
    if args.split == True:
        listing.file_split(args.output_file + ".DNV_to_remove.tmp.txt", args.output_file + ".DNV_to_remove.tmp.txt")
        listing.maf2maf2(args.output_file + ".DNV_to_remove.tmp.txt", args.output_file + ".DNV_to_add.txt", 
                         args.ref_genome_path, args.vcf2maf_path, args.vep_path, args.vep_data_path, args.filter_vcf)
        listing.file_merge(args.output_file + ".DNV_to_remove.tmp.txt", args.output_file + ".DNV_to_add.txt")
    if args.split == False:
        listing.maf2maf(args.output_file + ".DNV_to_remove.tmp.txt", args.output_file + ".DNV_to_add.txt", 
                 args.ref_genome_path, args.vcf2maf_path, args.vep_path, args.vep_data_path, args.filter_vcf)

    if args.sample_list != "":
        sample_check_DNV_one_three(args.output_file+ ".DNV_one_three_to_add.txt", args.output_file + ".DNV_one_three_to_add_read_number.txt", 
                                   args.output_file + ".DNV_one_three_to_remove.txt", args.output_file + ".DNV_one_three_to_remove_read_number.txt", 
                                   args.sample_list, args.debug)
        sample_check_TNV(args.output_file + ".TNV_to_add.txt", args.output_file + ".TNV_to_add_read_number.txt", 
                         args.output_file + ".TNV_to_remove.txt", args.output_file + ".TNV_to_remove_read_number.txt", 
                         args.sample_list, args.debug)
        if args.debug == False:
                subprocess.call("rm -f " + "*combined_cis.txt", shell = True)
                subprocess.call("rm -f " + "*combined_cis.txt.tmp1", shell = True)
                subprocess.call("rm -f " + "*combined_cis.txt.tmp2", shell = True)
                subprocess.call("rm -f " + "*combined_cis.txt.tmp3", shell = True)
        sample_check_DNV(args.output_file + ".DNV_to_add.txt", args.output_file + ".DNV_to_add_read_number.txt", 
                         args.output_file + ".DNV_to_remove.txt", args.output_file + ".DNV_to_remove_read_number.txt",
                         args.sample_list, args.debug)
        
    if args.debug == False:
        subprocess.call(["rm", "-f", args.output_file + ".tmp.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_one_three_to_remove.tmp.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_one_three_to_add.txt.vep.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".TNV_to_remove.tmp.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".TNV_to_add.txt.vep.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_to_remove.tmp.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_to_add.txt.vep.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_to_remove.tmp.txt"])
        subprocess.call(["rm", "-f", args.output_file + ".DNV_to_add.txt.vep.txt"])
        subprocess.call("rm -f " + args.output_file + ".DNV_to_remove.tmp.txt.split_*.txt", shell = True)
        subprocess.call("rm -f " + args.output_file + ".DNV_to_add.txt.split_*.txt", shell = True)
        subprocess.call("rm -f " + "*combined_cis.txt", shell = True)
        subprocess.call("rm -f " + "*combined_cis.txt.tmp1", shell = True)
        subprocess.call("rm -f " + "*combined_cis.txt.tmp2", shell = True)
        subprocess.call("rm -f " + "*combined_cis.txt.tmp3", shell = True)


def sample_check_DNV(input_file, output_file, input_file2, output_file2, sample_list, debug):

    # skip comment line(s)
    with open(sample_list) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    name_to_bam_df = pd.read_table(sample_list, sep='\t', skiprows=skip_rows, names=('name', 'BAM'))
    name_to_bam_dict = name_to_bam_df.set_index('name')['BAM'].to_dict()

    ref_ref = []
    alt_alt = []
    ref_alt = []
    alt_ref = []

    df = pd.read_table(input_file, sep='\t').fillna(0)
    for i in range(len(df)):
        if df.ix[i,'Tumor_Sample_Barcode'] not in name_to_bam_dict:
            ref_ref.append("NO BAM")
            ref_alt.append("NO BAM")
            alt_ref.append("NO BAM")
            alt_alt.append("NO BAM")
            continue

        is_ref_ref = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Reference_Allele'][1]])
        is_alt_alt = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Tumor_Seq_Allele2'][1]])
        is_ref_alt = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Tumor_Seq_Allele2'][1]])
        is_alt_ref = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Reference_Allele'][1]])

        cis_check.cis_check_two(name_to_bam_dict[df.ix[i,'Tumor_Sample_Barcode']], 
                                str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", 
                                df.ix[i,'Chromosome'], df.ix[i,'Start_Position'], df.ix[i,'End_Position'])
        cis_check.combine_two_results(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp1", 
                                      str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp2", 
                                      str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt")

        hhin = open(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", "r")
        head = hhin.readline().replace('\r', '').rstrip('\n').split('\n')
        ref_ref_num = 0
        alt_alt_num = 0
        ref_alt_num = 0
        alt_ref_num = 0
        for line in hhin:
            f = line.replace('\r', '').rstrip('\n').split('\t')
            if f[3] == is_ref_ref:
                ref_ref_num += 1
            if f[3] == is_alt_alt:
                alt_alt_num += 1
            if f[3] == is_ref_alt:
                ref_alt_num += 1
            if f[3] == is_alt_ref:
                alt_ref_num += 1
        hhin.close()

        ref_ref.append(ref_ref_num)
        alt_alt.append(alt_alt_num)
        ref_alt.append(ref_alt_num)
        alt_ref.append(alt_ref_num)

    df['ref_ref'] = ref_ref
    df['alt_alt'] = alt_alt
    df['ref_alt'] = ref_alt
    df['alt_ref'] = alt_ref
    df.to_csv(output_file, sep="\t", index=False)

    df2 = pd.read_table(input_file2, sep='\t')
    df2['ref_ref'] = sum([[i,i] for i in ref_ref],[])
    df2['alt_alt'] = sum([[i,i] for i in alt_alt],[])
    df2['ref_alt'] = sum([[i,i] for i in ref_alt],[])
    df2['alt_ref'] = sum([[i,i] for i in alt_ref],[])
    df2.to_csv(output_file2, sep="\t", index=False)   


def sample_check_DNV_one_three(input_file, output_file, input_file2, output_file2, sample_list, debug):

    # skip comment line(s)
    with open(sample_list) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    name_to_bam_df = pd.read_table(sample_list, sep='\t', skiprows=skip_rows, names=('name', 'BAM'))
    name_to_bam_dict = name_to_bam_df.set_index('name')['BAM'].to_dict()

    ref_ref = []
    alt_alt = []
    ref_alt = []
    alt_ref = []

    df = pd.read_table(input_file, sep='\t')
    for i in range(len(df)):
        if df.ix[i,'Tumor_Sample_Barcode'] not in name_to_bam_dict:
            ref_ref.append("NO BAM")
            ref_alt.append("NO BAM")
            alt_ref.append("NO BAM")
            alt_alt.append("NO BAM")
            continue

        is_ref_ref = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Reference_Allele'][2]])
        is_alt_alt = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Tumor_Seq_Allele2'][2]])
        is_ref_alt = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Tumor_Seq_Allele2'][2]])
        is_alt_ref = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Reference_Allele'][2]])

        cis_check.cis_check_two(name_to_bam_dict[df.ix[i,'Tumor_Sample_Barcode']], 
                                str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", 
                                df.ix[i,'Chromosome'], df.ix[i,'Start_Position'], df.ix[i,'End_Position'])
        cis_check.combine_two_results(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp1", 
                                      str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp2", 
                                      str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt")

        hhin = open(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", "r")
        head = hhin.readline().replace('\r', '').rstrip('\n').split('\n')
        ref_ref_num = 0
        alt_alt_num = 0
        ref_alt_num = 0
        alt_ref_num = 0
        for line in hhin:
            f = line.replace('\r', '').rstrip('\n').split('\t')
            if f[3] == is_ref_ref:
                ref_ref_num += 1
            if f[3] == is_alt_alt:
                alt_alt_num += 1
            if f[3] == is_ref_alt:
                ref_alt_num += 1
            if f[3] == is_alt_ref:
                alt_ref_num += 1
        hhin.close()

        ref_ref.append(ref_ref_num)
        alt_alt.append(alt_alt_num)
        ref_alt.append(ref_alt_num)
        alt_ref.append(alt_ref_num)

    df['ref_ref'] = ref_ref
    df['alt_alt'] = alt_alt
    df['ref_alt'] = ref_alt
    df['alt_ref'] = alt_ref
    df.to_csv(output_file, sep="\t", index=False)

    df2 = pd.read_table(input_file2, sep='\t')
    df2['ref_ref'] = sum([[i,i] for i in ref_ref],[])
    df2['alt_alt'] = sum([[i,i] for i in alt_alt],[])
    df2['ref_alt'] = sum([[i,i] for i in ref_alt],[])
    df2['alt_ref'] = sum([[i,i] for i in alt_ref],[])
    df2.to_csv(output_file2, sep="\t", index=False)   


def sample_check_TNV(input_file, output_file, input_file2, output_file2, sample_list, debug):

    # skip comment line(s)
    with open(sample_list) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    name_to_bam_df = pd.read_table(sample_list, sep='\t', skiprows=skip_rows, names=('name', 'BAM'))
    name_to_bam_dict = name_to_bam_df.set_index('name')['BAM'].to_dict()

    ref_ref_ref = []
    ref_ref_alt = []
    ref_alt_ref = []
    ref_alt_alt = []
    alt_ref_ref = []
    alt_ref_alt = []
    alt_alt_ref = []
    alt_alt_alt = []
    
    df = pd.read_table(input_file, sep='\t')
    for i in range(len(df)):
        if df.ix[i,'Tumor_Sample_Barcode'] not in name_to_bam_dict:
            ref_ref_ref.append("NO BAM")
            ref_ref_alt.append("NO BAM")
            ref_alt_ref.append("NO BAM")
            ref_alt_alt.append("NO BAM")
            alt_ref_ref.append("NO BAM")
            alt_ref_alt.append("NO BAM")
            alt_alt_ref.append("NO BAM")
            alt_alt_alt.append("NO BAM")
            continue

        is_ref_ref_ref = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Reference_Allele'][1], df.ix[i,'Reference_Allele'][2]])
        is_ref_ref_alt = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Reference_Allele'][1], df.ix[i,'Tumor_Seq_Allele2'][2]])
        is_ref_alt_ref = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Tumor_Seq_Allele2'][1], df.ix[i,'Reference_Allele'][2]])
        is_ref_alt_alt = ':'.join([df.ix[i,'Reference_Allele'][0], df.ix[i,'Tumor_Seq_Allele2'][1], df.ix[i,'Tumor_Seq_Allele2'][2]])
        is_alt_ref_ref = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Reference_Allele'][1], df.ix[i,'Reference_Allele'][2]])
        is_alt_ref_alt = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Reference_Allele'][1], df.ix[i,'Tumor_Seq_Allele2'][2]])
        is_alt_alt_ref = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Tumor_Seq_Allele2'][1], df.ix[i,'Reference_Allele'][2]])
        is_alt_alt_alt = ':'.join([df.ix[i,'Tumor_Seq_Allele2'][0], df.ix[i,'Tumor_Seq_Allele2'][1], df.ix[i,'Tumor_Seq_Allele2'][2]])

        cis_check.cis_check_three(name_to_bam_dict[df.ix[i,'Tumor_Sample_Barcode']], 
                                  str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", 
                                  df.ix[i,'Chromosome'], df.ix[i,'Start_Position'], df.ix[i,'Start_Position'] + 1, df.ix[i,'End_Position'])
        cis_check.combine_three_results(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp1", 
                                        str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp2", 
                                        str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt.tmp3", 
                                        str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt")

        hhin = open(str(df.ix[i,'Tumor_Sample_Barcode']) + str(df.ix[i,'Chromosome']) + str(df.ix[i,'Start_Position']) + "combined_cis.txt", "r")
        head = hhin.readline().replace('\r', '').rstrip('\n').split('\n')
        ref_ref_ref_num = 0
        ref_ref_alt_num = 0
        ref_alt_ref_num = 0
        ref_alt_alt_num = 0
        alt_ref_ref_num = 0
        alt_ref_alt_num = 0
        alt_alt_ref_num = 0
        alt_alt_alt_num = 0

        for line in hhin:
            f = line.replace('\r', '').rstrip('\n').split('\t')
            if f[4] == is_ref_ref_ref:
                ref_ref_ref_num += 1
            if f[4] == is_ref_ref_alt:
                ref_ref_alt_num += 1
            if f[4] == is_ref_alt_ref:
                ref_alt_ref_num += 1
            if f[4] == is_ref_alt_alt:
                ref_alt_alt_num += 1
            if f[4] == is_alt_ref_ref:
                alt_ref_ref_num += 1
            if f[4] == is_alt_ref_alt:
                alt_ref_alt_num += 1
            if f[4] == is_alt_alt_ref:
                alt_alt_ref_num += 1
            if f[4] == is_alt_alt_alt:
                alt_alt_alt_num += 1

        hhin.close()

        ref_ref_ref.append(ref_ref_ref_num)
        ref_ref_alt.append(ref_ref_alt_num)
        ref_alt_ref.append(ref_alt_ref_num)
        ref_alt_alt.append(ref_alt_alt_num)
        alt_ref_ref.append(alt_ref_ref_num)
        alt_ref_alt.append(alt_ref_alt_num)
        alt_alt_ref.append(alt_alt_ref_num)
        alt_alt_alt.append(alt_alt_alt_num)

    df['ref_ref_ref'] = ref_ref_ref
    df['ref_ref_alt'] = ref_ref_alt
    df['ref_alt_ref'] = ref_alt_ref
    df['ref_alt_alt'] = ref_alt_alt
    df['alt_ref_ref'] = alt_ref_ref
    df['alt_ref_alt'] = alt_ref_alt
    df['alt_alt_ref'] = alt_alt_ref
    df['alt_alt_alt'] = alt_alt_alt
    df.to_csv(output_file, sep="\t", index=False)

    df2 = pd.read_table(input_file2, sep='\t')
    df2['ref_ref_ref'] = sum([[i,i,i] for i in ref_ref_ref],[])
    df2['ref_ref_alt'] = sum([[i,i,i] for i in ref_ref_alt],[])
    df2['ref_alt_ref'] = sum([[i,i,i] for i in ref_alt_ref],[])
    df2['ref_alt_alt'] = sum([[i,i,i] for i in ref_alt_alt],[])
    df2['alt_ref_ref'] = sum([[i,i,i] for i in alt_ref_ref],[])
    df2['alt_ref_alt'] = sum([[i,i,i] for i in alt_ref_alt],[])
    df2['alt_alt_ref'] = sum([[i,i,i] for i in alt_alt_ref],[])
    df2['alt_alt_alt'] = sum([[i,i,i] for i in alt_alt_alt],[])
    df2.to_csv(output_file2, sep="\t", index=False)   