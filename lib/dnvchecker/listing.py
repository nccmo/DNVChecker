#! /usr/bin/env python

import pandas as pd
from numpy import mean
import numpy as np
import csv
import re
import os
import subprocess
import sys
import logging

logger = logging.getLogger(__name__)  # module logger

def modify_df(maf_file, output_file, gene, sample):

#####################
# function to modify MAF dataframe (select only necessary variables) and 
#####################

    # skip comment line(s)
    with open(maf_file) as handle:
        first_line = next(handle)
        skip_rows = 1 if first_line.startswith('#') else 0

    # read in data frame
    df = pd.read_table(maf_file, sep='\t', skiprows=skip_rows)

    # select only SNPs
    df = df[(df['Variant_Type'] == "SNP")]

    # select gene(s)
    if gene != "":
        gene_list = gene.split(':')
        df = df[(df['Hugo_Symbol'].isin(gene_list))]

    # select sample(s)
    if sample != "":
        sample_list = sample.split(':')
        df = df[(df['Tumor_Sample_Barcode'].isin(sample_list))]

    #rename column names
    df = df.rename(columns={'amino_acid_change': 'HGVSp_Short'})
    df = df.rename(columns={'aa_change': 'HGVSp_Short'})
    df = df.rename(columns={'aachange': 'HGVSp_Short'})
    df = df.rename(columns={'Protein_Change': 'HGVSp_Short'})
    df = df.rename(columns={'Start_position': 'Start_Position'})
    df = df.rename(columns={'End_position': 'End_Position'})
    if 'Tumor_Seq_Allele' in df.columns.values:
        df = df.rename(columns={'Tumor_Seq_Allele': 'Tumor_Seq_Allele2'})

    # error check
    # check if necessary columns are in the MAF file
    variables = ['Hugo_Symbol', 'Variant_Type', 'Variant_Type', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Chromosome', 
                 'Start_Position', 'End_Position', 'Strand', 'Reference_Allele', 'Tumor_Seq_Allele2', 'HGVSp_Short']
    if not all((s in df.columns.values) for s in variables):
        logger.error('necessary columns are not in the list')
        exit(1)

    df_s = df.sort(['Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Chromosome', 'Start_Position', 'Strand'])
    df_s.to_csv(output_file, sep='\t', index = False)


def DNV_filter(input_file, output_file):

#####################
# function to select DNV candidates
#####################

   df = pd.read_table(input_file, sep='\t').fillna(0)
   df['aa_pos'] = df['HGVSp_Short'].str.extract(r'([0-9]+)')

   df['pos_diff'] = (df['Start_Position'] - df['Start_Position'].shift(1))
   df['prev_chr'] = df['Chromosome'].shift(1)
   df['prev_sample'] = df['Tumor_Sample_Barcode'].shift(1)
   df['prev_norm_sample'] = df['Matched_Norm_Sample_Barcode'].shift(1)
   df['prev_aa_pos'] = df['aa_pos'].shift(1)
   
   df['pos_diff2'] = (df['Start_Position'] - df['Start_Position'].shift(2))
   df['prev_chr2'] = df['Chromosome'].shift(2)
   df['prev_sample2'] = df['Tumor_Sample_Barcode'].shift(2)
   df['prev_norm_sample2'] = df['Matched_Norm_Sample_Barcode'].shift(2)
   df['prev_aa_pos2'] = df['aa_pos'].shift(2)

   ############################
   # TNV
   x = df[(df['pos_diff'] == 1) &
          (df['pos_diff2'] == 2) & 
          (df['Chromosome'] == df['prev_chr']) & (df['prev_chr'] == df['prev_chr2']) & 
          (df['Tumor_Sample_Barcode'] == df['prev_sample']) & (df['prev_sample'] == df['prev_sample2']) &
          (df['Matched_Norm_Sample_Barcode'] == df['prev_norm_sample']) & (df['prev_norm_sample'] == df['prev_norm_sample2']) &
          (df['aa_pos'] == df['prev_aa_pos']) & (df['prev_aa_pos'] == df['prev_aa_pos2'])]

   indx = list(x.index)
   indx2 = [i-1 for i in indx]
   indx3 = [i-2 for i in indx]
   output_df = df[df.index.isin(indx) | df.index.isin(indx2) | df.index.isin(indx3)]

   output_df = output_df.drop(['aa_pos', 'pos_diff', 'prev_chr', 'prev_sample', 'prev_norm_sample', 'prev_aa_pos',
                               'pos_diff2', 'prev_chr2', 'prev_sample2', 'prev_norm_sample2', 'prev_aa_pos2'], axis=1)
   output_df.to_csv(output_file + ".TNV_to_remove.txt", sep = "\t", index = False)

   ############################
   # DNV one three
   y = df[(df['pos_diff'] == 2) &
          (df['Chromosome'] == df['prev_chr']) & 
          (df['Tumor_Sample_Barcode'] == df['prev_sample']) &
          (df['Matched_Norm_Sample_Barcode'] == df['prev_norm_sample']) &
          (df['aa_pos'] == df['prev_aa_pos'])]

   indy = list(y.index)
   indy2 = [i-1 for i in indy]
   output_df2 = df[df.index.isin(indy) | df.index.isin(indy2)]

   output_df2 = output_df2.drop(['aa_pos', 'pos_diff', 'prev_chr', 'prev_sample', 'prev_norm_sample', 'prev_aa_pos',
                                 'pos_diff2', 'prev_chr2', 'prev_sample2', 'prev_norm_sample2', 'prev_aa_pos2'], axis=1)
   output_df2.to_csv(output_file + ".DNV_one_three_to_remove.txt", sep = "\t", index = False)

   ############################
   # DNV
   z = df[(df['pos_diff'] == 1) & 
          (df['Chromosome'] == df['prev_chr']) & 
          (df['Tumor_Sample_Barcode'] == df['prev_sample']) &
          (df['Matched_Norm_Sample_Barcode'] == df['prev_norm_sample']) &
          (df['aa_pos'] == df['prev_aa_pos'])]

   indz = list(z.index)
   indz2 = [i-1 for i in indz]
   output_df3 = df[df.index.isin(indz) | df.index.isin(indz2)]

   # remove TNV
   output_df3 = output_df3[~output_df3.index.isin(indx) & ~output_df3.index.isin(indx2) & ~output_df3.index.isin(indx3)]
   output_df3 = output_df3.drop(['aa_pos', 'pos_diff', 'prev_chr', 'prev_sample', 'prev_norm_sample', 'prev_aa_pos',
                                'pos_diff2', 'prev_chr2', 'prev_sample2', 'prev_norm_sample2', 'prev_aa_pos2'], axis=1)
   output_df3.to_csv(output_file + ".DNV_to_remove.txt", sep = "\t", index = False)


def DNV_modification(input_file, output_file):

    df = pd.read_table(input_file, sep='\t')
    df['prev_ref'] = df['Reference_Allele'].shift(1)
    df['prev_alt'] = df['Tumor_Seq_Allele2'].shift(1)

    df['Start_Position'] = df['Start_Position']-1
    df['Reference_Allele'] = df['prev_ref'] + df['Reference_Allele'] 
    df['Tumor_Seq_Allele2'] = df['prev_alt'] + df['Tumor_Seq_Allele2']
   
    df['Variant_Type'] = ["DNV"] * len(df)
    df = df[1::2]

    df = df.drop(['prev_ref', 'prev_alt'], axis=1)
    df.to_csv(output_file, sep="\t", index=False)


def TNV_modification(input_file, output_file):

    df = pd.read_table(input_file, sep='\t')
    df['prev_ref'] = df['Reference_Allele'].shift(1)
    df['prev_alt'] = df['Tumor_Seq_Allele2'].shift(1)
    df['prev_ref2'] = df['Reference_Allele'].shift(2)
    df['prev_alt2'] = df['Tumor_Seq_Allele2'].shift(2)

    df['Start_Position'] = df['Start_Position']-2
    df['Reference_Allele'] = df['prev_ref2'] + df['prev_ref'] + df['Reference_Allele'] 
    df['Tumor_Seq_Allele2'] = df['prev_alt2'] + df['prev_alt'] + df['Tumor_Seq_Allele2']
   
    df['Variant_Type'] = ["TNV"] * len(df)
    df = df[2::3]

    df = df.drop(['prev_ref', 'prev_alt', 'prev_ref2', 'prev_alt2'], axis=1)
    df.to_csv(output_file, sep="\t", index=False)


def DNV_one_three_modification(input_file, output_file, ref_genome):

    df = pd.read_table(input_file, sep='\t')
    df['prev_ref'] = df['Reference_Allele'].shift(1)
    df['prev_alt'] = df['Tumor_Seq_Allele2'].shift(1)
    df['Start_Position'] = df['Start_Position']-2

    mid_aa = []

    for i in range(len(df)):
        inp = str(df.ix[i,'Chromosome']) +  ":" + str(int(df.ix[i, 'Start_Position'])+1) + "-" + str(int(df.ix[i, 'Start_Position'])+1)
        output = subprocess.Popen("samtools faidx %s %s" %(ref_genome, inp), shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        line = output.communicate()[0].rstrip().split("\n")
        mid_aa.append(line[1])

    df['mid_aa'] = mid_aa

    df['Reference_Allele'] = df['prev_ref'] + df['mid_aa'] + df['Reference_Allele'] 
    df['Tumor_Seq_Allele2'] = df['prev_alt'] + df['mid_aa'] + df['Tumor_Seq_Allele2']
    df['Variant_Type'] = ["TNV"] * len(df)
    df = df[1::2]

    df = df.drop(['prev_ref', 'prev_alt', 'mid_aa'], axis=1)
    df.to_csv(output_file, sep="\t", index=False)


def maf2maf(input_file, output_file, ref, path, vep_path, vep_data_path, filter_vcf):

    hin = open(input_file, 'r')
    cols = hin.readline().strip('\n').split('\t')
    
    input_df = pd.read_csv(input_file, sep='\t')
    cols = input_df.columns.values

    if len(input_df) != 0:    
      sret = subprocess.Popen("perl %s/maf2maf.pl --input-maf %s --output-maf %s --ref-fasta %s --vep-path %s --vep-data %s --filter-vcf %s" 
                              % (path, input_file, output_file + ".vep.txt", ref, vep_path, vep_data_path, filter_vcf), shell = True).wait()
      if sret != 0:
          print >> sys.stderr, "vcf2maf error, error code: " + str(sret)

      df = pd.read_table(output_file + ".vep.txt", sep='\t', skiprows=1)
      cols2 = [value for value in cols if value in df.columns] 
      df = df[cols2]

      dropped_input_df = input_df.drop(cols2, axis=1)
      df = pd.concat([df, dropped_input_df], axis=1)

      df.to_csv(output_file, sep="\t", index=False)

    else:
      input_df.to_csv(output_file, sep="\t", index=False)

def file_split(input_file, output_file):

    reader =  pd.read_csv(input_file, sep='\t', chunksize=100)
    chunk_number = -(-(sum(1 for row in open(input_file, 'r')) - 1) // 100)

    n = 0
    
    for i in range(int(chunk_number)):
        df = pd.DataFrame(reader.get_chunk(100))
        df.to_csv(output_file + ".split_" + str(i) + ".txt", sep="\t", index=False)


def maf2maf2(input_file, output_file, ref, path, vep_path, vep_data_path, filter_vcf):

    chunk_number = -(-(sum(1 for row in open(input_file, 'r')) - 1) // 100)

    ps = []
    for i in range(int(chunk_number)):
        sret = subprocess.Popen("perl %s/maf2maf.pl --input-maf %s --output-maf %s --ref-fasta %s --vep-path %s --vep-data %s --filter-vcf %s" 
                                % (path, input_file + ".split_" + str(i) + ".txt", output_file + ".split_" + str(i) + ".vep.txt", ref, vep_path, vep_data_path, filter_vcf), shell = True)
        ps.append(sret)
        if sret != 0:
          print >> sys.stderr, "vcf2maf error, error code: " + str(sret)

    for p in ps:
        p.wait()

def file_merge(input_file, output_file):

  chunk_number = -(-(sum(1 for row in open(input_file, 'r')) - 1) // 200)
  input_df = pd.read_csv(input_file, sep='\t')
  cols = input_df.columns.values

  merged_df = pd.read_table(output_file + ".split_0.vep.txt", sep='\t', skiprows=1)
  cols2 = [value for value in cols if value in merged_df.columns] 
  merged_df = merged_df[cols2]

  if chunk_number >= 2:
    for i in range(1, int(chunk_number)):
      df = pd.read_table(output_file + ".split_" + str(i) + ".vep.txt", sep='\t', skiprows=1)
      df = df[cols2]
      merged_df =merged_df.append(df,ignore_index=True)

  dropped_input_df = input_df.drop(cols2, axis=1)
  merged_df = pd.concat([merged_df, dropped_input_df], axis=1)
  merged_df.to_csv(output_file, sep="\t", index=False)
