# DNVChecker

## Introduction
DNVChecker is a software to count DNVs/TNVs from Mutation Annotation Format (MAF) input files. 
If a list of BAM files is provided, DNVChecker can return the read counts with ref/alt allele.

## Dependency
Python
Python (2.7), pysam (0.14.1), pandas(0.16.2) packages
vcf2maf (v1.6.16), Ensembl ’s Variant Effect Predictor (v95)
(Other versions are not checked)

## Input formats
### Mutations
Mutations are provided in a Mutation Annotation Format (MAF) file.

### Sample list of BAM files (optional)
If list of BAM files are provied, DNVChecker can check if DNV/TNV occurs in cis or in trans.
Sample list needs be tab-delimited.

example format
```
file1 /home/data/file1.bam
file2 /home/data/file2.bam
file3 /home/data/file3.bam
```

## How to use
Before usage, type the sentence below.
```
python setup.py install --user
```

```
$ dnvchecker -h
usage: dnvchecker [-h] [--version] [--gene GENE] [--sample SAMPLE] [--vcf2maf_path VCF2MAF_PATH] [--ref_genome_path REF_GENOME_PATH] [--vep_path VEP_PATH] [--vep_data_path VEP_DATA_PATH] [--filter_vcf FILTER_VCF] [--debug] [--sample_list SAMPLE_LIST] [--split SPLIT]
                  input.maf output.txt

positional arguments:
  input.maf             the path to input MAF file
  output.txt            the path to output file

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --gene GENE           specify the gene you want to check (if you want to choose multiple genes, separate them with colon) (default: all genes)
  --sample SAMPLE       specify the sample you want to check (if you want to choose multiple samples, separate them with colon) (default: all samples)
  --vcf2maf_path VCF2MAF_PATH
                        the path to vcf2maf folder (default: /home/ysaito/bin/vcf2maf-1.6.16)
  --ref_genome_path REF_GENOME_PATH
                        the path to reference genome (default: /home/ysaito/.vep/homo_sapiens/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz)
  --vep_path VEP_PATH   the path to Variant Effect Predictor(VEP) (default: /home/ysaito/vep)
  --vep_data_path VEP_DATA_PATH
                        the path to Variant Effect Predictor(VEP) data (default: /home/ysaito/.vep)
  --filter_vcf FILTER_VCF
                        A VCF for FILTER tag common_variant (default: /home/ysaito/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz)
  --debug               keep intermediate files
  --sample_list SAMPLE_LIST
                        the path to sample list
  --split SPLIT         use split mode (appropriate for large dataset)

```

Split mode is appropriate for large dataset (ex. TCGA pan-cancer dataset).

## Test running
```
$ dnvchecker examples/TCGA_sample_file.txt
```


## Output files
Six files will be provided.
```
output_file + ".DNV_to_remove.txt"
output_file + ".DNV_to_add.txt"
output_file + ".DNV_one_three_to_remove.txt"
output_file + ".DNV_one_three_to_add.txt"
output_file + ".TNV_to_remove.txt"
output_file + ".TNV_to_add.txt"
```

if sample list is provided, three additional files will be provided.
```
output_file + ".DNV_to_add_add_read_number.txt"
output_file + ".DNV_one_three_to_add_read_number.txt"
output_file + ".TNV_to_add_add_read_number.txt"
```