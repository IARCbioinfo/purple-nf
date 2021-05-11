# purple-nf
Nextflow pipeline for CNV calling with PURPLE

## Description
Pipeline using [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) for copy number calling from tumor/normal  or tumor-only sequencing data.

## Usage
  ```
  #using a tn_pairs file
  nextflow run iarcbioinfo/purple-nf -r v1.0 \
  -profile singularity  --tn_file tn_pairs..txt \
  --cohort_dir $PWD/CRAM \
  --ref hs38DH.fa --ref_dict hs38DH.dict \
  --output_folder PURPLE
  
  
  #activate BAM files mode
  nextflow run iarcbioinfo/purple-nf -r v1.0 \
  -profile singularity  --tn_file tn_pairs..txt \
  --cohort_dir $PWD/BAM \
  --bam \
  --ref hs38DH.fa --ref_dict hs38DH.dict \
  --output_folder PURPLE
 
  ```

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
	- [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
	- [COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt)
	- [AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber)
	
You can avoid installing all the external software by only installing Docker or singularity.
See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input (mandatory)

  | Type      | Description   |
  |-----------|---------------|
  | --cohort_dir    | Folder containing all BAM/CRAM files |  
  | --tn_file    | File containing the list of names of BAM files to be processed |
  |--ref         |  Fasta file of reference genome [hg38.fa], should be indexed [hg38.fa.fai]|
  | --ref_dict  | dict file for the reference genome [hg38.dict]|


### Example of Tumor/Normal pairs file (--tn_file)
A text file tabular separated, with the following header:

```
tumor_id	sample	tumor	normal
sample1_T1	sample1	sample1_T.cram	sample1_N.cram
sample2_T1	sample2	sample2_T.cram	sample2_N.cram
sample3_T1	sample3	sample3_T.cram	sample3_N.cram
``` 

### Optional parameters

| Name      | type | Description     |
|-----------|---------------|-----------------|
| --tumor_only |         [flag] | active tumor_only mode|
|      --bam     |       [flag] |active bam mode [def:cram]|
|     --output_folder |  [string] |name of output folder |
|      --cpu          |[Integer] | Number of CPUs[def:2] |
|      --mem |        [Integer] | Max memory [def:8Gb] |  



## Output

```
results
├── AMBER                               # AMBER result directory
│   ├── S00016_T_AMBER
│   │   ├── amber.version
│   │   ├── S00016_T_N.amber.snp.vcf.gz
│   │   ├── S00016_T_N.amber.snp.vcf.gz.tbi
│   │   ├── S00016_T_T.amber.baf.pcf
│   │   ├── S00016_T_T.amber.baf.tsv
│   │   ├── S00016_T_T.amber.baf.vcf.gz
│   │   ├── S00016_T_T.amber.baf.vcf.gz.tbi
│   │   ├── S00016_T_T.amber.contamination.tsv
│   │   ├── S00016_T_T.amber.contamination.vcf.gz
│   │   ├── S00016_T_T.amber.contamination.vcf.gz.tbi
│   │   └── S00016_T_T.amber.qc
├── COBALT									# COBALT result directory	
│   ├── S00016_T_COBALT
│   │   ├── cobalt.version
│   │   ├── S00016_T_N.cobalt.gc.median.tsv
│   │   ├── S00016_T_N.cobalt.ratio.median.tsv
│   │   ├── S00016_T_N.cobalt.ratio.pcf
│   │   ├── S00016_T_T.chr.len
│   │   ├── S00016_T_T.cobalt.gc.median.tsv
│   │   ├── S00016_T_T.cobalt.ratio.pcf
│   │   └── S00016_T_T.cobalt.ratio.tsv
│   ├── .....
├── PURPLE									# PURPLE result directory	
│   ├── S00016_T_PURPLE
│   │   ├── circos							# circos direcoty with files for plotting
│   │   ├── purple.version
│   │   ├── S00016_T_T.purple.cnv.gene.tsv
│   │   ├── S00016_T_T.purple.cnv.germline.tsv
│   │   ├── S00016_T_T.purple.cnv.somatic.tsv.           # Somatic Copy Number Segments
│   │   ├── S00016_T_T.purple.purity.range.tsv
│   │   ├── S00016_T_T.purple.purity.tsv
│   │   ├── S00016_T_T.purple.qc
│   │   ├── S00016_T_T.purple.segment.tsv
│   │   └── S00016_T_T.purple.somatic.clonality.tsv
│   ├── .....    
└── purple_summary.txt     # Summary file for all tumors
├── nf-pipeline_info		# Nextflow information directory
│   ├── purple_dag.html
│   ├── purple_report.html
│   ├── purple_timeline.html
│   ├── purple_trace.txt
│   └── run_parameters_report.txt # Custom file providing info for software versions and calling parameters
```


## Common errors

### Singularity
The first time that the container is built from the docker image, the TMPDIR  should be defined in a non parallel file-system, you can set this like:

```
export TMPDIR=/tmp
```

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    |            follm@iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Alex Di Genova | digenovaa@fellows.iarc.fr| Developer |
