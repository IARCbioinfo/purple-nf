.DELETE_ON_ERROR:

AMBER=/home/digenovaa/binaries/hmftools/amber-3.5.jar
COBALT=/home/digenovaa/binaries/hmftools/cobalt-1.11.jar
PURPLE=/home/digenovaa/binaries/hmftools/purple-2.52.jar

#additional files for run
DR=/home/digenovaa/binaries/hmftools/database/DiploidRegions.hg38.bed
GC=/home/digenovaa/binaries/hmftools/database/GC_profile.hg38.1000bp.cnp
GHPC=/home/digenovaa/binaries/hmftools/database/GermlineHetPon.hg38.snpcheck.vcf
GHP=/home/digenovaa/binaries/hmftools/database/GermlineHetPon.hg38.vcf
#REFERENCE GENOME
REF=/data/gcs/mesomics/files/references/hs38DH/hs38DH.fa
CPU=8

TC=CRAM/B00JAM2_alt_HMJV2CCXY.DUAL155.3-HMJV2CCXY.DUAL155.2_merged_BQSRecalibrated.cram
NC=CRAM/B00JAM1_alt_HMJV2CCXY.DUAL154.1_BQSRecalibrated.cram
PREFIX=MESO_001_T
#COBALT
#paired mode
${PREFIX}_COVALT/${PREFIX}_T.cobalt.ratio.pcf:
	java -jar ${COBALT} -gc_profile ${GC} -ref_genome ${REF} \
	 -reference ${PREFIX}_N -reference_bam ${NC} \
	 -tumor ${PREFIX}_T -tumor_bam ${TC} \
	 -output_dir ${PREFIX}_COVALT -threads ${CPU}

#tumor-only mode
${PREFIX}_COVALT_TON/${PREFIX}_T.cobalt.ratio.pcf: ${PREFIX}_COVALT/${PREFIX}_T.cobalt.ratio.pcf
	java -jar ${COBALT} -gc_profile ${GC} -ref_genome ${REF} -tumor_only -tumor_only_diploid_bed ${DR} \
	  -tumor  ${PREFIX}_T -tumor_bam ${TC} -output_dir ${PREFIX}_COVALT_TON -threads ${CPU}


#AMBER
#paired mode
${PREFIX}_AMBER/${PREFIX}_T.amber.baf.pcf: ${PREFIX}_COVALT_TON/${PREFIX}_T.cobalt.ratio.pcf
	java -jar ${AMBER} -loci ${GHP} -ref_genome ${REF} \
	 -reference ${PREFIX}_N -reference_bam ${NC} \
	 -tumor ${PREFIX}_T -tumor_bam ${TC} \
	 -output_dir ${PREFIX}_AMBER -threads ${CPU}

#t-only mode
${PREFIX}_AMBER_TON/${PREFIX}_T.amber.baf.pcf: ${PREFIX}_AMBER/${PREFIX}_T.amber.baf.pcf
	java -jar ${AMBER} -loci ${GHP} -ref_genome ${REF} -tumor_only \
	 -tumor ${PREFIX}_T -tumor_bam ${TC} -output_dir ${PREFIX}_AMBER_TON -threads ${CPU}

#PURPLE
#paired mode 
${PREFIX}_PURPLE/${PREFIX}_T.purple.cnv.gene.tsv:${PREFIX}_AMBER_TON/${PREFIX}_T.amber.baf.pcf
	java -jar ${PURPLE} \
	 -reference ${PREFIX}_N  -tumor ${PREFIX}_T \
	 -no_charts \
	 -output_dir ${PREFIX}_PURPLE \
	 -amber ${PREFIX}_AMBER \
	 -cobalt ${PREFIX}_COVALT \
	 -gc_profile ${GC} \
	 -threads ${CPU} \
	 -ref_genome ${REF}

#t-only mode 
${PREFIX}_PURPLE_TON/${PREFIX}_T.purple.cnv.gene.tsv:${PREFIX}_PURPLE/${PREFIX}_T.purple.cnv.gene.tsv
	java -jar ${PURPLE} \
	 -tumor ${PREFIX}_T \
	 -tumor_only -no_charts \
	 -output_dir ${PREFIX}_PURPLE_TON \
	 -amber ${PREFIX}_AMBER_TON \
	 -cobalt ${PREFIX}_COVALT_TON \
	 -gc_profile ${GC} \
	 -threads ${CPU} \
	 -ref_genome ${REF}
	

all: ${PREFIX}_PURPLE_TON/${PREFIX}_T.purple.cnv.gene.tsv
