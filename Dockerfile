################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="purple-nf"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for **purple-nf pipeline**"
LABEL about.home="http://github.com/IARCbioinfo/purple-nf"
LABEL about.documentation="http://github.com/IARCbioinfo/purple-nf/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/purple-nf/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER **digenovaa** <**digenovaa@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env update -n purple-nf -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/purple-nf/bin:$PATH
RUN conda env export --name purple-nf > purple-nf-v2.0.yml
RUN mkdir -p /hmftools/
COPY db /hmftools
RUN cd /hmftools/db/hg38/ && gzip -d DiploidRegions.38.bed.gz  && gzip -d  GC_profile.1000bp.38.cnp.gz && gzip -d  GermlineHetPon.38.vcf.gz  
#we get the code
#PURPLE
#ADD https://github.com/hartwigmedical/hmftools/releases/download/purple-v2.52/purple-2.52.jar /hmftools/jars
#COVALT
#ADD https://github.com/hartwigmedical/hmftools/releases/download/cobalt-v1.11/cobalt-1.11.jar /hmftools/jars
#RUN cd /home/programs && unzip master.zip
# amber
#ADD https://github.com/hartwigmedical/hmftools/releases/download/amber-v3.5/amber-3.5.jar /hmftools/jars
# databases
#ADD 


