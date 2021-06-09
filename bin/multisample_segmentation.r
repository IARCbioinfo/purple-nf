############################################
## Perform multi-sample segmentation      ##
############################################

#load libraries
library(readr)
library(copynumber)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)
library(data.table)

# get command line arguments
args = commandArgs(trailingOnly=TRUE) #ID, input folder location, output folder location, gamma param

# get chromosome arms location

if(args[5]=="NO_ARMS_FILE"){
  x <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", col_names = c("chrom","chromStart","chromEnd","name","gieStain"))
  x = x %>% filter(chrom %in% paste0("chr",c(1:22,"X","Y","M")))
  x$arm = sapply(x$name,function(a) str_split(a,"")[[1]][1])
  chromarms = x %>% group_by(chrom,arm) %>% summarize(armStart=min(chromStart),armEnd=max(chromEnd))
  write_tsv(chromarms,file = paste0(args[3],"chrarms.tsv") )
}else{
  chromarms = read_tsv(args[5])
}

# read baf
AMBERs       = lapply( list.files(args[2],recursive = T,pattern = "amber.baf.tsv",full.names = T), read_tsv)
AMBERs.names = str_remove(sapply(list.files(args[2],recursive = T,pattern = "amber.baf.tsv"), function(x){str_split(x,"/")[[1]][2]} ),".amber.baf.tsv")
for(i in 1:length(AMBERs)) colnames(AMBERs[[i]])[3:5] = paste0(colnames(AMBERs[[i]])[3:5],".",AMBERs.names[i])

# read RDR
COBALTs = lapply( list.files(args[2],recursive = T,pattern = "cobalt.ratio.tsv",full.names = T), read_tsv) 
COBALTs.names = str_remove(sapply(list.files(args[2],recursive = T,pattern = "cobalt.ratio.tsv"), function(x){str_split(x,"/")[[1]][2]} ),".cobalt.ratio.tsv")
for(i in 1:length(COBALTs)) colnames(COBALTs[[i]])[6:7] = paste0(colnames(COBALTs[[i]])[6:7],".",COBALTs.names[i])

# shape input
AMBER.pcfinput = as.data.frame(lapply(AMBERs,function(x) x[,c(1,2,4)]) %>% purrr::reduce(inner_join) )
COBALT.pcfinput = as.data.frame( lapply(COBALTs,function(x) x[,c(1,2,6,7)]) %>% purrr::reduce(inner_join) )

#determine chr arms
chromarms.GR       = GRanges(seqnames = chromarms$chrom,IRanges(start=chromarms$armStart,end=chromarms$armEnd),arm=chromarms$arm)
AMBER.pcfinput.GR  = GRanges(seqnames = AMBER.pcfinput$chromosome, IRanges(start=AMBER.pcfinput$position,end=AMBER.pcfinput$position))
COBALT.pcfinput.GR = GRanges(seqnames = COBALT.pcfinput$chromosome, IRanges(start=COBALT.pcfinput$position,end=COBALT.pcfinput$position))

AMBER.pcfinput.arms = rep(NA,nrow(AMBER.pcfinput))
AMBER.pcfinput.arms[queryHits(findOverlaps(AMBER.pcfinput.GR,chromarms.GR))] = chromarms.GR$arm[subjectHits(findOverlaps(AMBER.pcfinput.GR,chromarms.GR))]
COBALT.pcfinput.arms = rep(NA,nrow(COBALT.pcfinput))
COBALT.pcfinput.arms[queryHits(findOverlaps(COBALT.pcfinput.GR,chromarms.GR))] = chromarms.GR$arm[subjectHits(findOverlaps(COBALT.pcfinput.GR,chromarms.GR))]

if(! (all(!is.na(COBALT.pcfinput.arms)) & all(!is.na(AMBER.pcfinput.arms))) ) warning("Segments outside hg38 assembly!")

# format input
AMBER.pcfinput$chromosome = str_remove(AMBER.pcfinput$chromosome,"chr")
COBALT.pcfinput$chromosome = str_remove(COBALT.pcfinput$chromosome,"chr")

COBALT.pcfinput.arms = COBALT.pcfinput.arms[apply(COBALT.pcfinput[,str_detect(colnames(COBALT.pcfinput),"tumorGCRatio")]>0,1,all)]
COBALT.pcfinput = COBALT.pcfinput[apply(COBALT.pcfinput[,str_detect(colnames(COBALT.pcfinput),"tumorGCRatio")]>0,1,all),]
COBALT.pcfinputvals = COBALT.pcfinput[,str_detect(colnames(COBALT.pcfinput),"tumorGCRatio")]#/COBALT.pcfinput[,str_detect(colnames(COBALT.pcfinput),"referenceGCDiploidRatio")]
COBALT.pcfinputvals[COBALT.pcfinputvals<0.001] = 0.001
COBALT.pcfinputvals = log(COBALT.pcfinputvals,2)
COBALT.pcfinput = cbind(COBALT.pcfinput[,1:2] , COBALT.pcfinputvals )

#fit
AMBER.multipcfoutput  = multipcf(AMBER.pcfinput, gamma=as.numeric(args[4]),arms = AMBER.pcfinput.arms)
COBALT.multipcfoutput = multipcf(COBALT.pcfinput,gamma=as.numeric(args[4]),arms = COBALT.pcfinput.arms)
#COBALT.pcfoutput = pcf(COBALT.pcfinput[,-4],gamma=as.numeric(args[4]),kmin = 1)
#AMBER.pcfoutput = pcf(AMBER.pcfinput[,-4],gamma=as.numeric(args[4]),kmin = 1)

# filter size 0 segments
#AMBER.multipcfoutput = AMBER.multipcfoutput %>% filter(start.pos!=end.pos)
#COBALT.multipcfoutput = COBALT.multipcfoutput %>% filter(start.pos!=end.pos)

# format output
#AMBER.multipcfoutput$chrom  = paste0("chr",AMBER.multipcfoutput$chrom)
#COBALT.multipcfoutput$chrom = paste0("chr",COBALT.multipcfoutput$chrom)

# combine 
#AMBER.multipcfoutput.GR  = GRanges(seqnames = AMBER.multipcfoutput$chrom , 
#                                   ranges = IRanges(start=AMBER.multipcfoutput$start.pos,end=AMBER.multipcfoutput$end.pos), 
#                                   BAF = AMBER.multipcfoutput[,-(1:5)],n.probes=AMBER.multipcfoutput$n.probes )
#COBALT.multipcfoutput.GR = GRanges(seqnames = COBALT.multipcfoutput$chrom , 
#                                   ranges = IRanges(start=COBALT.multipcfoutput$start.pos,end=COBALT.multipcfoutput$end.pos),
#                                   RDR = 2**COBALT.multipcfoutput[,-(1:5)] , 
#                                   n.probes=COBALT.multipcfoutput$n.probes )

#ACintersect.GR = BiocGenerics::intersect(AMBER.multipcfoutput.GR,COBALT.multipcfoutput.GR)
#ovCint = findOverlaps(COBALT.multipcfoutput.GR,ACintersect.GR)
#ovAint = findOverlaps(AMBER.multipcfoutput.GR,ACintersect.GR)
#mcols(ACintersect.GR) = DataFrame(matrix(NA,length(ACintersect.GR),ncol(mcols(COBALT.multipcfoutput.GR)) + ncol(mcols(AMBER.multipcfoutput.GR))) )
#colnames(mcols(ACintersect.GR)) = c( colnames(mcols(COBALT.multipcfoutput.GR)) , colnames(mcols(AMBER.multipcfoutput.GR)) )
#mcols(ACintersect.GR[subjectHits(ovCint)])[,1:ncol(mcols(COBALT.multipcfoutput.GR))] = mcols(COBALT.multipcfoutput.GR)[queryHits(ovCint),]
#mcols(ACintersect.GR[subjectHits(ovAint)])[,(ncol(mcols(COBALT.multipcfoutput.GR))+1):(ncol(mcols(COBALT.multipcfoutput.GR))+ncol(mcols(AMBER.multipcfoutput.GR))) ] = mcols(AMBER.multipcfoutput.GR)[queryHits(ovAint),]
#ACintersect.GR = ACintersect.GR[!seqnames(ACintersect.GR)%in%c("chrX","chrY")]

# plot segments along the genome
#ggchrom <- ggplot(COBALT.multipcfoutput %>% pivot_longer(cols=colnames(COBALT.multipcfoutput)[str_detect(colnames(COBALT.multipcfoutput),"GCRatio")]) 
#                  %>% mutate(chrom=factor(chrom,levels=paste0("chr",c(1:22,"X","Y")))) , 
#       aes(x=start.pos,xend=end.pos,y=value,yend=value,col=name)) + geom_segment() + 
#  facet_wrap(.~chrom) + coord_cartesian(ylim=c(-1,1)) 

#ggsave(filename = paste0(args[3],args[1],"_segments.pdf"), ggchrom, height = 4*3, width = 4*3)

#write output in purple format
#out_formatted0 = bind_cols( "#CHR"=factor(as.character(seqnames(ACintersect.GR)),levels=paste0("chr",c(1:22,"X","Y") )) , 
#           START=start(ACintersect.GR), END=end(ACintersect.GR),
#           as_tibble(mcols(ACintersect.GR) )) 

#out_formatted0_RDR = out_formatted0[,c(1:(3+length(COBALTs.names)))] %>% 
#  pivot_longer(cols=colnames(out_formatted0)[str_detect(colnames(out_formatted0),"RDR")],names_pattern="[A-Za-z]*[.]([A-Za-z0-9_]*[0-9]+)",
#               names_to="SAMPLE",values_to="RD") 

#out_formatted0_BAF = out_formatted0[,c(1:3,(5+length(COBALTs.names)):(4+2*length(COBALTs.names)) )] %>% 
#  pivot_longer(cols=colnames(out_formatted0)[str_detect(colnames(out_formatted0),"BAF")],names_pattern="[A-Za-z]*[.]([A-Za-z0-9_]*[0-9]+)",
#               names_to="SAMPLE",values_to="BAF")
 
#out_formatted = left_join(out_formatted0_RDR,out_formatted0_BAF) %>% arrange(`#CHR`,START,END) %>% 
#  mutate("#SNPS" = 50,COV	= 60,ALPHA	= 10,BETA	= 1, BAF = 1-BAF) %>% relocate( `#SNPS`:BETA, .before=BAF)
for(i in 1:length(AMBERs)){
  dir.create(paste0(AMBERs.names[i],"_AMBER_multisampleseg") )
  dir.create(paste0(COBALTs.names[i],"_COBALT_multisampleseg") )
  # format AMBER output
  out_AMBER_tmp = AMBER.multipcfoutput %>% mutate(sampleID=AMBERs.names[i]) %>% 
    dplyr::select(sampleID,chrom:n.probes,colnames(AMBER.multipcfoutput)[str_detect(colnames(AMBER.multipcfoutput),AMBERs.names[i])])
  colnames(out_AMBER_tmp)[7] = "mean"
  write_tsv( out_AMBER_tmp, paste0(AMBERs.names[i],"_AMBER_multisampleseg/",AMBERs.names[i],".amber.baf.pcf") )
  
  # format COBALT output
  out_COBALT_tmp = COBALT.multipcfoutput %>% mutate(sampleID=COBALTs.names[i]) %>% 
    dplyr::select(sampleID,chrom:n.probes,colnames(COBALT.multipcfoutput)[str_detect(colnames(COBALT.multipcfoutput),COBALTs.names[i])])
  colnames(out_COBALT_tmp)[7] = "mean"
  write_tsv( out_COBALT_tmp, paste0(COBALTs.names[i],"_COBALT_multisampleseg/",COBALTs.names[i],".cobalt.ratio.pcf") )
}