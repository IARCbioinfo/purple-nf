#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/purple-nf -singularity [OPTIONS]

    Mandatory arguments:
      --tn_file		         [file] File containing list of T/N bam/cram files to be processed
      --cohort_dir         [dir]  directory where the BAM or CRAM  file are stored
      --ref                [file] fasta file of reference genome [hg38.fa]
      --ref_dict           [file] dict file for the reference genomep [hg38.dict]
    Optional arguments:
      --tumor_only         [flag] active tumor_only mode
      --bam                 [flag] active bam mode [def:cram]
      --output_folder       [string] name of output folder

      """.stripIndent()
}





//we display help information
if (params.help){ show_help(); exit 0;}
//we display the header of the tool
log.info IARC_Header()
log.info tool_header()
//Check mandatory parameters
assert (params.ref != null) : "please specify --ref reference.fasta"
assert (params.ref_dict != null) : "please specify --ref_dict reference.dict"
assert (params.tn_file != null ) : "please specify --tn_file"
assert (params.cohort_dir != null ) : "please specify --cohort_dir"

//if(params.tn_file != null && params.cohort_dir == null){ println "--cohort_dir shold be specified when using the --tn_file variable"; exit 1;}
//function that read the tumors to process from a tn_file
if(params.tn_file){
  def cram = params.bam ? false:true
 tn_pairs = parse_tn_file(params.tn_file,params.cohort_dir,cram)
 //we duplicate the tn_pairs channel
 tn_pairs.into { tn_pairs_cobalt; tn_pairs_amber}
}
//chanel for VCF file
//ch_vcf = Channel.value(file(params.dbsnp_vcf_ref)).ifEmpty{exit 1, "VCF file not found: ${params.dbsnp_vcf_ref}"}
ref_fasta = Channel.value(file(params.ref)).ifEmpty{exit 1, "reference file not found: ${params.ref}"}
ref_dict = Channel.value(file(params.ref_dict)).ifEmpty{exit 1, "Dict reference file not found: ${params.ref_dict}"}


print_params()

//PATHS in the container for databases
// /hmftools/hg38

// /hmftools/hg38/DiploidRegions.38.bed
// /hmftools/hg38/GC_profile.1000bp.38.cnp
// /hmftools/hg38/GermlineHetPon.38.vcf

process COBALT {

  publishDir params.output_folder+'/COBALT/', mode: 'copy'
  input:
  set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_cobalt
  file(ref) from ref_fasta
  file(dict) from ref_dict
  output:
  set val(tumor_id), path("${tumor_id}_COBALT") into cobalt
  script:
     if(params.tumor_only == false){
       """
      COBALT -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \
              -ref_genome ${ref} -tumor_only -tumor_only_diploid_bed /hmftools/hg38/DiploidRegions.38.bed \
     	        -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_COBALT -threads 1
       """
     }else{
       """
      COBALT -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \
              -ref_genome ${ref} -reference ${tumor_id}_N -reference_bam ${normal} \
     	        -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_COBALT -threads 1
      """
     }

}

process AMBER {

  publishDir params.output_folder+'/AMBER/', mode: 'copy'
  input:
  set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_amber
  file(ref) from ref_fasta
  file(dict) from ref_dict
  output:
  set val(tumor_id), path("${tumor_id}_AMBER") into amber
  script:
     if(params.tumor_only == false){
       """
      AMBER  -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} -tumor_only \
              -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads 1
      """
     }else{
       """
       AMBER  -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} \
               -reference ${tumor_id}_N -reference_bam ${normal}  \
               -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads 1
        """
     }

}
//we merge previous results from amber and cobalt
amber_cobalt=amber.join(cobalt, remainder: true)

process PURPLE {

  publishDir params.output_folder+'/PURPLE/', mode: 'copy'
  input:
  //set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_amber
  set val(tumor_id), path(amber_dir), path(cobalt_dir) from amber_cobalt
  file(ref) from ref_fasta
  file(dict) from ref_dict
  output:
  set val(tumor_id), path("${tumor_id}_PURPLE") into purple
  //MESO_071_T_T.purple.purity.tsv
  //set val(tumor_id), file("${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv") into stats_purple
  file("${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.sample.tsv") into stats_purple

  script:
     if(params.tumor_only == false){
       """
       PURPLE -reference ${tumor_id}_N  -tumor ${tumor_id}_T \
              -no_charts \
              -output_dir ${tumor_id}_PURPLE \
              -amber ${amber_dir} \
              -cobalt ${cobalt_dir} \
              -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \
              -threads 1 \
              -ref_genome ${ref}

              awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.sample.tsv
       """
     }else{
       """
       PURPLE  -tumor_only  -tumor ${tumor_id}_T \
               -no_charts \
               -output_dir ${tumor_id}_PURPLE \
               -amber ${amber_dir} \
               -cobalt ${cobalt_dir} \
               -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \
               -threads 1 \
               -ref_genome ${ref}

        awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.sample.tsv

       """
     }

}

stats_purple.collectFile(name: 'purple_summary.txt', storeDir: params.output_folder, seed: 'tumor_id\tpurity\tnormFactor\tscore\tdiploidProportion\tploidy\tgender\tstatus\tpolyclonalProportion\tminPurity\tmaxPurity\tminPloidy\tmaxPloidy\tminDiploidProportion\tmaxDiploidProportion\tversion\tsomaticPenalty\twholeGenomeDuplication\tmsIndelsPerMb\tmsStatus\ttml\ttmlStatus\ttmbPerMb\ttmbStatus\tsvTumorMutationalBurden\n', newLine: false, skip: 1)


/*
*
* Functions to create channels from TSV or directories containing BAM/CRAM
*
*/

//we read the pairs from tn_file
def parse_tn_file (tn_file,path,cram){
	    // FOR INPUT AS A TAB DELIMITED FILE
			def file_ext = cram ? '.crai':'.bai'
			//[sample t[.bam,cram] t[.bai,crai] n[.bam,.cram] n[.bai,.crai]]
    def tn_pairs=Channel.fromPath(tn_file)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map{row -> [ row.tumor_id,
               file(path + "/" + row.tumor),
               file(path + "/" + row.tumor+file_ext),
               file(path + "/" + row.normal),
               file(path + "/" + row.normal+file_ext)]}
      .ifEmpty{exit 1, "${tn_file} was empty - no tumor/normal supplied" }
	//we return the channel
  return tn_pairs
}

// print the calling parameter to the log and a log file
def print_params () {
  //software versions for v2.0
  def software_versions = ['hmftools-cobalt' : '1.11',
                          'hmftools-amber'   : '2.52',
                          'hmftools-purple' : '3.5' ]
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
}


//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        F\u001b[31;1mA\u001b[32;1mC\u001b[33;1mE\u001b[0mT\u001b[33;1ms\u001b[31;1m : Somatic\u001b[32;1m Copy\u001b[33;1m Number\u001b[31;1m Variant\u001b[33;1m caller\u001b[32;1m (${workflow.manifest.version})
        """
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}
