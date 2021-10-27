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
      --ref                [file] fasta file of reference genome [hg38.fa], should be indexed [hg38.fa.fai]
      --ref_dict           [file] dict file for the reference genomep [hg38.dict]
    Optional arguments:
      --tumor_only         [flag] active tumor_only mode
      --bam                [flag] active bam mode [def:cram]
      --output_folder      [string] name of output folder
      --cpu                [Integer]  Number of CPUs[def:2]
      --mem 		           [Integer] Max memory [def:8Gb]
      --somatic_vcfs       [file] file containing list of somatic VCF variants with tumor_id and vcf_path, by default assumes that tumor_sample is the second listed sample of the VCF file.
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

//function that read the tumors to process from a tn_file
if(params.tn_file){
  def cram = params.bam ? false:true
 tn_pairs = parse_tn_file(params.tn_file,params.cohort_dir,cram)
 //we duplicate the tn_pairs channel
 tn_pairs.into { tn_pairs_cobalt; tn_pairs_amber}
}
//function to read the somatic vcfs

if(params.somatic_vcfs){
  svcfs=parse_vcf_file(params.somatic_vcfs)
}else{
  svcfs=channel.empty()
}


//chanel for reference genome
ref_fasta = Channel.value(file(params.ref)).ifEmpty{exit 1, "reference file not found: ${params.ref}"}
ref_dict = Channel.value(file(params.ref_dict)).ifEmpty{exit 1, "Dict reference file not found: ${params.ref_dict}"}
ref_fai = Channel.value(file(params.ref+'.fai')).ifEmpty{exit 1, "index file not found: ${params.ref}.fai"}


print_params()

//PATHS in the container for databases
// /hmftools/hg38

// /hmftools/hg38/DiploidRegions.38.bed
// /hmftools/hg38/GC_profile.1000bp.38.cnp
// /hmftools/hg38/GermlineHetPon.38.vcf
//run preprocesing with bcftools
process HQ_VCF{
  cpus '1'
  memory '1G'
  tag "hq_vcf"

  publishDir params.output_folder+'/VCF_highconf/', mode: 'copy'
  input:
  set val(tumor_id), file(vcf) from svcfs
  output:
  set val(tumor_id), file("${vcf.baseName}_highconf.vcf.gz") into hc_vcfs
  when:
  params.somatic_vcfs != 'null'

  script:
  if(params.debug == false){
  """
  #we assume that tumor WGS is the second sample in the VCF file
  # SNVs variants with a minimum read support of 5 are kept as high quality
  bcftools view -Oz -e 'FORMAT/AD[1:0]<5 | TYPE!="snps"' ${vcf} -o ${vcf.baseName}_filter.vcf.gz
  #we reheader the vcf file to match purple, cobalt, and amber sample names
  zcat ${vcf.baseName}_filter.vcf.gz | egrep "^#CHR" |\
  awk -v tid=${tumor_id} '{print \$(NF-1)" "tid"_N"; print \$NF" "tid"_T"}' > rename_sample.txt
  bcftools reheader -s rename_sample.txt -o ${vcf.baseName}_highconf.vcf.gz ${vcf.baseName}_filter.vcf.gz
  """
  }else{
  """
    touch ${vcf.baseName}_highconf.vcf.gz
  """
  }
}


process COBALT {

 cpus params.cpu
 memory params.mem+'G'


  publishDir params.output_folder+'/COBALT/', mode: 'copy'
  input:
  set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_cobalt
  file(ref) from ref_fasta
  file(fai) from ref_fai

  output:
  set val(tumor_id), path("${tumor_id}_COBALT") into cobalt
  script:
  if(params.debug == false){
     if(params.tumor_only){
       """
     COBALT  -Xms1g -Xmx${params.mem}g -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
              -ref_genome ${ref} -tumor_only -tumor_only_diploid_bed /hmftools/hg38/DiploidRegions.38.bed \\
     	        -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_COBALT -threads ${params.cpu}
       """
     }else{
       """
     COBALT  -Xms1g -Xmx${params.mem}g -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
              -ref_genome ${ref} -reference ${tumor_id}_N -reference_bam ${normal} \\
     	        -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_COBALT -threads ${params.cpu}
      """
     }
  }else{
    //debug code with options
    """
    mkdir ${tumor_id}_COBALT
    touch ${tumor_id}_COBALT/${tumor_id}.cobalt
    """
  }
}

process AMBER {

 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/AMBER/', mode: 'copy'
  input:
  set val(tumor_id), file(tumor), file(tumor_index), file(normal), file(normal_index) from tn_pairs_amber
  file(ref) from ref_fasta
  file(fai) from ref_fai
  output:
  set val(tumor_id), path("${tumor_id}_AMBER") into amber
  script:
  if(params.debug == false){
     if(params.tumor_only){
       """
      AMBER  -Xms1g -Xmx${params.mem}g  -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} -tumor_only \\
              -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads ${params.cpu}
      """
     }else{
       """
       AMBER   -Xms1g -Xmx${params.mem}g -loci /hmftools/hg38/GermlineHetPon.38.vcf -ref_genome ${ref} \\
               -reference ${tumor_id}_N -reference_bam ${normal}  \\
               -tumor  ${tumor_id}_T -tumor_bam ${tumor} -output_dir ${tumor_id}_AMBER -threads ${params.cpu}
        """
     }
  }else{
    """
    mkdir ${tumor_id}_AMBER
    touch ${tumor_id}_AMBER/${tumor_id}.amber
    """
  }
}

//we merge previous results from amber and cobalt
amber_cobalt=amber.join(cobalt, remainder: true)
//we ask if somatics where given and join them to the amber_cobalt struct
if(params.somatic_vcfs){
  amber_cobalt_vcf=amber_cobalt.join(hc_vcfs,remainder:true)
}else{
  //when no vcfs files are given
  vcf_dump = file('NO_VCF')
  amber_cobalt_vcf=amber_cobalt.spread([vcf_dump])
}

//dumping some channels
//amber_cobalt_vcf.into{print_vcfs;amber_cobalt_vcf2}
//print_vcfs.view()

//we have to merge the different inputs
process PURPLE {

 cpus params.cpu
 memory params.mem+'G'

  publishDir params.output_folder+'/PURPLE/', mode: 'copy'

  input:
  set val(tumor_id), path(amber_dir), path(cobalt_dir), file(hcvcf) from amber_cobalt_vcf
  file(ref) from ref_fasta
  file(fai) from ref_fai
  file(dict) from ref_dict
  output:
  set val(tumor_id), path("${tumor_id}_PURPLE") into purple
  file("${tumor_id}_T.purple.purity.sample.tsv") into stats_purple

  script:
  def include_vcf_purple = hcvcf.name != 'NO_VCF' ? "-somatic_vcf ${hcvcf}" : ''
  if(params.debug == false){
     if(params.tumor_only){
       """
       PURPLE  -Xms1g -Xmx${params.mem}g -tumor_only  -tumor ${tumor_id}_T \\
               -no_charts \\
               -output_dir ${tumor_id}_PURPLE \\
               -amber ${amber_dir} \\
               -cobalt ${cobalt_dir} \\
               -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
               -threads ${params.cpu} \\
               -ref_genome ${ref} ${include_vcf_purple}

        awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_T.purple.purity.sample.tsv

       """
     }else{
       """
        PURPLE  -Xms1g -Xmx${params.mem}g -reference ${tumor_id}_N  -tumor ${tumor_id}_T \\
               -no_charts \\
               -output_dir ${tumor_id}_PURPLE \\
               -amber ${amber_dir} \\
               -cobalt ${cobalt_dir} \\
               -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
               -threads ${params.cpu} \\
               -ref_genome ${ref} ${include_vcf_purple}

         awk -v tumor=${tumor_id} '{print tumor"\t"\$0}' ${tumor_id}_PURPLE/${tumor_id}_T.purple.purity.tsv > ${tumor_id}_T.purple.purity.sample.tsv
       """
     }
   }else{
     """
     echo PURPLE  -Xms1g -Xmx${params.mem}g -reference ${tumor_id}_N  -tumor ${tumor_id}_T \\
            -no_charts \\
            -output_dir ${tumor_id}_PURPLE \\
            -amber ${amber_dir} \\
            -cobalt ${cobalt_dir} \\
            -gc_profile /hmftools/hg38/GC_profile.1000bp.38.cnp \\
            -threads ${params.cpu} \\
            -ref_genome ${ref} ${include_vcf_purple}
      mkdir ${tumor_id}_PURPLE
      touch ${tumor_id}_T.purple.purity.sample.tsv
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

//we parse the somatics vcf files
def parse_vcf_file (vcf_file){
	    //[sample_id vcf_path]
    def svcfs=Channel.fromPath(vcf_file)
      .splitCsv(header: true, sep: '\t', strip: true)
      .map{row -> [ row.tumor_id,
               file(row.vcf_path)]}
      .ifEmpty{exit 1, "${svcfs} was empty - no vcf files supplied" }
	//we return the channel
  return svcfs
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
        PURPLE: Somatic CNV caller (${workflow.manifest.version})
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
