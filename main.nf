#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '1.0'

params.help            = false
params.resume          = false

log.info """

╔╗ ┬┌─┐┌─┐┌─┐┬─┐┌─┐╔═╗╦═╗╔═╗  ┬┌┐┌┌┬┐┬─┐┌─┐┌─┐╔═╗╔═╗╔═╗ 
╠╩╗││ ││  │ │├┬┘├┤ ║  ╠╦╝║ ╦  ││││ ││├┬┘│ │├─┘╚═╗║╣ ║═╬╗
╚═╝┴└─┘└─┘└─┘┴└─└─┘╚═╝╩╚═╚═╝  ┴┘└┘─┴┘┴└─└─┘┴  ╚═╝╚═╝╚═╝╚
                                                                                       
====================================================
BIOCORE@CRG indropSEQ - N F  ~  version ${version}
====================================================
genome                        : ${params.genome}
annotation                    : ${params.annotation}
config                        : ${params.config}
barcode_list                  : ${params.barcode_list}
output (output folder)        : ${params.output}
"""

if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

annotationFile      = file(params.annotation) 
configFile        	= file(params.config[params.indrop_version]) 
barcodeFile        	= file(params.barcode_list[params.indrop_version]) 

outputfolder    = "${params.output}"
outputQC		= "${outputfolder}/QC"
outputMultiQC	= "${outputfolder}/multiQC"
outputMapping   = "${outputfolder}/Alignments";
filt_folder		= "${outputfolder}/Tagged_reads";
est_folder		= "${outputfolder}/Estimated_counts";
rep_folder		= "${outputfolder}/Reports";

if( !barcodeFile.exists() ) exit 1, "Missing barcode file: ${barcodeFile}"
if( !annotationFile.exists() ) exit 1, "Missing annotation file: ${annotationFile}"

/*
* if (params.strand == "yes") qualiOption = "strand-specific-forward"
* else if (params.strand != "no") qualiOption = "non-strand-specific"
*/

/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing
 * three elements: the pair ID, the first read-pair file and the second read-pair file
 */
Channel
    .from( params.input_files )
//    .ifEmpty { error "Please specify a list of input files using input_files parameter" }  
    .into { ch_input_files2; ch_input_files3 }

Channel
    .from( params.input_files )
//    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { reads_for_fastqc; fastq_files_for_size_est}    

ch_star_index = Channel.fromPath(params.star_index)
      .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }

Channel.fromPath(params.gtf)
  .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
  .into { ch_gtf_star; ch_gtf_featureCounts; }

/*
 * Step 0. Run FastQC on raw data
*/
process fastqc {
	publishDir outputQC

	tag { read }

    input:
    val read from reads_for_fastqc

    output:
   	file("*_fastqc.*") into raw_fastqc_files

    script:
	"""
    fastqc -t ${task.cpus} -q $read
	"""
}

    
/*
 * Step 1. Launch droptag for tagging your files
 */
process dropTag_inDrop_v2 {    
	publishDir filt_folder
	label 'indrop'

    when:
    params.indrop_version == 'v2'

    input:
    val inputs_names from ch_input_files2
    file configFile
    
    output:
    set params.tag, file("${params.tag}_tagged.fastq.gz") into tagged_files_for_alignment2
    file("${params.tag}_tagged.fastq.gz") into tagged_files_for_fastqc2
    
    script:
    """
	droptag -S -p ${task.cpus} -c ${configFile} ${inputs_names}
	zcat *.tagged.*.gz >> ${params.tag}_tagged.fastq
	gzip ${params.tag}_tagged.fastq
	rm 	*.fastq.gz.tagged.*.gz
    """
}   

/*
 * Step 1. Launch droptag for tagging your files
 */
process dropTag_inDrop_v3 {    
	publishDir filt_folder
	label 'indrop'

    when:
    params.indrop_version == 'v3'

    input:
    val inputs_names from ch_input_files3
    file configFile
    
    output:
    set params.tag, file("${params.tag}_tagged.fastq.gz") into tagged_files_for_alignment3
    file("${params.tag}_tagged.fastq.gz") into tagged_files_for_fastqc3
    
    script:
    """
	droptag -S -p ${task.cpus} -c ${configFile} ${inputs_names}
	zcat *.tagged.*.gz >> ${params.tag}_tagged.fastq
	gzip ${params.tag}_tagged.fastq
	rm 	*.fastq.gz.tagged.*.gz
    """
}   

tagged_files_for_alignment2.mix(tagged_files_for_alignment3).set{tagged_files_for_alignment}
tagged_files_for_fastqc2.mix(tagged_files_for_fastqc3).set{tagged_files_for_fastqc}

/*
 * Step 2. FastQC of your trimmed files
 */

process QCFiltReads {
	publishDir outputQC

	tag { filtered_read }

   	input:
    file(filtered_read) from tagged_files_for_fastqc.flatten()

    output:
   	file("*_fastqc.*") into trimmed_fastqc_files

    script:
    """
    fastqc -t ${task.cpus} -q $filtered_read
	"""
   }

/*
 * Step 3 extract read length of filtered reads?
*/

process getReadLength {   
	tag { fastq_file_for_size_est }

    input: 
    file(fastq_file_for_size_est) from fastq_files_for_size_est.first()
 
	output: 
	stdout into read_length

 	script:
 	"""
	if [ `echo ${fastq_file_for_size_est} | grep "gz"` ]; then cat="zcat"; else cat="cat"; fi
	\$cat ${fastq_file_for_size_est} | awk '{num++}{if (num%4==2){line++; sum+=length(\$0)} if (line==100) {printf "%.0f", sum/100; exit} }'
    """
} 

/*
 * Step 5. Mapping with STAR
 */

process star {
	label 'big_mem_cpus'
	publishDir outputMapping
	tag { pair_id }

	input:
    file index from ch_star_index.collect()
    file gtf from ch_gtf_star.collect()
	set pair_id, file(reads) from tagged_files_for_alignment	

	output:
	set pair_id, file('*.bam') into STARmappedTags_for_est

    script:
    """
    STAR --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --runDirPerm All_RWX \\
        --quantMode GeneCounts
    """

}


process dropEst {
	label 'indrop'
	publishDir est_folder
	tag { pair_id }

	input:
	set pair_id, file(tags) from STARmappedTags_for_est
	file annotationFile
    file configFile

	output:
	set pair_id, file ("*.rds")  into estimates_rds
	set pair_id, file ("*.tsv")  into estimates_tsv
	set pair_id, file ("*.mtx")  into estimates_mtx_for_plots
	set pair_id, file ("*.mtx")  into estimates_mtx

	script:		
    """
    ln -s ${barcodeFile} .
	dropest -m -w -g ${annotationFile} -c ${configFile} -o ${pair_id}.est ${tags} 
    """
	

}

/*
*
*/
process dropReport {
	label 'indrop'
	publishDir rep_folder
	tag { pair_id }
    errorStrategy 'ignore'

	input:
	set pair_id, file(estimate) from estimates_rds

	output:
	set pair_id, file ("${pair_id}_report.html")  into outreport

	script:		
    """
    dropReport.Rsc -o ${pair_id}_report.html ${estimate}
    """
}

workflow.onComplete {
    def subject = 'indropSeq execution'
    def recipient = "${params.email}"
    def attachment = "${outputMultiQC}/multiqc_report.html"

    ['mail', '-s', subject, '-a', attachment, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}




