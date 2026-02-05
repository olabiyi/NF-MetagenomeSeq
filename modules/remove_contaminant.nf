#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SAMTOOLS_FILTER_FASTQ as  NANO_REMOVE_CONTAMINANT} from "./samtools.nf"

// max_mem = 100e9 // 100GB

process SPADES {

    tag "Assembling contaminants..."


    input:
        path(forward)
        path(reverse)
        val(isPaired)
        val(type) // illumina or pacbio or nanopore

    output:
        path('blank-scaffolds.fasta'), optional:true, emit: assembly
        path('blank-warnings.log'), optional:true, emit: warnings
        path('blank-assembly.log'), emit: log
        path("versions.txt"), emit: version

    script:
        def maxmem = task.memory.toGiga()
    """
    if [ ${type} ==  "illumina" ];then

         if [ ${isPaired} == 'true' ];then

              zcat ${forward} > merged_R1.fastq
              zcat ${reverse} > merged_R2.fastq

              INPUT=' -1 merged_R1.fastq -2 merged_R2.fastq'
         else
              
              zcat ${forward} > merged.fastq
              INPUT=' -s merged.fastq'
         fi

    fi

    if [ ${type} ==  "pacbio" ];then

       zcat ${forward} > merged.fastq

       INPUT=' --pacbio merged.fastq'

    fi


    if [ ${type} ==  "nanopore" ];then
 
       zcat ${forward} > merged.fastq

       INPUT=' --nanopore merged.fastq'

    fi


    spades.py --meta --threads ${task.cpus} --memory ${maxmem} \${INPUT} -o .


    # Renaming output files
    mv scaffolds.fasta blank-scaffolds.fasta
    mv spades.log blank-assembly.log

    [ -f warnings.log ] && mv warnings.log blank-warnings.log


    VERSION=`spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p'`
    echo "spades \${VERSION}" > versions.txt
    """
}
    


process FLYE {

    tag "Assembling contaminants.."
    label "flye"
 
    input:
       path(reads)

    output:
        path("blank-assembly.fasta"), emit: assembly
        path("blank-flye.log"), emit: log
        path("versions.txt"), emit: version

    script:
    """
    flye --meta \\
         --threads ${task.cpus} \\
         --out-dir . \\
         --nano-raw ${reads}

    mv assembly.fasta  blank-assembly.fasta
    mv flye.log blank-flye.log

    VERSION=`flye --version`
    echo "flye \${VERSION}" > versions.txt
    """
}


process BUILD_CONTAMINANT_DB {

    tag "Building a contaminant Database from blank samples"
    label "bowtie2"

    input:
        path(FASTA)

    output:
        path("blank-index/"), emit: index
        path("versions.txt"), emit: version

    script:
        """
        mkdir blank-index/
        bowtie2-build ${FASTA} blank-index/blanks

        bowtie2 --version  | \\
             head -n 1 | \\
             sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
        """
}


process BUILD_CONTAMINANT_INDEX {

    tag "Building a contaminant index from blank samples"
    label "minimap2"

    input:
        path(FASTA) // blank assembly
    output:
        path("blanks.mmi"), emit: index
        path("versions.txt"), emit: version
    script:
        """
        minimap2 -ax splice -t ${task.cpus} -d blanks.mmi ${FASTA}
        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}



process REMOVE_CONTAMINANT {

    tag "Removing reads mapping to the blanks assembly.."
    label "bowtie2"

    input:
        each path(BLANKS_DB)
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("*_decontam*.fastq.gz"), val(isPaired), emit: reads
        tuple val(sample_id), path("${sample_id}-mapping-info.txt"), emit: info
        path("versions.txt"), emit: version

    script:
        def input =  isPaired == "true" ?  
                 "-1 ${reads[0]} -2 ${reads[1]} --un-conc-gz" : 
                 "-U ${reads[0]} --un-gz"
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
        [ -z "\${INDEX}" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
        [ -z "\${INDEX}" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

        bowtie2 -p ${task.cpus} -x \${INDEX} --very-sensitive-local \\
               ${input} ${sample_id}_decontam.fastq.gz \\
               > ${sample_id}.sam 2> ${sample_id}-mapping-info.txt 

        # Rename Fastq Files
        if [ -f ${sample_id}_decontam.fastq.1.gz ]; then
            mv ${sample_id}_decontam.fastq.1.gz ${sample_id}_decontam_R1.fastq.gz
        fi

        if [ -f ${sample_id}_decontam.fastq.2.gz ]; then
            mv ${sample_id}_decontam.fastq.2.gz ${sample_id}_decontam_R2.fastq.gz
        fi


        # Delete sam file to free up memory
        rm -rf ${sample_id}.sam 

        bowtie2 --version  | \\
             head -n 1 | \\
             sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
    """
}


// This process builds the bowtie2 index and runs the mapping for each sample
process MAPPING_TO_CONTAMINANT {

    tag "mapping reads to the blanks assembly..."
    label "minimap2"
    

    input:
        each path(BLANKS_DB)
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
        INDEX=`find -L ./ -name "*.mmi"`
        minimap2 -t ${task.cpus} -ax splice \${INDEX} ${reads} \\
                     > ${sample_id}.sam  2> ${sample_id}-mapping-info.txt 
        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}



// A function to delete white spaces from an input string and covert it to lower case
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}


workflow illumina_remove_contaminants {


    take:
        file_ch
        filtered_ch

    main:

        file_ch.map{row -> if(deleteWS(row.NTC) == "true") [row.sample_id] }.set{blank_samples}

        isPaired  = filtered_ch.map{sample_id, reads, paired -> paired}.first()
         
         // Retain only blank samples
         filtered_ch.join(blank_samples).set{reads_ch}

         forward   = reads_ch.map{sample_id, reads, paired -> reads[0]}.toSortedList()
         reverse   = reads_ch.map{sample_id, reads, paired -> reads[1]}.toSortedList()

         
        SPADES(forward, reverse, isPaired, Channel.of("illumina"))
        BUILD_CONTAMINANT_DB(SPADES.out.assembly)
        REMOVE_CONTAMINANT(BUILD_CONTAMINANT_DB.out.index, filtered_ch)

        // Collect software versions
       software_versions_ch = Channel.empty()
       SPADES.out.version | mix(software_versions_ch) | set{software_versions_ch}
       BUILD_CONTAMINANT_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
       REMOVE_CONTAMINANT.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
       clean_reads = REMOVE_CONTAMINANT.out.reads
       versions = software_versions_ch

}



workflow nano_remove_contaminants {


    take:
        file_ch
        trimmed_ch

    main:

        // Get blank samples names from input file
        file_ch.map{row -> if(deleteWS(row.NTC) == "true") [row.sample_id] }.set{blank_samples}

        // Retain only blank samples
        trimmed_ch.join(blank_samples).set{reads_ch}
        // Get only the blank reads
        blank_reads   = reads_ch.map{sample_id, reads, paired -> reads}.toSortedList()

        // Assemble contaminants / blank reads
        FLYE(blank_reads)
        // Build contaminant index for mapping
        BUILD_CONTAMINANT_INDEX(FLYE.out.assembly)
        // Map all samples reads to contaminant assembly
        MAPPING_TO_CONTAMINANT(BUILD_CONTAMINANT_INDEX.out.index, trimmed_ch)
        
       /*
       - Sort and convert sam to bam
       - Index bam
       - Collect sample mapping stats
       - Remove unmapped reads/contaminants using samtools fastq
       */
       MAPPING_TO_CONTAMINANT.out.sam.map{ sample_id, sam, mapping_info ->
                         tuple([sample_id: sample_id, suffix: '_decontam', isPaired: 'false'],
                                sam, mapping_info)
                  }.set{mod_sam_ch}

       NANO_REMOVE_CONTAMINANT(mod_sam_ch)

       NANO_REMOVE_CONTAMINANT.out.log.map{ sample_id, stats, logs ->
                                   tuple(sample_id, [stats,logs])
                                   }.set{log_ch}

        // Collect software versions
       software_versions_ch = Channel.empty()
       FLYE.out.version | mix(software_versions_ch) | set{software_versions_ch}
       BUILD_CONTAMINANT_INDEX.out.version | mix(software_versions_ch) | set{software_versions_ch}
       MAPPING_TO_CONTAMINANT.out.version| mix(software_versions_ch) | set{software_versions_ch}
       NANO_REMOVE_CONTAMINANT.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
       clean_reads = NANO_REMOVE_CONTAMINANT.out.reads
       logs = log_ch
       versions = software_versions_ch

}
