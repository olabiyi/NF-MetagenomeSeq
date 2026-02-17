#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.paired = false
//params.max_mem = 100e9

/**************************************************************************************** 
**************************  Sequence assembly and summary *******************************
****************************************************************************************/

// This process handles running the assembly for each individual sample.
process ASSEMBLE {

    tag "Assembling ${sample_id}-s reads using megahit..."

    input:
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}_final.contigs.fa"), emit: contigs
        path("${sample_id}-assembly.log"), emit: log
        path("versions.txt"), emit: version
    script:
        """
        # Removing output directory if exists already but process still needs to be 
        # run (because there is no --force option to megahit i dont't think):        
        [ -d ${sample_id}-megahit-out/ ] && rm -rf ${sample_id}-megahit-out/

        if [ ${isPaired} == true ]; then
       
            BASENAME_FORWARD=`basename -s '.gz' ${reads[0]}`
            BASENAME_REVERSE=`basename -s '.gz' ${reads[1]}` 

            zcat  ${reads[0]} > \${BASENAME_FORWARD}
            zcat  ${reads[1]} > \${BASENAME_REVERSE}

            megahit -1 \${BASENAME_FORWARD} -2 \${BASENAME_REVERSE} \\
               -m ${params.max_mem} -t ${task.cpus} \\
               --min-contig-len 500 -o ${sample_id}-megahit-out > ${sample_id}-assembly.log 2>&1
         
        else

            BASENAME=`basename -s '.gz' ${reads[0]}`
            zcat ${reads[0]} > \${BASENAME}
            megahit -r \${BASENAME} -m ${params.max_mem} -t  ${task.cpus} \\
                --min-contig-len 500 -o ${sample_id}-megahit-out > ${sample_id}-assembly.log 2>&1
        fi
        
        mv ${sample_id}-megahit-out/final.contigs.fa ${sample_id}_final.contigs.fa
        megahit -v > versions.txt
        """
}


process FLYE {

    tag "Assembling ${sample_id}-s reads.."
    label "flye"
    label "assembly"

    input:
    tuple val(sample_id), path(reads), val(isPaired)

    output:
    tuple val(sample_id), path("${sample_id}-assembly.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}-assembly.log"), emit: log
    path("versions.txt"), emit: version

    script:
    """
    flye --meta \\
        --out-dir . \\
        --threads ${task.cpus} \\
        --nano-hq  ${reads}  || touch assembly.fasta 

    mv assembly.fasta ${sample_id}-assembly.fasta
    mv flye.log ${sample_id}-assembly.log

    VERSION=`flye --version`
    echo "flye \${VERSION}" > versions.txt
    """
}

// Assembly polishing with medaka
process POLISH_ASSEMBLY {

    tag "Polishing ${sample_id}-s assembly with medaka.."

    input:
    tuple val(sample_id), path(assembly), path(reads), val(isPaired)

    output:
    tuple val(sample_id), path("${sample_id}_polished.fasta"), emit: contigs
    tuple val(sample_id), path("${sample_id}-medaka.log"), emit: log
    path("versions.txt"), emit: version

    script:
    """
    # Check if contig assembly was successul before attempting to polish with medaka
    if [ -s ${assembly} ]; then

        medaka_consensus \\
            -t ${task.cpus} \\
            -i ${reads} \\
            -d ${assembly} \\
            -o .  >  ${sample_id}-medaka.log

    else 

        printf "${sample_id}\\tNo contig assembled, hence, contig polishing wasn't performed.\\n"
        touch consensus.fasta ${sample_id}-medaka.log
    
   fi
    mv consensus.fasta ${sample_id}_polished.fasta

    VERSION=`medaka --version 2>&1 | sed 's/medaka //g'`
    echo "medaka \${VERSION}" > versions.txt
    """
}

process SPADES {

    tag "Assembling ${sample_id}-s reads.."

    input:
        tuple val(sample_id), path(reads), val(isPaired) 
        val(type) // illumina or pacbio or nanopore

    output:
    tuple val(sample_id), path('*.scaffolds.fa'), optional:true, emit: scaffolds
    tuple val(sample_id), path('*warnings.log'), optional:true, emit: warnings
    tuple val(sample_id), path('*spades.log'), emit: log
    path("versions.txt"), emit: version

    script:
        def maxmem = task.memory.toGiga()
    """
    if [ ${type} ==  "illumina" ]; then

         if [ isPaired == 'true' ]; then

              INPUT=' -1 ${reads[0]} -2 ${reads[1]}'
         else

              INPUT=' -s ${reads[0]}'
         fi

    fi

 
    if [ ${type} ==  "pacbio" ]; then

       INPUT=' --pacbio ${reads[0]}' 
   
    fi


    if [ ${type} ==  "nanopore" ]; then

       INPUT=' --nanopore ${reads[0]}'

    fi    


    spades.py --meta --threads ${task.cpus} --memory ${maxmem} \${INPUT} -o .


    # Renaming output files    
    mv spades.log ${sample_id}-spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${sample_id}.scaffolds.fa
    fi

    if [ -f warnings.log ]; then
        mv warnings.log ${sample_id}-warnings.log
    fi

 
    VERSION=`spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p'`

    echo "spades \${VERSION}" > versions.txt

    """
}


process RENAME_HEADERS {

    tag "Renaming ${sample_id}-s assembly fasta file-s headers..."
    label "bit"
    label "assembly"

    input:
        tuple val(sample_id), path(assembly)
    output:
        tuple val(sample_id), path("${sample_id}-assembly.fasta"), emit: contigs
        path("versions.txt"), emit: version
        path("Failed-assemblies.tsv"), optional: true, emit: failed_assembly
    script:
        """
        bit-rename-fasta-headers -i ${assembly} \\
                                 -w c_${sample_id} \\
                                 -o ${sample_id}-assembly.fasta

        # Checking the assembly produced anything (megahit can run, produce 
        # the output fasta, but it will be empty if no contigs were assembled)
        if [ ! -s ${sample_id}-assembly.fasta ]; then
            printf "${sample_id}\\tNo contigs assembled\\n" > Failed-assemblies.tsv
        fi
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


// This process summarizes and reports general stats for all individual sample assemblies in one table.
process SUMMARIZE_ASSEMBLIES {

    tag "Generating a summary of all the assemblies..."
    label "bit"
    label "assembly"

    input:
        path(assemblies)      
    output:
        path("${params.additional_filename_prefix}assembly-summaries${params.assay_suffix}.tsv"), emit: summary
        path("versions.txt"), emit: version
    script:
        """
        bit-summarize-assembly \\
                 -o ${params.additional_filename_prefix}assembly-summaries${params.assay_suffix}.tsv \\
                 ${assemblies}
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


