#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Read mapping to contig assembly using Bowtie2 ********************
****************************************************************************************/

// process to build minimap index
process MINIMAP_INDEX {

    tag "Building ${sample_id}-s index..."
    label "minimap2"

    input:
        tuple val(sample_id), path(assembly)
    output:
        tuple val(sample_id), path("${sample_id}.mmi"), emit: index
        path("versions.txt"), emit: version
    script:
        """
        minimap2 -ax splice -t ${task.cpus} -d ${sample_id}.mmi ${assembly}
        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}



// This process builds the bowtie2 index and runs the mapping for each sample
process LONG_MAPPING {

    tag "Mapping ${sample_id}-s reads to its assembly ${assembly}..."
    label "minimap2"
    label "mapping"
    

    input:
        tuple val(sample_id), path(assembly), path(reads)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

               minimap2 -ax map-ont -t ${task.cpus} ${assembly} ${reads} \\
                     > ${sample_id}.sam  2> ${sample_id}-mapping-info.txt
                        
     
            else

                touch ${sample_id}.sam
                echo "Mapping not performed for ${sample_id} because the assembly didn't produce anything."  > ${sample_id}-mapping-info.txt
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi

        VERSION=`minimap2 --version`
        echo "minimap2 \${VERSION}" > versions.txt
        """
}
