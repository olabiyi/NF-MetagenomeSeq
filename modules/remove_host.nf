#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
params.host_name  = null //"host"
params.host_url   = null //"http://www."
params.host_fasta = null // "/path/to/host_sequences.fasta"
params.host_db    = null // "path/to/host/database"
*/


/* Build kraken2 database in one of 3 modes
 1. User supplied host url
 2. User supplied custume sequence(s) and hostname 
 3. User supplied hostname
*/

process BUILD_HOSTDB {

    tag "Building kraken database..."
    label "kraken2"
    label "db_setup"

    input:
        tuple val(host_name), val(host_url) , path(host_fasta)

    output:
        path("kraken2-${host_name}-db/"), emit: krakendb_dir
        path("versions.txt"), emit: version

    script:
        """
        # Download and unpack database from URL
        if [ "${host_url}" != 'null' ]; then

             echo "Downloading and unpacking database from ${host_url}"
             wget -O ${host_name}.tar.gz --timeout=3600 --tries=0 --continue  ${host_url}

             mkdir kraken2-${host_name}-db/ && tar -zxvf ${host_name}.tar.gz -C kraken2-${host_name}-db/ && \\

            # Cleaning up
            [ -f  ${host_name}.tar.gz ] && rm -rf  ${host_name}.tar.gz

        # Build custome host reference
        elif [ "${host_fasta.name}" != 'empty.txt' ]; then

            echo "Attempting to build a custome ${host_name} reference database from ${host_fasta}"

            if [ ! -e "${host_fasta}" ]; then

                echo "${host_fasta} does not exist. please supply a valid fasta file and try again"
                exit 1

            fi

            # Install taxonomy
            kraken2-build --download-taxonomy --db kraken2-${host_name}-db/
            # Add sequence to your database's genomic library
            kraken2-build --add-to-library ${host_fasta} \\
                          --db kraken2-${host_name}-db/ --no-masking
            # Once your library is finalized, you need to build the database
            kraken2-build --build --db kraken2-${host_name}-db/

        # Build reference from named host e.g human
        elif [ "${host_name}" != 'null' ];then

            echo "Build kraken reference from ${host_name}"
            kraken2-build --download-library ${host_name} \\
                          --db kraken2-${host_name}-db/ \\
                          --threads ${task.cpus} --no-masking
            kraken2-build --download-taxonomy \\
                          --db kraken2-${host_name}-db/
            kraken2-build --build --db kraken2-${host_name}-db/ \\
                          --threads ${task.cpus}
            kraken2-build --clean --db kraken2-${host_name}-db/


        else

            echo "Input error! host_name, host_url and host_fasta are all set to null"
            echo "You must supply at least a valid host_name for database creation"
            exit 1

        fi

    VERSION=`echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//'`
    echo "kraken2 \${VERSION}"  > versions.txt
    """
}


// Remove host sequence using kranken2
process REMOVE_HOST {

    tag "Removing host reads"
    label "kraken2"

    input:
        val(host_suffix)
        each path(HOST_DB)
        tuple val(sample_id), path(reads), val(isPaired)

    output:
        tuple val(sample_id), path("*_${host_suffix}.fastq.gz"), val(isPaired), emit: reads
        tuple val(sample_id), path("*-kraken2-report_${host_suffix}.tsv"), emit: report
        path("versions.txt"), emit: version

    script:
    """
        if [ ${isPaired} == 'true' ]; then

        kraken2 --db ${HOST_DB} --gzip-compressed \\
            --threads ${task.cpus} \\
            --use-names --paired \\
            --output ${sample_id}-kraken2-output_${host_suffix}.txt \\
            --report ${sample_id}-kraken2-report_${host_suffix}.tsv \\
            --unclassified-out ${sample_id}_R#.fastq ${reads[0]} ${reads[1]}

        # Rename and gzip output files
        mv  ${sample_id}_R_1.fastq ${sample_id}_R1_${host_suffix}.fastq && \\
        gzip ${sample_id}_R1_${host_suffix}.fastq

        mv  ${sample_id}_R_2.fastq ${sample_id}_R2_${host_suffix}.fastq && \\
        gzip ${sample_id}_R2_${host_suffix}.fastq

    else

        # Single end
        kraken2 --db ${HOST_DB} --gzip-compressed \\
            --threads ${task.cpus} --use-names \\
            --output ${sample_id}-kraken2-output_${host_suffix}.txt \\
            --report ${sample_id}-kraken2-report_${host_suffix}.tsv \\
            --unclassified-out ${sample_id}.fastq ${reads[0]}

	# Rename and gzip output file	 
        mv ${sample_id}.fastq ${sample_id}_${host_suffix}.fastq && \\
        gzip ${sample_id}_${host_suffix}.fastq

    fi

    VERSION=`echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//'`
    echo kraken2 \${VERSION}  > versions.txt
    """
}




workflow remove_host {

    take:
       host_suffix // "HostRM", "HRrm"
       host_name
       host_url
       host_fasta
       host_db_dir
       reads_ch


    main:

       // Initialize software versions channel
       software_versions_ch = Channel.empty()

       // Use user supplied database path
       if ( host_db_dir ){

            host_db = Channel.fromPath(host_db_dir, checkIfExists: true)

       // Build database if path to existing database isn't supplied
       }else{

           fasta = file_path = host_fasta ?: file('empty.txt')
           input_ch = Channel.of([host_name, host_url , fasta])
           BUILD_HOSTDB(input_ch)
           host_db =  BUILD_HOSTDB.out.krakendb_dir
           BUILD_HOSTDB.out.version | mix(software_versions_ch) | set{software_versions_ch}
       }

       REMOVE_HOST(host_suffix, host_db, reads_ch)

       // Collect software versions
       REMOVE_HOST.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        clean_reads  = REMOVE_HOST.out.reads
        logs         =  REMOVE_HOST.out.report
        versions     = software_versions_ch
}
    
