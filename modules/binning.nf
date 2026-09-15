#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Assembly binning *************************************************
****************************************************************************************/

/*
 * ========================================================================================
 * PROCESS: METABAT_BINNING
 * ========================================================================================
 *
 * SUMMARY:
 *   Binning sample contigs with metabat2
 *
 * INPUTS:
 *   1. tuple: tuple val(sample_id), path(assembly), path(bam)
 *      Cardinality: one
 *      Description: Tuple input combining multiple channel elements
 *                 - sample_id: string specifying the input sample name
 *                 - assembly: path to sample assembly/contigs
 *                 - bam: path to sample bam file
 *
 * OUTPUTS:
 *   1. tuple: tuple val(sample_id), path("${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv") (emit: depth)
 *
 *   2. tuple: tuple val(sample_id), path("${sample_id}-bin*"), optional: true (emit: bins) [OPTIONAL]
 *
 *   3. path: versions.txt (emit: version)
 *
 * SOFTWARE & CONTAINERS:
 *   Primary Tool: MetaBAT2
 *   Container: [Defined in config/default.config]
 *   Conda: envs/metabat.yaml
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

// This process runs metabat2 for binning contigs.
process METABAT_BINNING {

    tag "Binning ${sample_id}-s contigs with metabat2..."

    input:
        tuple val(sample_id), path(assembly), path(bam) 
    output:
        tuple val(sample_id), path("${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv"), emit: depth
        tuple val(sample_id), path("${sample_id}-bin*"), emit: bins, optional: true
        path("versions.txt"), emit: version
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${assembly} ]; then

            jgi_summarize_bam_contig_depths \\
                    --outputDepth ${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv \\
                    --percentIdentity 97 \\
                    --minContigLength 1000 \\
                    --minContigDepth 1.0  \\
                    --referenceFasta ${assembly} ${bam} 

            # only running if there are contigs with coverage 
            # information in the coverage file we just generated
            if [ `wc -l ${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv | sed 's/^ *//' | cut -f 1 -d " "` -gt 1 ]; then 

                metabat2  \\
                    --inFile ${assembly} \\
                    --outFile ${sample_id}-bin \\
                    --abdFile ${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv \\
                    -t ${task.cpus}

            else

                printf "\\n\\nThere was no coverage info generated in ${sample_id}-metabat-assembly-depth.tsv, so no binning with metabat was performed.\\n\\n" 

            fi

            # changing extensions from .fa to .fasta to match nt fasta extension elsewhere in GeneLab
            find . -name '${sample_id}*.fa' > ${sample_id}-bin-files.tmp
            
            if [ -s ${sample_id}-bin-files.tmp ]; then
                paste -d " " <( sed 's/^/mv /' ${sample_id}-bin-files.tmp ) \\
                             <( sed 's/.fa/.fasta/' ${sample_id}-bin-files.tmp ) \\
                             > ${sample_id}-rename.tmp
                bash ${sample_id}-rename.tmp
            fi

            rm -rf ${sample_id}-bin-files.tmp ${sample_id}-rename.tmp

        else

            touch ${sample_id}-metabat-assembly-depth${params.assay_suffix}.tsv
            printf "Binning not performed because the assembly didn't produce anything.\\n" 
        fi
        echo metabat2 \$(metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/') > versions.txt
        """
}



workflow binning {

    take:
        assembly_ch
        read_mapping_ch


    main:
        binning_ch = METABAT_BINNING(assembly_ch.join(read_mapping_ch))


    emit:
    binning_results = binning_ch.out.bins

}

