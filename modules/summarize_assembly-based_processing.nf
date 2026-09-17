#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Summarize Assembly based metagenomics processing **************************
****************************************************************************************/

/*
 * ========================================================================================
 * PROCESS: GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE
 * ========================================================================================
 *
 * SUMMARY:
 *   Summarizing the results of assembly-based processing
 *
 * INPUTS:
 *   1. path: sample_IDs_file
 *      Cardinality: one
 *      Description: Input file: one column sample IDs file
 *
 *   2. path: MAGs_dir
 *      Cardinality: one
 *      Description: Input file: directory of MAGs
 *
 *   3. path: assemblies
 *      Cardinality: one
 *      Description: Input file: sample assemblies
 *
 *   4. path: genes_aa
 *      Cardinality: one
 *      Description: Input file: sample genes amino acid sequences
 *
 *   5. path: metabat_assembly_depth_files
 *      Cardinality: one
 *      Description: Input file: metabat assembly depth files
 *
 *   6. path: bins
 *      Cardinality: one
 *      Description: Input file: bins
 *
 *   7. path: bam_files
 *      Cardinality: one
 *      Description: Input file: sample bam files
 *
 * OUTPUTS:
 *   1. path: ${params.additional_filename_prefix}Assembly-based-processing-overview${params.assay_suffix}.tsv
 *
 * SOFTWARE & CONTAINERS:
 *   Container: [Defined in config/default.config]
 *   Conda: envs/bit.yaml
 *   Labels: bit
 *
 * RESOURCE REQUIREMENTS:
 *   - CPU cores: task.cpus
 *   - Memory: task.memory
 *
 * ========================================================================================
 */

process GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE {

    tag "Summarizing the results of assembly-based processing...."
    label "bit"

    input:
        path(sample_IDs_file)
        path(MAGs_dir)
        path(assemblies)
        path(genes_aa)
        path(metabat_assembly_depth_files)
        path(bins)
        path(bam_files)
    output:
        path("${params.additional_filename_prefix}Assembly-based-processing-overview${params.assay_suffix}.tsv")
    script:
        """
        mkdir assemblies_dir/ && mv *-assembly${params.assay_suffix}.fasta assemblies_dir/
        mkdir genes_dir/ && mv *-genes${params.assay_suffix}.faa genes_dir/ 
        mkdir mapping_dir/ && mv *-metabat-assembly-depth${params.assay_suffix}.tsv *.bam  mapping_dir/

        mkdir bins_dir/
        if compgen -G *-bin*.fasta > /dev/null; then
            mv  *-bin*.fasta  bins_dir/
        fi
        bash generate-assembly-based-overview-table.sh \\
                ${sample_IDs_file} \\
                assemblies_dir/ \\
                genes_dir/ \\
                mapping_dir/ \\
                bins_dir/ \\
                ${MAGs_dir}/ \\
                ${params.additional_filename_prefix}Assembly-based-processing-overview${params.assay_suffix}.tsv
        """
}

