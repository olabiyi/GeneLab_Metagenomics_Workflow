#!/usr/bin/env python
"""
Metagenomics Assay Table Updater for the NASA GeneLab Data Processing Pipeline

This script processes and updates Metagenomics assay tables from ISA.zip files within a GLDS dataset.
It adds or updates various parameter columns required for the GeneLab Metagenomics processing pipeline,
including paths to processed files.

Usage:
    python update_assay_table.py --assay_suffix <suffix> --glds_accession <GLDS-XXX> [--mode <mode>]

Parameters:
    --outdir: Directory containing the Metadata folder with ISA.zip
    --assay_suffix: Suffix to append to output filenames (e.g., "_GLmetagenomics")
    --glds_accession: GLDS accession number (e.g., "GLDS-123")
    --technology: Technology used for sequencing, either "illumina" or "nanopore"
    --sample-type: Sample type sequenced, either "low_biomass" or "standard"
    --host-removed: Indicates if host sequences were removed from the dataset
    --single-ended: Add this flag if illumina data are single-end sequencing.
 
 The script automatically:
1. Extracts and locates the assay table from the ISA.zip file
2. Loads the runsheet if available for additional metadata
3. Detects if data is paired-end or single-end
4. Updates the assay table with appropriate file paths for all processing outputs
5. Saves the updated assay table using the original filename
"""

import sys
import argparse
from typing import Any
import zipfile
import json
import re
from pathlib import Path
import pandas as pd
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        prog="update_assay_table",
        description="Update Metagenomics assay table from ISA.zip with processed data file information.")
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-s', '--assay_suffix', required=True, action='store', 
        choices=["_GLmetagenomics", "_GLlblMetag", "_GLlbsMetag", ""],
        help="Specifies the assay suffix for the dataset (empty string if no suffix).")
    required.add_argument('-r', '--runsheet', required=True, help="Runsheet used to process the data")
    required.add_argument('-g', "--glds_accession", required=True, 
                          help="GLDS accession number (e.g. GLDS-123)")
    isa_file_input = required.add_mutually_exclusive_group(required=True)
    isa_file_input.add_argument('-i', '--isa_zip', action='store', default='',
                                help='Appropriate ISA file for the dataset (a zip archive, providing this instead of '
                                     'an assay table directly will attempt to extract the correct assay table given '
                                     'the provided assay and technology types.)')
    isa_file_input.add_argument('-a', '--assay_table', action='store', default='',
                                help='Assay table for the dataset provided directly instead of extracting from an ISA '
                                     'zip file')
    parser.add_argument("--sample-type", default="standard", action="store", choices=["standard", "low_biomass"], 
                        help="Sample type for the dataset ('standard' or 'low_biomass').")
    parser.add_argument("--technology", default="illumina", action="store", choices=["illumina", "nanopore"], 
                        help="Sequencing technology used for the dataset ('illumina' or 'nanopore').")
    parser.add_argument("--output-prefix", default="", type=str, action="store", help="Dataset file prefix.")
    parser.add_argument("--validation-json", type=Path, action="store", help="Path to validation manifest json")
    parser.add_argument("--read-stats", type=Path, action="store", help="Path to read stats (if available)")
    return parser.parse_args()


tty_colors = {
    "green": "\033[0;32m%s\033[0m",
    "yellow": "\033[0;33m%s\033[0m",
    "red": "\033[0;31m%s\033[0m",
}


def color_text(text:str, color:str="green"):
    """
    Colors text for output in terminal

    Args:
        text (str): input text
        color (str): a valid tty color ('red', 'yellow', or 'green')

    Returns:
        str: colored text

    """
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def report_failure_and_exit(message:str, color:str="red"):
    """
    Reports a failure and exits with status '1'.

    Args:
        message (str): Error message to report.
        color (str): Color in which to render the error message, default = 'red'
    """
    print("")
    print(color_text(f"Error: {message}", color))
    sys.exit("\nAssay table update failed.\n")


def report_warning(message: str, color:str="yellow"):
    """
    Reports a warning message.

    Args:
        message (str): Warning message to report
        color (str): Color in which to render the warning message, default = 'yellow'
    """
    print("")
    print(color_text(f"Warning: {message}", color))


def report_info(message:str, color:str="green"):
    """
    Reports an info message.

    Args:
        message (str): Info message to report
        color (str): Color in which to render the info message, default = 'green'
    """
    print("")
    print(color_text(f"Info: {message}", color))


def load_runsheet(runsheet_file: Path):
    """
    Load the runsheet as a pandas.DataFrame.

    Args:
        runsheet_file (Path): a file containing the assay runsheet used to generate the processed data

    Returns:
        pandas.DataFrame: sample information from the runsheet
    """
    try:
        runsheet_df = pd.read_csv(runsheet_file)
        print(f"Runsheet has {len(runsheet_df)} rows and {len(runsheet_df.columns)} columns")
        return runsheet_df
    except Exception as e:
        report_warning(f"Cannot read runsheet, proceeding without it: {e}")
        return None


def load_validation_json(validation_json_file: Path):
    """
    Load the validation manifest json file.

    Args:
        validation_json_file (Path): a file containing the validation manifest json

    Returns:
        dict: validation manifest data
    """
    try:
        with open(validation_json_file, "r") as f:
            validation_data = json.load(f)
        return dict(validation_data)
    except Exception as e:
        report_warning(f"Cannot read validation manifest json, proceeding without it: {e}")
        return dict()


def get_verified_folders(validation_data: dict[Any, Any], common_root: Path):
    """
    Extracts the list of validation folders from the validation manifest data.

    Args:
        validation_data (dict): validation manifest data
    Returns: 
        dict: dictionary of output folders present in the validation manifest
    """
    folders = defaultdict(dict)
    if "verified" in validation_data:
        for file in validation_data["verified"]:
            p = Path(file).relative_to(common_root)
            pieces = list(p.parts)
            if pieces[len(pieces)-1] == "Fastq":
                pieces.remove('Fastq')
            subfolder = "-".join(pieces[1:len(pieces)-1])
            filename = p.name

            if subfolder in folders[pieces[0]]:
                folders[pieces[0]][subfolder].append(filename)
            else:
                folders[pieces[0]][subfolder] = [filename]
    
    return folders


def get_runsheet_sample_name_map(runsheet_df, assay_sample_names):
    """
    Generates a mapping of sample names in the assay table to the samplenames in the runsheet

    Args:
        runsheet_df (pandas.DataFrame): runsheet sample information
        assay_sample_names (list): sample names from the assay table

    Returns:
        dict: sample name mapping
    """
    sample_name_map = {}
    if 'Sample Name' in runsheet_df.columns:
        # Check for 'Original Sample Name' column to map between assay table and runsheet
        if 'Original Sample Name' in runsheet_df.columns:
            for _, row in runsheet_df.iterrows():
                orig_name = row['Original Sample Name']
                rs_name = row['Sample Name']
                if orig_name in assay_sample_names:
                    sample_name_map[orig_name] = rs_name
    return sample_name_map


def is_paired_end_data(runsheet_df):
    """
    Determine if this is paired-end data based on information the runsheet.

    Args:
        runsheet_df (pandas.DataFrame): runsheet sample information

    Returns:
        bool: True if paired-end, False if single-end
    """
    if runsheet_df is None:
        report_warning("Warning: No runsheet provided, assuming single-end data")
        return False

    # Check if there's a paired_end column
    if "paired" in runsheet_df.columns:
        paired_end_values = runsheet_df["paired"].unique()
        if len(paired_end_values) == 1:
            # If all values are the same, use that
            value = paired_end_values[0]
            # Handle different types of values (string or boolean)
            if isinstance(value, bool):
                return value
            elif isinstance(value, str):
                return value.lower() == "true"
            else:
                # Try to convert to string if it's not a boolean or string
                return str(value).lower() == "true"

    # Check based on R1/R2 file presence in first row
    if any(val for val in runsheet_df.iloc[0].values if "_R2_" in str(val)):
        print("Detected paired-end data based on _R2 files in runsheet")
        return True
    else:
        print("Assuming single-end data (no _R2 files in runsheet)")
        return False


def get_assay_table_from_isa(isa_file, assay, technology):
    """
    tries to find an assay table in an ISA zip file that matches the type expected for the provided assay

    Args:
        isa_file (PathLike[str]): path to ISA zip file
        assay (str): valid OSDR methylseq assay type ('MethylSeq' or 'RNAMethylSeq')
        technology (str): valid OSDR methylseq technology type

    Returns:
        pandas.DataFrame: assay table from extracted from ISA zip
    """

    zip_file = zipfile.ZipFile(isa_file)
    isa_files = zip_file.namelist()

    valid_measurement = "Metagenomic sequencing"

    # Parse investigation file to build STUDY ASSAYS table
    study_assays_table = {}
    study_assays_section = False

    with zip_file.open("i_Investigation.txt", "r") as f:
        for line in f:
            line = line.decode("utf-8").strip()

            # Track STUDY ASSAYS section
            if line == "STUDY ASSAYS":
                study_assays_section = True
                continue
            elif study_assays_section and not line:
                study_assays_section = False
                continue

            # Extract data from section
            if study_assays_section and line:
                parts = line.split("\t")
                if parts and parts[0]:
                    key = parts[0]
                    values = [v.strip() for v in parts[1:] if v.strip()]
                    study_assays_table[key] = values
    # Check if we have all required keys
    required_keys = [
        "Study Assay Measurement Type",
        "Study Assay Technology Type",
        "Study Assay File Name",
    ]
    if not all(key in study_assays_table for key in required_keys):
        report_failure_and_exit("Missing required keys in STUDY ASSAYS section")

    # Get the values from the table
    measurement_types = study_assays_table["Study Assay Measurement Type"]
    technology_types = study_assays_table["Study Assay Technology Type"]
    file_names = study_assays_table["Study Assay File Name"]

    # Ensure all lists have equal length
    if not (len(measurement_types) == len(technology_types) == len(file_names)):
        report_failure_and_exit(
            "Measurement types, technology types, and file names have different lengths"
        )

    # Find matching assay file
    matched_file = ""
    for i in range(len(measurement_types)):
        if (
            measurement_types[i].lower() == valid_measurement.lower()
            and technology_types[i].lower() == technology.lower()
        ):
            matched_file = file_names[i]
            break

    if not matched_file:
        report_failure_and_exit(
            f"No assay file matched for {assay}. "
            f"Measurement types: {measurement_types}, Technology types: {technology_types}"
        )
    elif matched_file not in isa_files:
        # Load the matched assay file
        report_failure_and_exit(
            f"Matched assay file doesn't exist in ISA zip: {matched_file}"
        )
    else:
        return pd.read_csv(zip_file.open(matched_file), sep="\t"), matched_file

    return pd.DataFrame(), matched_file


def add_data_column_to_dataframe(df: pd.DataFrame, column: pd.Series, column_name: str, unit: str):
    if column_name not in df.columns:
        print(f"Adding column: '{column_name}' with unit: '{unit}'")
        df[column_name] = column
        df["Unit"] = unit
    else:
        print(f"Updating column: '{column_name}' from '{df[column_name]}' to '{column}'")
        df[column_name] = column
    return df


def add_read_counts_from_read_stats_file(df: pd.DataFrame, read_stats_file: Path):

    print("Adding read stats")
    raw = pd.Series()
    hrrm_perc = pd.Series()
    decontam_perc = pd.Series()
    read_stats = pd.read_csv(Path(read_stats_file), sep="\t", index_col=0)
    raw = read_stats["Raw"].astype(int)
    hrrm_perc = read_stats["Percent_human_reads_removed"].round(2)
    decontam_perc = read_stats["Percent_blank_reads_removed"].round(2)

    if not raw.empty:
        print("Adding raw read depth")
        add_data_column_to_dataframe(df, raw, 
                                     "Parameter Value[Read Depth]", "reads")
    if not hrrm_perc.empty:
        print("Adding human reads removed")
        add_data_column_to_dataframe(df, hrrm_perc, 
                                     "Parameter Value[Human Reads Removed]", "percent")
    if not decontam_perc.empty:
        print("Adding contaminant reads removed")
        add_data_column_to_dataframe(df, decontam_perc, 
                                     "Parameter Value[Contaminant Reads Removed]", "percent")
    return df


def add_parameter_column(df: pd.DataFrame, column_name: str, value: str, glds_prefix:str=""):
    """
    Add a parameter column to the dataframe if it doesn't exist already. Use the 
    same value for all rows in the table.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each value (only for 
                           values that are filenames)
        column_name (str): parameter column name to add (e.g., "Parameter Value[Entry]")
        value (str): value to set for all rows in the table

    Returns:
        pandas.DataFrame: update assay table
    """
    # Apply prefix to value if provided
    if glds_prefix and isinstance(value, str):
        # Check if value already has the prefix
        if not value.startswith(glds_prefix):
            prefixed_value = f"{glds_prefix}{value}"
        else:
            prefixed_value = value
    else:
        prefixed_value = value

    if column_name not in df.columns:
        print(f"Adding new column: {column_name}")
        df[column_name] = prefixed_value
    else:
        print(f"Column {column_name} already exists, updating values")
        df[column_name] = prefixed_value

    return df


def add_fastq_data_column(df, osdr_prefix, file_prefix, assay_suffix, assay_sample_names,
                          verified_fastq_files, fastq_type="filtered", is_paired_end=True,):
    """
    Add the Trimmed Sequence Data column and the optional 
    Merged Sequence Data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to 
                                sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)
        fastq_type (str): 'filtered' or 'HRrm' or 'HostRm'

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = ""
    if fastq_type == "filtered":
        column_name = "Parameter Value[Filtered Sequence Data]"
    elif fastq_type == "HRrm":
        column_name = "Parameter Value[HR-removed Sequence Data]"
    elif fastq_type == "HostRm":
        column_name = "Parameter Value[Host-removed Sequence Data]"
    elif fastq_type == "decontam":
        column_name = "Parameter Value[Decontaminated Sequence Data]"
    else:
        report_failure_and_exit(
            f"Unrecognized fastq type: '{fastq_type}'. Must be one of ['filtered', 'HRrm', 'HostRm', 'decontam']"
        )

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        if is_paired_end:
            # For paired-end data, create entries with both R1 and R2 files, comma-separated without spaces
            # check if file exists
            fastq_file_r1 = f"{file_prefix}{sample}{assay_suffix}_R1_{fastq_type}{assay_suffix}.fastq.gz"
            fastq_file_r2 = f"{file_prefix}{sample}{assay_suffix}_R2_{fastq_type}{assay_suffix}.fastq.gz"
            if fastq_file_r1 in verified_fastq_files and fastq_file_r2 in verified_fastq_files: 
                values.append(f"{osdr_prefix}{fastq_file_r1},{osdr_prefix}{fastq_file_r2}")
            else:
                report_warning(f"Fastq files for sample {sample} not found: \n\t{fastq_file_r1}\n\t{fastq_file_r2}")
                values.append("")
        else:
            # For single-end data
            fastq_file = f"{file_prefix}{sample}{assay_suffix}_{fastq_type}{assay_suffix}.fastq.gz"
            if fastq_file in verified_fastq_files:
                values.append(f"{osdr_prefix}{fastq_file}")
            else:
                report_warning(f"Fastq file for sample {sample} not found: \n\t{fastq_file}")
                values.append("")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding new column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_assemblies(df, assay_sample_names: list, failed_assemblies: list, file_prefix: str = "", 
                   osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/assemblies]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = [f"{osdr_prefix}{file_prefix}assembly-summaries{assay_suffix}.tsv"]
        if sample not in failed_assemblies:
            files.append(f"{osdr_prefix}{file_prefix}{sample}-assembly{assay_suffix}.fasta")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_predicted_genes(df, assay_sample_names: list, failed_assemblies: list, file_prefix: str = "", 
                        osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/predicted-genes]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = []
        if sample not in failed_assemblies:
            files.append(f"{osdr_prefix}{file_prefix}{sample}-genes{assay_suffix}.fasta")
            files.append(f"{osdr_prefix}{file_prefix}{sample}-genes{assay_suffix}.faa")
            files.append(f"{osdr_prefix}{file_prefix}{sample}-genes{assay_suffix}.gff")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_read_mapping(df, assay_sample_names: list, failed_assemblies: list,
                     file_prefix: str = "", osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/read-mapping]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = []
        if sample not in failed_assemblies:
            files.append(f"{osdr_prefix}{file_prefix}{sample}{assay_suffix}.bam")
            files.append(f"{osdr_prefix}{file_prefix}{sample}-mapping-info{assay_suffix}.txt")
            files.append(f"{osdr_prefix}{file_prefix}{sample}-metabat-assembly-depth{assay_suffix}.tsv")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_annotation_and_tax(df, assay_sample_names: list, failed_assemblies: list, file_prefix: str = "", 
                           osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/annotations-and-taxonomy]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = []
        if sample not in failed_assemblies:
            files.append(f"{osdr_prefix}{file_prefix}{sample}-contig-coverage-and-tax{assay_suffix}.bam")
            files.append(f"{osdr_prefix}{file_prefix}{sample}-gene-coverage-annotation-and-tax{assay_suffix}.txt")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_bins(df, assay_sample_names: list, failed_assemblies: list, bins_files: list,
             file_prefix: str = "", osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/bins]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = [f"{osdr_prefix}{file_prefix}bins-overview{assay_suffix}.tsv"]
        if sample not in failed_assemblies:
            sample_file = f"{file_prefix}{sample}-bins{assay_suffix}.zip"
            if sample_file in bins_files:
                files.append(f"{osdr_prefix}{sample_file}")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_MAGs(df, assay_sample_names: list, failed_assemblies: list, mags_files: list, 
             file_prefix: str = "", osdr_prefix: str = "", assay_suffix: str = ""):
    """
    Add the assemblies column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Assembly-based Processing/MAGs]"

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        files = [f"{osdr_prefix}{file_prefix}MAGs-overview{assay_suffix}.tsv"]
        if sample not in failed_assemblies:
            sample_file = f"{file_prefix}{sample}-MAGs{assay_suffix}.zip"
            if sample_file in mags_files:
                files.append(f"{osdr_prefix}{sample_file}")
        values.append(",".join(files))
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_combined_outputs(df: pd.DataFrame, verified_assembly_based: dict[str, list], osdr_prefix: str = ""):
    # Contig-level taxonomy
    column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Contig-level]"
    counts_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Contig-level/Counts Tables]"
    plots_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Contig-level/Heatmaps]"

    table_dir = "combined-outputs-Contig-level"
    table_files = []
    counts_files = []
    plot_files = []
    if table_dir in verified_assembly_based:
        for file in verified_assembly_based[table_dir]:
            if file.endswith(".tsv"):
                if 'filtered' in file or 'decontam' in file:
                    counts_files.append(f"{osdr_prefix}{file}")
                else:
                    table_files.append(f"{osdr_prefix}{file}")
            elif file.endswith(".png"):
                plot_files.append(f"{osdr_prefix}{file}")
    else:
        print("Assembly-based processing: combined-outputs/Contig-level folder does not exist.")

    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=",".join(table_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=counts_column_name, value=",".join(counts_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=plots_column_name, value=",".join(plot_files))

    # Gene-level taxonomy
    column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/Taxonomy]"
    counts_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/Taxonomy/Counts Tables]"
    plots_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/Taxonomy/Heatmaps]"

    table_dir = "combined-outputs-Gene-level-Taxonomy"
    table_files = []
    counts_files = []
    plot_files = []
    if table_dir in verified_assembly_based:
        for file in verified_assembly_based[table_dir]:
            if file.endswith(".tsv"):
                if 'filtered' in file or 'decontam' in file:
                    counts_files.append(f"{osdr_prefix}{file}")
                else:
                    table_files.append(f"{osdr_prefix}{file}")
            elif file.endswith(".png"):
                plot_files.append(f"{osdr_prefix}{file}")
    else:
        print("Assembly-based processing: combined-outputs/Gene-level/Taxonomy folder does not exist.")

    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=",".join(table_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=counts_column_name, value=",".join(counts_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=plots_column_name, value=",".join(plot_files))

    # Gene-level KO
    column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/KO]"
    counts_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/KO/Counts Tables]"
    plots_column_name = "Parameter Value[Assembly-based Processing/combined-outputs/Gene-level/KO/Heatmaps]"

    table_dir = "combined-outputs-Gene-level-KO"
    table_files = []
    counts_files = []
    plot_files = []
    if table_dir in verified_assembly_based:
        for file in verified_assembly_based[table_dir]:
            if file.endswith(".tsv"):
                if 'filtered' in file or 'decontam' in file:
                    counts_files.append(f"{osdr_prefix}{file}")
                else:
                    table_files.append(f"{osdr_prefix}{file}")
            elif file.endswith(".png"):
                plot_files.append(f"{osdr_prefix}{file}")
    else:
        print("Assembly-based processing: combined-outputs/Gene-level/KO folder does not exist.")

    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=",".join(table_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=counts_column_name, value=",".join(counts_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=plots_column_name, value=",".join(plot_files))

    return df


def add_nanoplot_data_column(df, osdr_prefix, file_prefix, assay_suffix, 
                             assay_sample_names, seq_type, verified_nanoplot_files):
    """
    Add the NanoPlot report data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
    """
    # Title Case column name

    column_name = ""
    if seq_type == "filtered":
        column_name = "Parameter Value[Filtered Sequence Data/NanoPlots]"
    elif seq_type == "raw":
        column_name = "Parameter Value[Merged Sequence Data/NanoPlots]"
    elif seq_type == "HRrm":
        column_name = "Parameter Value[HR-removed Sequence Data/NanoPlots]"
    elif seq_type == "decontam":
        column_name = "Parameter Value[Decontaminated Sequence Data/NanoPlots]"
    elif seq_type == "trimmed":
        column_name = "Parameter Value[Trimmed Sequence Data/NanoPlots]"
    elif seq_type == "HostRm":
        column_name = "Parameter Value[Host-removed Sequence Data/NanoPlots]"

    else:
        report_failure_and_exit(
            f"Unrecognized sequence type: '{seq_type}'. Must be one of ['raw', 'trimmed', 'filtered', 'decontam', 'HRrm', 'kraken2', 'HostRm']"
        )

    # Generate file paths using the appropriate sample names
    values = []
    for sample in assay_sample_names:
        filename = (f"{file_prefix}{sample}_{seq_type}_NanoPlot-report{assay_suffix}.html")
        if filename in verified_nanoplot_files:
            values.append(f"{osdr_prefix}{filename}")
        else:
            report_warning(f"Nanoplot file: {filename} does not exist")
            values.append("")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_multiqc_reports_column(df, glds_prefix, assay_suffix, multiqc_type):
    """
    Add the raw/trimmed/aligned sequence MultiQC reports data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        multiqc_type (str): specifies if the processing step for which these multiqc reports were generated

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = ""
    if multiqc_type == "filtered":
        column_name = "Parameter Value[Filtered Sequence Data/MultiQC Reports]"
    elif multiqc_type == "raw":
        column_name = "Parameter Value[Merged Sequence Data/MultiQC Reports]"
    elif multiqc_type == "HRrm":
        column_name = "Parameter Value[HR-removed Sequence Data/MultiQC Reports]"
    elif multiqc_type == "decontam":
        column_name = "Parameter Value[Decontaminated Sequence Data/MultiQC Reports]"
    elif multiqc_type == "trimmed":
        column_name = "Parameter Value[Trimmed Sequence Data/MultiQC Reports]"
    elif multiqc_type == "HostRm":
        column_name = "Parameter Value[Host-removed Sequence Data/MultiQC Reports]"
    elif multiqc_type == "kraken2":
        column_name = (
            "Parameter Value[Read-based Processing/Kraken2 Outputs/MultiQC Reports]"
        )

    else:
        report_failure_and_exit(
            f"Unrecognized multiqc type: '{multiqc_type}'. Must be one of ['raw', 'trimmed', 'filtered', 'decontam', 'HRrm', 'kraken2', 'HostRm']"
        )

    # Create the alignment multiqc report filename - same for all samples
    multiqc_report = (
        f"{glds_prefix}{multiqc_type}_multiqc{assay_suffix}.html,"
        f"{glds_prefix}{multiqc_type}_multiqc{assay_suffix}_data.zip"
    )

    # Add the column to the dataframe with the same value for all rows
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = multiqc_report
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = multiqc_report

    return df


def add_read_based_taxonomy_plots(df: pd.DataFrame, verified_read_based: dict[str, list], 
                                  output_type: str, osdr_prefix: str = ""):
    column_name = (
        f"Parameter Value[Read-based Processing/{output_type} Outputs/Barplots]"
    )

    plot_dir = f"{output_type}_Outputs-Barplots"
    plot_files = []
    if plot_dir in verified_read_based:
        for plotname in verified_read_based[plot_dir]:
            plot_files.append(f"{osdr_prefix}{plotname}")
    else:
        print(f"Read-based processing: barplots folder for {output_type} does not exist.")

    combined_files = ",".join(plot_files)

    # add column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=combined_files)
    return df


def add_read_based_taxonomy_counts(df: pd.DataFrame, verified_read_based: dict[str, list], 
                                   output_type: str, osdr_prefix: str = ""):
    column_name = (
        f"Parameter Value[Read-based Processing/{output_type} Outputs/Counts Tables]"
    )

    table_dir = f"{output_type}_Outputs-Count_Tables"
    table_files = []
    if table_dir in verified_read_based:
        for tablename in verified_read_based[table_dir]:
            table_files.append(f"{osdr_prefix}{tablename}")
    else:
        print(f"Read-based processing: Count_Tables folder for {output_type} does not exist.")

    combined_files = ",".join(table_files)

    # add column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=combined_files)
    return df


def add_read_based_krona(df: pd.DataFrame, verified_read_based: dict[str, list], 
                                   output_type: str, osdr_prefix: str = ""):
    column_name = (
        f"Parameter Value[Read-based Processing/{output_type} Outputs/Krona Reports]"
    )

    table_dir = f"{output_type}_Outputs-Krona_Reports"
    table_files = []
    if table_dir in verified_read_based:
        for tablename in verified_read_based[table_dir]:
            table_files.append(f"{osdr_prefix}{tablename}")
    else:
        print(f"Read-based processing: Krona_Reports folder for {output_type} does not exist.")

    combined_files = ",".join(table_files)

    # add column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=combined_files)
    return df


def add_read_based_humann_output(df: pd.DataFrame, verified_read_based: dict[str, list], osdr_prefix: str = ""):
    column_name = "Parameter Value[Read-based Processing/Humann Outputs/Pathway Coverage]"

    table_files = []

    # Pathway Coverage tables
    path_coverage_dir = "Humann_Outputs-Pathway_Coverages"
    if path_coverage_dir in verified_read_based:
        for file in verified_read_based[path_coverage_dir]:
            table_files.append(f"{osdr_prefix}{file}")

    else:
        print("Read-based processing: Humann_Outputs Pathway_Coverages files do not exist")

    combined_files = ",".join(table_files)
    # add column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=combined_files)

    column_name = "Parameter Value[Read-based Processing/Humann Outputs/Pathway Abundances]"
    counts_column_name = "Parameter Value[Read-based Processing/Humann Outputs/Pathway Abundances/Counts Tables]"
    plots_column_name = "Parameter Value[Read-based Processing/Humann Outputs/Pathway Abundances/Heatmaps]"
    table_files = []
    counts_files = []
    plot_files = []

    # Pathway abundances
    if "Humann_Outputs-Pathway_Abundances" in verified_read_based:
        for file in verified_read_based["Humann_Outputs-Pathway_Abundances"]:
            if file.endswith(".tsv"):
                if 'filtered' in file or 'decontam' in file:
                    counts_files.append(f"{osdr_prefix}{file}")
                else:
                    table_files.append(f"{osdr_prefix}{file}")
            elif file.endswith(".png"):
                plot_files.append(f"{osdr_prefix}{file}")
    else:
        print("Read-based processing: Humann_Outputs Pathway-abundances files do not exist")
        
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=",".join(table_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=counts_column_name, value=",".join(counts_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=plots_column_name, value=",".join(plot_files))

    # Gene Families
    column_name = "Parameter Value[Read-based Processing/Humann Outputs/Gene Families]"
    counts_column_name = "Parameter Value[Read-based Processing/Humann Outputs/Gene Families/Counts Tables]"
    plots_column_name = "Parameter Value[Read-based Processing/Humann Outputs/Gene Families/Heatmaps]"

    table_files = []
    counts_files = []
    plot_files = []

    if "Humann_Outputs-Gene_Families" in verified_read_based:
        for file in verified_read_based["Humann_Outputs-Gene_Families"]:
            if file.endswith(".tsv"):
                if 'filtered' in file or 'decontam' in file:
                    counts_files.append(f"{osdr_prefix}{file}")
                else:
                    table_files.append(f"{osdr_prefix}{file}")
            elif file.endswith(".png"):
                plot_files.append(f"{osdr_prefix}{file}")
    else:
        print("Read-based processing: Humann_Outputs Gene-Families files do not exist")
        
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=column_name, value=",".join(table_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=counts_column_name, value=",".join(counts_files))
    # add table column to dataframe with same value for all rows
    df = add_parameter_column(df=df, column_name=plots_column_name, value=",".join(plot_files))

    return df


def clean_comma_space(df):
    """
    Remove spaces after commas in all string columns of the dataframe.

    Args:
        df (pandas.DataFrame): current assay table

    Returns:
        pandas.DataFrame: updated assay table
    """
    # Loop through all columns in the dataframe
    for col in df.columns:
        # Only process string (object) columns
        if df[col].dtype == "object":
            # Replace comma-space with just comma
            df[col] = df[col].str.replace(", ", ",", regex=False)

    print("Removed spaces after commas in all string columns")
    return df


def clean_column_names(df):
    """Clean column names by removing any .# suffixes pandas adds to duplicates.

    Args:
        df: The DataFrame to clean column names

    Returns:
        The DataFrame with cleaned column names
    """
    # Create a mapping of old_name -> new_name (without .# suffix)
    name_mapping = {}
    for col in df.columns:
        # Use regex to match column names with .digits suffix
        if re.search(r"\.\d+$", col):
            # Remove the .# suffix
            base_name = re.sub(r"\.\d+$", "", col)
            name_mapping[col] = base_name

    # Rename columns using the mapping if any found
    if name_mapping:
        print(f"Cleaning {len(name_mapping)} column names by removing .# suffixes:")
        for old_name, new_name in name_mapping.items():
            print(f"  - {old_name} -> {new_name}")
        df = df.rename(columns=name_mapping)

    return df


def main():
    args = parse_args()

    # Find and parse the runsheet
    runsheet_df = load_runsheet(args.runsheet)
    validation_data = load_validation_json(args.validation_json)
    processed_outdir = Path()
    failed_assemblies = []
    if validation_data:
        processed_outdir = Path(str(validation_data.get("outdir")))
        verified_folders = get_verified_folders(validation_data, processed_outdir)
        failed_assemblies = validation_data["failed_assemblies"]
    else:
        processed_outdir = Path()
        verified_folders = defaultdict(dict)
        

    read_stats_file = None
    if args.read_stats:
        read_stats_file = Path(args.read_stats)
        if read_stats_file.exists() and read_stats_file.is_file():
            report_info(f"Using Read Stats file ({read_stats_file}) to add read counts.")
        else:
            report_warning(f"Read Stats file ({args.read_stats}) is not a valid file. No read stats will be added to the assay table.")

    suffix = args.assay_suffix
    prefix = args.output_prefix
    if prefix != "" and not prefix.endswith("_"):
        prefix = f"{prefix}_"

    assay = "Metagenomic sequencing"

    resource_category = ""
    osdr_technology_str = ""
    tech_type = "short"
    if args.sample_type == "low_biomass" and args.technology == "nanopore":
        tech_type = "long"
        if suffix != "" and suffix != "_GLlblMetag":
            report_warning(
                f"Incorrect assay suffix provided. metagenomics-lowbiomass-longread requires either an empty suffix or '_GLlblMetag', not {suffix}"
            )
        resource_category = "metagenomics-lowbiomass-longread"
        osdr_technology_str = "long read low-biomass DNA sequencing"
    elif args.sample_type == "low_biomass" and args.technology == "illumina":
        tech_type = "short"
        if suffix != "" and suffix != "_GLlbsMetag":
            report_warning(
                f"Incorrect assay suffix provided. metagenomics-lowbiomass-shortread requires either an empty suffix or '_GLlbsMetag', not {suffix}"
            )
        resource_category = "metagenomics-lowbiomass-shortread"
        osdr_technology_str = "short read low-biomass DNA sequencing"
    elif args.technology == "illumina":
        tech_type = "short"
        if suffix != "" and suffix != "_GLmetagenomics":
            report_warning(
                f"Incorrect assay suffix provided. metagenomics-lowbiomass-shortread requires either an empty suffix or '_GLmetagenomics', not {suffix}"
            )
        resource_category = "metagenomics"
        osdr_technology_str = "Whole-Genome Shotgun Sequencing"
    else:
        report_failure_and_exit(
            f"Assay type {args.sample_type} and technology {args.technology} do not have a defined assay file format."
        )

    # Create GLDS prefix for all filenames
    glds_id = args.glds_accession.upper()
    osdr_prefix = f"{glds_id}_{resource_category}_"
    glds_prefix = f"{osdr_prefix}{prefix}"

    # Find Metagenomics assay file and get its contents
    if args.isa_zip != "":
        print(f"Extracting assay table from {args.isa_zip} based on assay type: {assay} and technology type: {osdr_technology_str}")
        assay_df, assay_filename = get_assay_table_from_isa(args.isa_zip, assay, osdr_technology_str)
        print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
    else:
        print(f"Reading user-provided assay table from file: {args.assay_table}")
        assay_filename = args.assay_table
        assay_df = pd.read_csv(open(assay_filename), sep="\t")
        print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")

    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")

    # Create a mapping from assay table sample names to runsheet sample names, exit if no sample column found
    sample_col = next((col for col in assay_df.columns if "Sample Name" in col), None)
    if sample_col is None:
        report_failure_and_exit(f"Could not find 'Sample Name' column in assay table '{assay_filename}'")

    assay_sample_names = assay_df[sample_col].tolist()
    assay_df.index = assay_sample_names
    if runsheet_df is not None and "Sample Name" in runsheet_df.columns:
        sample_name_map = get_runsheet_sample_name_map(runsheet_df, assay_sample_names)
    else:
        sample_name_map = {}

    assay_sample_dropouts = []

    # Process and save assay file
    try:

        # Following the original workflow order:

        # Read counts
        if read_stats_file is not None:
            assay_df = add_read_counts_from_read_stats_file(df=assay_df, 
                                                            read_stats_file=read_stats_file)

        # Merged MultiQC Data (if it exists)
        if "Merged_Sequence_Data" in verified_folders:
            if "long" in tech_type:
                assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "raw")
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "raw", 
                                                    verified_folders["Merged_Sequence_Data"]["NanoPlots"])
            else:
                assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "HRrm")

        # Trimmed MultiQC Data (if it exists)
        if "Trimmed_Sequence_Data" in verified_folders:
            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "trimmed")
            if "long" in tech_type:
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "trimmed", 
                                                    verified_folders["Trimmed_Sequence_Data"]["NanoPlots"])

        # Filtered Sequence Data/MultiQC Reports column
        if "Filtered_Sequence_Data" in verified_folders:
            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "filtered")
            # add filtered fastq if not long-read
            if "long" not in osdr_technology_str:
                assay_df = add_fastq_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, 
                                                 verified_folders["Filtered_Sequence_Data"]["Fastq-Fastp"],
                                                 "filtered", is_paired_end)
                # gather samples with no fastq files
                assay_sample_dropouts = list(assay_df["Parameter Value[Filtered Sequence Data]"][assay_df["Parameter Value[Filtered Sequence Data]"] == ""].index)

            else:
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "filtered", 
                                                    verified_folders["Filtered_Sequence_Data"]["NanoPlots"])
        else:
            report_failure_and_exit("Filtered Sequence Data not found.")

        # HRrm if present
        if "HR-removed_Sequence_Data" in verified_folders:
            assay_df = add_fastq_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, 
                                             verified_folders["HR-removed_Sequence_Data"]["Fastq"],
                                             "HRrm", is_paired_end)
            # gather samples with no fastq files
            assay_sample_dropouts = list(assay_df["Parameter Value[HR-removed Sequence Data]"][assay_df["Parameter Value[HR-removed Sequence Data]"] == ""].index)

            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "HRrm")
            if "long" in tech_type:
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "HRrm", 
                                                    verified_folders["HR-removed_Sequence_Data"]["NanoPlots"])

        # HostRm if present
        if "Host-removed_Sequence_Data" in verified_folders:
            assay_df = add_fastq_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, 
                                             verified_folders["Host-removed_Sequence_Data"]["Fastq"],
                                             "HostRm", is_paired_end)
            # gather samples with no fastq files
            assay_sample_dropouts = list(
                assay_df["Parameter Value[Host-removed Sequence Data]"][assay_df["Parameter Value[Host-removed Sequence Data]"] == ""].index)

            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "HostRm")
            if "long" in tech_type:
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "HostRm", 
                                                    verified_folders["Host-removed_Sequence_Data"]["NanoPlots"])

        # Decontaminated if present
        if "Decontaminated_Sequence_Data" in verified_folders:
            assay_df = add_fastq_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, 
                                             verified_folders["Decontaminated_Sequence_Data"]["Fastq"],
                                             "decontam", is_paired_end)
            # gather samples with no fastq files
            assay_sample_dropouts = list(
                assay_df["Parameter Value[Decontaminated Sequence Data]"][assay_df["Parameter Value[Decontaminated Sequence Data]"] == ""].index)

            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "decontam")
            if "long" in tech_type:
                assay_df = add_nanoplot_data_column(assay_df, osdr_prefix, prefix, suffix, assay_sample_names, "decontam", 
                                                    verified_folders["Decontaminated_Sequence_Data"]["NanoPlots"])


        # Read-based Processing
        if "Read-based_Processing" in verified_folders:
            verified_read_based = verified_folders["Read-based_Processing"]
            # Kaiju
            # Counts Tables
            assay_df = add_read_based_taxonomy_counts(df=assay_df, verified_read_based=verified_read_based, 
                                                      output_type="Kaiju", osdr_prefix=osdr_prefix)
            # Barplots
            assay_df = add_read_based_taxonomy_plots(df=assay_df, verified_read_based=verified_read_based, 
                                                     output_type="Kaiju", osdr_prefix=osdr_prefix)
            # Krona Reports
            assay_df = add_read_based_krona(df=assay_df, verified_read_based=verified_read_based, 
                                            output_type="Kaiju", osdr_prefix=osdr_prefix)
            # Kraken2
            # Counts Tables
            assay_df = add_read_based_taxonomy_counts(df=assay_df, verified_read_based=verified_read_based, 
                                                      output_type="Kraken2", osdr_prefix=osdr_prefix)
            # Barplots
            assay_df = add_read_based_taxonomy_plots(df=assay_df, verified_read_based=verified_read_based, 
                                                     output_type="Kraken2", osdr_prefix=osdr_prefix)
            # Krona Reports
            assay_df = add_read_based_krona(df=assay_df, verified_read_based=verified_read_based, 
                                            output_type="Kraken2", osdr_prefix=osdr_prefix)
            # Kraken2 multiqc
            assay_df = add_multiqc_reports_column(assay_df, glds_prefix, suffix, "kraken2")

            # Metaphlan
            if not "long" in osdr_technology_str:
                # Counts Tables
                assay_df = add_read_based_taxonomy_counts(df=assay_df, verified_read_based=verified_read_based, 
                                                          output_type="Metaphlan", osdr_prefix=osdr_prefix)
                # Barplots
                assay_df = add_read_based_taxonomy_plots(df=assay_df, verified_read_based=verified_read_based, 
                                                         output_type="Metaphlan", osdr_prefix=osdr_prefix)
                # Krona Reports
                assay_df = add_read_based_krona(df=assay_df, verified_read_based=verified_read_based, 
                                                output_type="Kraken2", osdr_prefix=osdr_prefix)


            # Humann
            assay_df = add_read_based_humann_output(df=assay_df, verified_read_based=verified_read_based, 
                                                    osdr_prefix=osdr_prefix)
        else:
            report_failure_and_exit("Read-based Processing data not found.")

        if "Assembly-based_Processing" in verified_folders:
            verified_assembly_based = verified_folders["Assembly-based_Processing"]
            # add assembly based processing overview file (if it exists)
            assembly_overview_file = f"{prefix}Assembly-based-processing-overview{suffix}.tsv"

            if assembly_overview_file in verified_assembly_based[""]:
                assay_df = add_parameter_column(df=assay_df, column_name="Parameter Value[Assembly-based Processing]",
                                                value=assembly_overview_file, glds_prefix=osdr_prefix)
            else:
                report_warning("Assembly based overview not found.")

            # assemblies
            if "assemblies" in verified_assembly_based:
                # add assay table dropouts to failed assemblies list
                failed_assemblies.extend(assay_sample_dropouts)
                assay_df = add_assemblies(df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                        file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix)

            # predicted genes
            if "predicted-genes" in verified_assembly_based:
                assay_df = add_predicted_genes(df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                            file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix)

            # read-mapping
            if "read-mapping" in verified_assembly_based:
                assay_df = add_read_mapping(df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                            file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix)

            # annotations-and-taxonomy
            if "annotations-and-taxonomy" in verified_assembly_based:
                assay_df = add_annotation_and_tax( df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                                  file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix)
            # bins
            if "bins" in verified_assembly_based:
                assay_df = add_bins(df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                    bins_files=verified_assembly_based["bins"],
                                    file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix,)

            # MAGs
            if "MAGs" in verified_assembly_based:
                assay_df = add_MAGs(df=assay_df, assay_sample_names=assay_sample_names, failed_assemblies=failed_assemblies,
                                    mags_files=verified_assembly_based["MAGs"],
                                    file_prefix=prefix, osdr_prefix=osdr_prefix, assay_suffix=suffix)
                mags_ko_files = [
                    f"{osdr_prefix}{prefix}MAG-KEGG-Decoder-out{suffix}.html",
                    f"{osdr_prefix}{prefix}MAG-KEGG-Decoder-out{suffix}.tsv",
                    f"{osdr_prefix}{prefix}MAG-level-KO-annotations{suffix}.tsv",
                ]
                assay_df = add_parameter_column(df=assay_df, column_name="Parameter Value[Assembly-based Processing/MAGs/KEGG Outputs]", 
                                                value=",".join(mags_ko_files))

            # combined-outputs
            assay_df = add_combined_outputs(df=assay_df, verified_assembly_based=verified_assembly_based, 
                                            osdr_prefix=osdr_prefix)

        else:
            report_failure_and_exit("Assembly-based Processing data not found.")

        # Clean comma-space in all string columns
        assay_df = clean_comma_space(assay_df)

        # Clean column names by removing any .# suffixes pandas adds
        assay_df = clean_column_names(assay_df)

        # Use the filename we found in extract_and_find_assay
        orig_filename = assay_filename

        # Create both original and modified output files
        # Original file (preserving the original name)
        # assay_df.to_csv(orig_filename, sep='\t', index=False)
        # print(f"Original assay table saved as: {orig_filename}")

        # Modified file with GLDS prefix
        orig_filename_base = Path(orig_filename).name
        if not orig_filename_base.startswith(glds_prefix):
            mod_filename = f"{glds_prefix}{orig_filename_base}"
        else:
            mod_filename = orig_filename_base

        print(mod_filename)
        assay_df.to_csv(mod_filename, sep="\t", index=False)
        print(f"Modified assay table saved as: {mod_filename}")

    except Exception as e:
        report_failure_and_exit(str(e))


if __name__ == "__main__":
    main()
