#!/usr/bin/env python3

# Package imports
import argparse
import os
import re
import math

# Aliased imports
import pandas as pd
import numpy as np

# Plotting imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib_inline.backend_inline
from matplotlib.colors import TABLEAU_COLORS

# Plotting settings
matplotlib_inline.backend_inline.set_matplotlib_formats('svg')
default_dpi = mpl.rcParams['figure.dpi']
mpl.rcParams['figure.dpi'] = 200

# Output root
TLD_PATH = 'evidence_qc'

# PED file validation
ID_TYPE_SAMPLE = "sample"
ID_TYPE_FAMILY = "family"
ID_TYPE_PARENT = "parent"
FIELD_NUMBER_ID = 1
FIELD_NUMBER_SEX = 4
ILLEGAL_ID_SUBSTRINGS = ["chr", "name", "DEL", "DUP", "CPX", "CHROM"]

# Helper functions
def generate_file_path(tld_path, file_type, file_name):
    """
    Calculate the output file path in the local filesystem.

    Args:
        tld_path (str): Top-level directory path.
        file_type (str): Enables generation of the specific sub-directory which a file should live in.
        file_name (str): File name to chain at the end of the path.

    Returns:
        str: Path to file as it should be saved, per file system outline.
    """
    return tld_path + '/' + file_type + '/' + file_name

def save_df(file_path, df):
    """
    Save a dataframe to the specified file path, creating directories if they don't exist.
    
    Args:
        file_path (str): The path where the figure should be saved.
        df (pandas.DataFrame): The dataframe to save.
    
    Returns:
        None.
    """
    dir_path = os.path.dirname(file_path)
    
    if dir_path:
        os.makedirs(dir_path, exist_ok=True)
    
    df.to_csv(file_path, sep='\t', index=False)

    print(f"File saved to: {file_path}")

def save_figure(file_path, fig=None):
    """
    Save a figure to the specified file path, creating directories if they don't exist.
    
    Args:
        file_path (str): The path where the figure should be saved.
        fig (matplotlib.figure.Figure, optional): The figure to save. If None, uses the current figure.
    
    Returns:
        None.
    """
    dir_path = os.path.dirname(file_path)
    
    if dir_path:
        os.makedirs(dir_path, exist_ok=True)
    
    if fig is None:
        plt.savefig(file_path)
    else:
        fig.savefig(file_path)

def write_batch_assignments(batches):
    """
    Creates a two-column DataFrame with batch assignments and writes it to a file.

    Args:
        batches (dict): Dictionary containing batch assignment information.

    Returns:
        None.
    """
    # Create batches dataframe
    batch_assignments = []
    for batch_name in sorted(batches.keys()):
        sample_list = batches[batch_name]
        batch_df = pd.DataFrame({
            'batch': [batch_name] * len(sample_list),
            'sample': sample_list
        })
        batch_assignments.append(batch_df)
    
    # Create and write the dataframe
    batch_assignment_df = pd.concat(batch_assignments, ignore_index=True)
    file_path = generate_file_path(TLD_PATH, "batching", "batch_assignments.tsv")
    save_df(file_path, batch_assignment_df)

# Validation
def validate_numeric_inputs(input_vals, log=True):
    """
    Validates user input to check whether it is all numeric or not. 
    
    Args:
        input_vals (list): List of user input to validate.
    
    Returns:
        None.
    """
    for i in input_vals:
        if not isinstance(i, (int, float)) or not i:
            raise Exception('Value input must be numeric.')
    
    if log:
        print("Inputs are valid - please proceed to the next cell.")

def validate_string_inputs(input_vals, log=True):
    """
    Validates user input to check whether it is all strings or not. 
    
    Args:
        input_vals (list): List of user input to validate.
    
    Returns:
        None.
    """
    for i in input_vals:
        if not isinstance(i, str):
            raise Exception('Value input must be a string.')
    
    if log:
        print("Inputs are valid - please proceed to the next cell.")

def validate_include_metrics(include_metrics, include_cols):
    """
    Validates user input to check whether they are valid metrics to include in batching. 
    
    Args:
        include_metrics (list): List of metrics to validate.
        include_cols (list): List of available metrics.
    
    Returns:
        None.
    """
    validate_string_inputs(include_metrics, log=False)
    
    if (not len(include_metrics) > 0):
        raise Exception(f"Length of INCLUDE_METRICS must be at least 1.")
    
    for i in include_metrics:
        if not i in include_cols:
            raise Exception(f"Metric '{i}' not available - please only use metrics from above.")
    
    print("Inputs are valid - please proceed to the next cell.")

def validate_include_bins(include_metrics, include_bins):
    """
    Validates user input to check whether they are valid bins to use when batching. 
    
    Args:
        include_metrics (list): List of metrics to validate.
        include_bins (list): List of metric bins to validate.
    
    Returns:
        None.
    """
    validate_numeric_inputs(include_bins, log=False)
    
    if (len(include_metrics) != len(include_bins) + 1):
        raise Exception(f"Length of INCLUDE_BINS must be 1 less than the length of INCLUDE_METRICS.")
    
    print("Inputs are valid - please proceed to the next cell.")

def validate_batch_sizes(df, target_batch_size, min_batch_size, max_batch_size):
    """
    Validates user input to check whether they are valid bins to use when batching. 
    
    Args:
        df (pd.DataFrame): Dataframe with sample data.
        target_batch_size (int): Target size of each batch.
        min_batch_size (int): Minimum size of each batch.
        max_batch_size (int): Maximum size of each batch.
    
    Returns:
        None.
    """
    validate_numeric_inputs([target_batch_size, min_batch_size, max_batch_size], log=False)
    
    if (target_batch_size < min_batch_size):
        raise Exception("TARGET_BATCH_SIZE must exceed MIN_BATCH_SIZE.")
        
    if (max_batch_size < target_batch_size):
        raise Exception("MAX_BATCH_SIZE must exceed TARGET_BATCH_SIZE.")
    
    if (len(df) < min_batch_size):
        raise Exception("MIN_BATCH_SIZE must exceed the number of samples.")
    
    print("Inputs are valid - please proceed to the next cell.")

def validate_id(identifier, id_type, source_file):
    """
    Validates sample IDs provided based on a source file of samples.
    
    Args:
        identifier (str): ID for a given sample.
        id_type (str): Type of ID provided.
        source_file (str): File that contains all samples.
    
    Returns:
        None.
    """
    # Check for empty IDs
    if identifier is None or identifier == "":
        raise ValueError(f"Empty {id_type} ID in {source_file}.")

    # Check all characters are alphanumeric or underscore
    if not re.match(r'^\w+$', identifier):
        raise ValueError(f"Invalid {id_type} ID in {source_file}: '{identifier}'." + 
                         "IDs should only contain alphanumeric and underscore characters.")

    # Check for all-numeric IDs, besides maternal & paternal ID (can be 0) and all-numeric family IDs
    if id_type != ID_TYPE_FAMILY and not (id_type == ID_TYPE_PARENT and identifier == "0") and identifier.isdigit():
        raise ValueError(f"Invalid {id_type} ID in {source_file}: {identifier}. " +
                         "IDs should not contain only numeric characters.")

    # Check for illegal substrings
    for sub in ILLEGAL_ID_SUBSTRINGS:
        if sub in identifier:
            raise ValueError(f"Invalid {id_type} ID in {source_file}: {identifier}. " +
                             f"IDs cannot contain the following substrings: {', '.join(ILLEGAL_ID_SUBSTRINGS)}.")

def validate_ped(ped_file, samples):
    """
    Validates structure and data within PED file based on series of samples provided.
    
    Args:
        ped_file (str): Path to local PED file.
        samples (set): Set of sample IDs to validate against.
    
    Returns:
        None
    """
    seen_sex_1 = False
    seen_sex_2 = False
    samples_found = set()

    # Read PED file
    try:
        df = pd.read_table(ped_file, dtype=str, header=None, comment='#', names=[
            'Family_ID', 'Sample_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype'
        ])
    except Exception as e:
        raise ValueError(f"Error reading PED file: {str(e)}")
    
    # Ensure column count
    if len(df.columns) != 6:
        raise ValueError("PED file must have 6 columns: Family_ID, Sample_ID, " +
                         "Paternal_ID, Maternal_ID, Sex, Phenotype.")

    # Iteratively validate each row
    for _, row in df.iterrows():
        # Validate ID
        for identifier, id_type in zip(row[:FIELD_NUMBER_SEX],
                                       [ID_TYPE_FAMILY, ID_TYPE_SAMPLE, ID_TYPE_PARENT, ID_TYPE_PARENT]):
            validate_id(identifier, id_type, "PED file")

        # Assign main information to variables
        sample_id = row['Sample_ID']
        sex = int(row['Sex'])

        # Check for appearance of each sex
        if sex == 1:
            seen_sex_1 = True
        elif sex == 2:
            seen_sex_2 = True
        elif sex != 0:
            raise ValueError(f"Sample {sample_id} has an invalid value for sex: {sex}. " +
                             "PED file must use the following values for sex: " + 
                             "1 for Male, 2 for Female, 0 for Unknown/Other.")

        # Verify no duplications
        if sample_id in samples_found:
            raise ValueError(f"Duplicate entries for sample {sample_id}.")
        elif sample_id in samples:
            samples_found.add(sample_id)

    # Check if all samples in the sample list are present in PED file
    if len(samples_found) < len(samples):
        missing = samples - samples_found
        raise ValueError(f"PED file is missing sample(s): {','.join(missing)}.")

    # Raise error if at least one of either sex is not found
    if not (seen_sex_2 and seen_sex_1):
        raise ValueError("Did not find existence of multiple sexes in file. "  +
                         "PED file must use the following values for sex: " + 
                         "1 for Male, 2 for Female, 0 for Unknown/Other.")
    
    print("PED file is valid - please proceed to the next cell.")

# Batching
def plot_batch_metrics(df, batches, include_metrics, display=False, prefix=''):
    """
    Plots graphs per batch regarding the distribution of included metrics.
    
    Args:
        df (pandas.DataFrame): Dataframe containing sample data.
        batches (list): Contains names of batches generated.
        include_metrics (list): Contains metrics included in batching process.
        prefix (str): Prefix to use when naming output files.
    
    Returns:
        None
    """
    fig, ax = plt.subplots(nrows=len(batches), ncols=len(include_metrics), 
                           figsize=(5*len(include_metrics), 5*len(batches)), sharex='col', sharey='col')
    titlesize = 18
    labelsize = 16
    
    # Determine global min and max for each metric using df
    global_min = df[include_metrics].min()
    global_max = df[include_metrics].max()
    
    for i, (batch_name, batch_samples) in enumerate(batches.items()):
        batch_samples = set(batch_samples)
        batch_meta = df[df.sample_id.isin(batch_samples)]
        
        for j, metric in enumerate(include_metrics):
            data = batch_meta[metric].values
            bins = np.linspace(global_min[metric], global_max[metric], 30)
            
            if len(batches) == 1:
                current_ax = ax[j]
            else:
                current_ax = ax[i][j]
            
            current_ax.hist(data, bins=bins, color='gray', edgecolor='black')
            current_ax.set_xlabel(metric, fontsize=labelsize)
            current_ax.set_ylabel("Sample count", fontsize=labelsize)
            current_ax.set_xlim(global_min[metric], global_max[metric])
            current_ax.xaxis.set_tick_params(labelbottom=True, size=labelsize)
        
        if len(batches) == 1:
            ax[0].set_title(batch_name, fontsize=titlesize)
        else:
            ax[i][0].set_title(batch_name, fontsize=titlesize)
    
    plt.tight_layout()
    
    # Save the plot as an image
    file_name = f"{prefix}_distribution_{len(batches)}_batches.png"
    file_path = generate_file_path(TLD_PATH, "batching/metric_plots", file_name)
    save_figure(file_path)
    plt.close()
    
    # Supplementary Plot
    plot_cluster_distribution(df, include_metrics, batches, prefix)
    
    # Supplementary Plot
    plot_panelled_cluster_distribution(df, include_metrics, batches, prefix)
        
    # Show figure if display = True
    if display:
        img = mpimg.imread(file_path)
        plt.figure(figsize=(10, 20))
        plt.imshow(img)
        plt.axis('off')
        plt.show()

def plot_cluster_distribution(df, include_metrics, batches, prefix=''):
    """
    Create a plot showing the batch distribution based on the first two metrics in include_metrics.
    
    Args:
        df (pd.DataFrame): DataFrame containing sample data.
        include_metrics (list): List of metric names considered for batching.
        batches (dict): Dictionary of batch names and lists of sample IDs.
        prefix (str): Prefix to use when naming output files.
    """
    x_metric, y_metric = include_metrics[:2]
    
    # Create a color map for batches
    color_list = list(TABLEAU_COLORS.values())
    batch_colors = {batch: color_list[i % len(color_list)] for i, batch in enumerate(batches.keys())}
    
    # Original combined plot
    plt.figure(figsize=(10, 8))
    
    for batch_name, sample_ids in batches.items():
        batch_data = df[df['sample_id'].isin(sample_ids)]
        plt.scatter(batch_data[x_metric], batch_data[y_metric], label=batch_name, 
                    color=batch_colors[batch_name], alpha=0.7)
    
    plt.xlabel(x_metric)
    plt.ylabel(y_metric)
    plt.title(f"Batch Distribution: {x_metric} vs {y_metric}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    
    # Save the combined plot as an image
    file_name = f"{prefix}_distribution_{len(batches)}_batches.png"
    file_path = generate_file_path(TLD_PATH, "batching/cluster_plots", file_name)
    save_figure(file_path)
    plt.close()

def plot_panelled_cluster_distribution(df, include_metrics, batches, prefix=''):
    """
    Create a pannelled plot showing the batch distribution based on the first two metrics in include_metrics.
    
    Args:
        df (pd.DataFrame): DataFrame containing sample data.
        include_metrics (list): List of metric names considered for batching.
        batches (dict): Dictionary of batch names and lists of sample IDs.
        prefix (str): Prefix to use when naming output files.
    """
    x_metric, y_metric = include_metrics[:2]
    
    # Create a color map for batches
    color_list = list(TABLEAU_COLORS.values())
    batch_colors = {batch: color_list[i % len(color_list)] for i, batch in enumerate(batches.keys())}
    
    # New panel plot
    plt.figure(figsize=(10, 8))
    n_batches = len(batches)
    n_cols = min(3, n_batches)  # Maximum 3 columns
    n_rows = math.ceil(n_batches / n_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
    fig.suptitle(f"Individual Batch Distributions: {x_metric} vs {y_metric}", fontsize=16)
    
    axes = axes.flatten() if n_batches > 1 else [axes]
    
    # Find global min and max for consistent axis limits
    x_min, x_max = df[x_metric].min(), df[x_metric].max()
    y_min, y_max = df[y_metric].min(), df[y_metric].max()
    
    for i, (batch_name, sample_ids) in enumerate(batches.items()):
        batch_data = df[df['sample_id'].isin(sample_ids)]
        axes[i].scatter(batch_data[x_metric], batch_data[y_metric], 
                        color=batch_colors[batch_name], alpha=0.7)
        axes[i].set_title(batch_name)
        axes[i].set_xlabel(x_metric)
        axes[i].set_ylabel(y_metric)
        axes[i].set_xlim(x_min, x_max)
        axes[i].set_ylim(y_min, y_max)
    
    # Remove any unused subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    
    # Save the panel plot as an image
    file_name = f"{prefix}_distribution_{len(batches)}_batches.png"
    file_path = generate_file_path(TLD_PATH, "batching/pannelled_cluster_plots", file_name)
    save_figure(file_path)
    plt.close()

def group_related_samples(df, ped_df):
    """
    Groups samples that are related to other samples within a cohort.
    
    Args:
        df (pandas.DataFrame): Dataframe that contains sample data.
        ped_df (pandas.DataFrame): Dataframe that contains family structure data derived from the PED file.
        
    Returns:
        Dictionary that maps samples to family members that are related to them.
    """
    if ped_df is None:
        return {sample_id: {sample_id} for sample_id in df.sample_id}
    
    related_groups = {}
    for _, row in ped_df.iterrows():
        sample_id = row['Sample_ID']
        if sample_id not in df.sample_id.values:
            continue
        related_ids = {sample_id, row['Paternal_ID'], row['Maternal_ID']}
        related_ids = {id for id in related_ids if id in df.sample_id.values}
        
        for id in related_ids:
            if id in related_groups:
                related_groups[id].update(related_ids)
            else:
                related_groups[id] = related_ids
    
    merged_groups = []
    for group in related_groups.values():
        for merged in merged_groups:
            if group & merged:
                merged.update(group)
                break
        else:
            merged_groups.append(group)
    
    return {sample_id: next(group for group in merged_groups if sample_id in group) for sample_id in df.sample_id}

def select_family_representatives(df, ped_df, related_groups):
    """
    Selects a set of representatives for each family group part of a cohort.
    
    Args:
        ped_df (pandas.DataFrame): Dataframe that contains family structure data.
        related_groups (dict): Dictionary mapping samples to their families.
        
    Returns:
        Set of representatives to use as a proxy for unique families when batching.
    """
    representatives = set()
    parents = set(ped_df['Paternal_ID']).union(set(ped_df['Maternal_ID']))
    parents.discard('0')
    
    for group in related_groups.values():
        lowest_descendants = [member for member in group if member not in parents]
        if lowest_descendants:
            group_df = ped_df[ped_df['Sample_ID'].isin(lowest_descendants)]
            phenotype_2 = group_df[group_df['Phenotype'] == 2]
            rep = group_df['Sample_ID'].iloc[0] if phenotype_2.empty else phenotype_2['Sample_ID'].iloc[0]
        else:
            rep = min(group, key=lambda x: df[df.sample_id == x].index[0])
        representatives.add(rep)
    
    return representatives

def split_dataframe(df, metrics, bins, batches_per_level):
    """
    Splits samples hierarchically.
    
    Args:
        df (pandas.DataFrame): Dataframe that contains sample data.
        metrics (list): Metrics to split on.
        bins (list): Number of bins to split into for each level except the last.
        batches_per_level (list): Number of batches for each level.
        
    Returns:
        List of split dataframes.
    """
    if not metrics:
        return [df]

    metric = metrics[0]
    if len(metrics) == 1:
        return np.array_split(df.sort_values(by=metric), batches_per_level[-1])

    n_bins = bins[0]
    splits = np.array_split(df.sort_values(by=metric), n_bins)
    return [subsplit for split in splits for subsplit in split_dataframe(split, metrics[1:], bins[1:], batches_per_level[1:])]

def generate_hierarchical_batches(df, include_metrics, include_bins, target_batch_size, min_batch_size, max_batch_size, batch_prefix, batch_suffix, reference_ped=None):
    """
    Generate hierarchical batches based on specified metrics and parameters, with family-based batching.
    
    Args:
        df (pandas.DataFrame): Dataframe that contains sample data.
        include_metrics (list): Metrics to consider for batching.
        include_bins (list): List of number of bins for each metric, except the last.
        target_batch_size (int): Target size for each batch.
        min_batch_size (int): Minimum allowable batch size.
        max_batch_size (int): Maximum allowable batch size.
        batch_prefix (str): Prefix for batch names.
        batch_suffix (str): Suffix for batch names.
        reference_ped (pandas.DataFrame): Dataframe that contains family structure data.
        
    Returns:
        Tuple of batches dictionary and batches metadata dictionary.
    """
    # Validation
    estimated_samples_per_split = len(df)
    for bins in include_bins:
        estimated_samples_per_split /= bins
        if estimated_samples_per_split < min_batch_size:
            raise Exception(f"Based on INCLUDE_BINS, expected samples per split ({estimated_samples_per_split}) " +
                            "is less than MIN_BATCH_SIZE ({min_batch_size}). Please adjust parameters accordingly.")
            
    # Group related samples and select family representatives
    related_groups = group_related_samples(df, reference_ped)
    family_representatives = select_family_representatives(df, reference_ped, related_groups) if reference_ped is not None else set(df['sample_id'])
    df_unrelated = df[df['sample_id'].isin(family_representatives)]
    
    # Split by gender
    isfemale = (df_unrelated.chrX_CopyNumber_rounded >= 2)
    male_df = df_unrelated[~isfemale]
    female_df = df_unrelated[isfemale]
    
    # Calculate batches_per_level
    n_samples = len(df_unrelated)
    batches_per_level = [max(1, int(np.round(n_samples / target_batch_size)))]
    for bins in include_bins:
        n_samples = int(np.round(n_samples / bins))
        batches_per_level.append(max(1, int(np.round(n_samples / target_batch_size))))
    
    # Create hierarchical splits
    male_splits = split_dataframe(male_df, include_metrics, include_bins, batches_per_level)
    female_splits = split_dataframe(female_df, include_metrics, include_bins, batches_per_level)
    combined_splits = [pd.concat([m, f]) for m, f in zip(male_splits, female_splits)]
    
    # Create batches
    batches = {}
    batches_meta = {}
    batch_num = 0
    for split_index, split in enumerate(combined_splits):
        batch_samples = split['sample_id'].tolist()
        while batch_samples:
            batch_num += 1
            batch_name = f"{batch_prefix}{batch_num}{batch_suffix}"
            batch_size = min(max(min_batch_size, len(batch_samples)), max_batch_size)
            batches[batch_name] = batch_samples[:batch_size]
            batch_samples = batch_samples[batch_size:]
            batches_meta[batch_name] = tuple(
                [split_index // (len(combined_splits) // bins) + 1 for bins in include_bins] 
                + [split_index % batches_per_level[-1] + 1]
            )
    
    # Add related individuals to proband's batch
    unbatched_samples = set()
    for batch_name, sample_ids in batches.items():
        for sample_id in sample_ids.copy():
            if sample_id in related_groups:
                related_samples = related_groups[sample_id] - set(sample_ids)
                space_left = max_batch_size - len(batches[batch_name])
                batches[batch_name].extend(list(related_samples)[:space_left])
                unbatched_samples.update(list(related_samples)[space_left:])
    
    # Rebalance batches if they exceed max_batch_size
    while unbatched_samples:
        eligible_batches = [b for b in batches if len(batches[b]) < max_batch_size]
        if not eligible_batches:
            batch_num += 1
            new_batch_name = f"{batch_prefix}{batch_num}{batch_suffix}"
            batches[new_batch_name] = []
            eligible_batches = [new_batch_name]
        smallest_batch = min(eligible_batches, key=lambda x: len(batches[x]))
        sample_to_add = unbatched_samples.pop()
        batches[smallest_batch].append(sample_to_add)
    
    # Final adjustment to ensure within batch size bounds
    while True:
        batches_to_adjust = [b for b in batches if len(batches[b]) < min_batch_size or len(batches[b]) > max_batch_size]
        if not batches_to_adjust:
            break
        
        sorted_batches = sorted(batches.items(), key=lambda x: len(x[1]))
        for batch_name in batches_to_adjust:
            if len(batches[batch_name]) < min_batch_size:
                samples_needed = min_batch_size - len(batches[batch_name])
                for large_batch_name, large_batch_samples in reversed(sorted_batches):
                    if len(large_batch_samples) > min_batch_size:
                        samples_to_move = min(samples_needed, len(large_batch_samples) - min_batch_size)
                        batches[batch_name].extend(large_batch_samples[-samples_to_move:])
                        batches[large_batch_name] = large_batch_samples[:-samples_to_move]
                        samples_needed -= samples_to_move
                        if samples_needed == 0:
                            break
            
            elif len(batches[batch_name]) > max_batch_size:
                overflow = batches[batch_name][max_batch_size:]
                batches[batch_name] = batches[batch_name][:max_batch_size]
                for small_batch_name, small_batch_samples in sorted_batches:
                    if len(small_batch_samples) < max_batch_size:
                        space_available = max_batch_size - len(small_batch_samples)
                        samples_to_move = min(space_available, len(overflow))
                        batches[small_batch_name].extend(overflow[:samples_to_move])
                        overflow = overflow[samples_to_move:]
                        if not overflow:
                            break
                if overflow:
                    batch_num += 1
                    new_batch_name = f"{batch_prefix}{batch_num}{batch_suffix}"
                    batches[new_batch_name] = overflow
        sorted_batches = sorted(batches.items(), key=lambda x: len(x[1]))
    
    return batches, batches_meta



def _parse_csv_list(value):
    return [x for x in (value or "").split(",") if x]


def main():
    parser = argparse.ArgumentParser(
        description="Batch generation from EvidenceQC metadata")
    parser.add_argument("--pass-metadata", required=True, help="Path to passing_samples_metadata.tsv from EvidenceQC")
    parser.add_argument("--include-metrics", required=True, help="Comma-separated metric columns")
    parser.add_argument("--include-bins", default="", help="Comma-separated bins (length = metrics-1)")
    parser.add_argument("--target-batch-size", required=True, type=int, help="Target batch size")
    parser.add_argument("--min-batch-size", required=True, type=int, help="Minimum batch size")
    parser.add_argument("--max-batch-size", required=True, type=int, help="Maximum batch size")
    parser.add_argument("--batch-prefix", default="batch", help="Batch prefix")
    parser.add_argument("--batch-suffix", default="", help="Batch suffix")
    parser.add_argument("--ped-file", default=None, help="Path to reference PED file (only to batch related samples)")
    parser.add_argument("--outdir", default="evidence_qc", help="Output top-level directory")
    parser.add_argument("--plot-prefix", default="fig_hierarchical")
    parser.add_argument("--display-plots", action="store_true")
    args = parser.parse_args()

    global TLD_PATH
    TLD_PATH = args.outdir

    include_metrics = _parse_csv_list(args.include_metrics)
    include_bins = [int(x) for x in _parse_csv_list(args.include_bins)]

    pass_df = pd.read_table(args.pass_metadata)

    # Duplicate sample_id guard
    id_counts = pass_df['sample_id'].value_counts()
    duplicates_dict = id_counts[id_counts > 1].to_dict()
    if duplicates_dict:
        for sample_id, count in duplicates_dict.items():
            print(f"Sample ID: {sample_id}, Count: {count}")
        raise Exception("Batching requires unique 'sample_id' - please resolve duplicates first.")

    include_cols = [col for col in pass_df.columns if col not in ('sample_id', 'chrX_CopyNumber_rounded')]
    validate_include_metrics(include_metrics, include_cols)
    validate_include_bins(include_metrics, include_bins)
    validate_batch_sizes(pass_df, args.target_batch_size, args.min_batch_size, args.max_batch_size)
    validate_string_inputs([args.batch_prefix, args.batch_suffix])

    reference_ped = None
    if args.ped_file:
        validate_ped(args.ped_file, set(pass_df['sample_id']))
        reference_ped = pd.read_table(args.ped_file, dtype=str, header=None, names=[
            "Family_ID", "Sample_ID", "Paternal_ID", "Maternal_ID", "Sex", "Phenotype"
        ])
        reference_ped = reference_ped[reference_ped['Sample_ID'].isin(pass_df['sample_id'])]

    batches, batches_meta = generate_hierarchical_batches(
        pass_df,
        include_metrics,
        include_bins,
        args.target_batch_size,
        args.min_batch_size,
        args.max_batch_size,
        args.batch_prefix,
        args.batch_suffix,
        reference_ped
    )

    for batch_name, sample_ids in batches.items():
        batch_df = pass_df[pass_df['sample_id'].isin(sample_ids)]
        female_count = (batch_df.chrX_CopyNumber_rounded >= 2).sum()
        male_count = (batch_df.chrX_CopyNumber_rounded < 2).sum()
        print(f"{batch_name}: {len(sample_ids)} samples ({male_count} male, {female_count} female)")

    plot_batch_metrics(pass_df, batches, include_metrics, display=args.display_plots, prefix=args.plot_prefix)

    tsv_meta = pd.DataFrame.from_dict(batches_meta, orient='index', columns=include_metrics)
    tsv_meta['n_samples'] = tsv_meta.index.map(lambda x: len(batches[x]))
    tsv_meta = tsv_meta.reset_index().rename(columns={'index': 'batch'})
    file_path = generate_file_path(TLD_PATH, 'batching', "batching_metadata.tsv")
    save_df(file_path, tsv_meta)
    write_batch_assignments(batches)


if __name__ == "__main__":
    main()