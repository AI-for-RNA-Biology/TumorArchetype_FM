#!/usr/bin/env python3
"""
Memory-efficient processing of RDS files to create concatenated count matrix
Processes samples one by one and saves intermediate results to stay under 128GB
"""

import os
import sys
import glob
import pandas as pd
import numpy as np
import rpy2.robjects as ro
import rpy2.robjects.vectors as rvec
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import gc
import argparse
from gtfparse import read_gtf

def convert_to_df(data):
    """Convert R objects to pandas DataFrame"""
    if isinstance(data, rvec.Matrix):
        return pandas2ri.rpy2py(ro.r['as.data.frame'](data))
    elif isinstance(data, rvec.Vector):
        return pandas2ri.rpy2py(data)
    elif isinstance(data, rvec.DataFrame):
        return pandas2ri.rpy2py(data)
    else:
        return None

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Memory-efficient processing of RDS files to create concatenated count matrix'
    )
    
    parser.add_argument(
        '--rds-dir', 
        type=str,
        default="/storage/research/dbmr_luisierlab/database/ST_TNBC_v3/Robjects/countsNonCorrected",
        help='Directory containing RDS files (default: %(default)s)'
    )
    
    parser.add_argument(
        '--gtf-path',
        type=str, 
        default="/storage/research/dbmr_luisierlab/group/genomics/annotation/hg38/GENCODE/gencode.v38.annotation.gtf",
        help='Path to GTF annotation file (default: %(default)s)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default="/storage/research/dbmr_luisierlab/temp/lfournier/repositories/TumorArchetype-FM/results/molecular/TNBC",
        help='Output directory (default: %(default)s)'
    )
    
    parser.add_argument(
        '--batch-size',
        type=int,
        default=10,
        help='Number of samples to process in each batch (default: %(default)d)'
    )
    
    parser.add_argument(
        '--output-name',
        type=str,
        default="TNBC_all_samples_coding_genes",
        help='Base name for output files (default: %(default)s)'
    )
    
    parser.add_argument(
        '--save-csv',
        action='store_true',
        help='Also save output as CSV file (warning: can be very large)'
    )
    
    parser.add_argument(
        '--keep-intermediates',
        action='store_true', 
        help='Keep intermediate batch files after processing'
    )
    
    return parser.parse_args()

def load_annotation(gtf_path):
    """Load and filter annotation for protein-coding genes"""
    print("Loading GTF annotation...")
    annotation = read_gtf(gtf_path)
    
    # Convert to pandas if it's Polars
    if hasattr(annotation, 'to_pandas'):
        annotation = annotation.to_pandas()
    
    # Get unique protein-coding gene names
    protein_coding = annotation[annotation['gene_type'] == 'protein_coding']
    coding_genes = set(protein_coding['gene_name'].unique())
    
    print(f"Found {len(coding_genes)} unique protein-coding genes")
    return coding_genes

def process_single_rds(rds_path, coding_genes):
    """Process a single RDS file and return filtered count matrix"""
    sample_name = os.path.basename(rds_path).replace('.RDS', '')
    print(f"Processing {sample_name}...")
    
    try:
        # Load RDS file
        df = ro.r['readRDS'](rds_path)
        
        # Extract count matrix
        counts_matrix = df.rx2('cnts')
        counts_df = convert_to_df(counts_matrix)
        
        # Clear R objects from memory
        del df, counts_matrix
        ro.r('gc()')  # R garbage collection
        
        if counts_df is None:
            print(f"Failed to convert {sample_name}")
            return None
            
        print(f"  Original shape: {counts_df.shape}")
        
        # Filter for protein-coding genes only
        available_coding_genes = [gene for gene in counts_df.columns if gene in coding_genes]
        counts_df_coding = counts_df[available_coding_genes].copy()
        
        # Clear original dataframe
        del counts_df
        gc.collect()
        
        print(f"  Filtered shape: {counts_df_coding.shape}")
        
        # Rename index to include sample name
        counts_df_coding.index = [f"{sample_name}_{idx}" for idx in counts_df_coding.index]
        
        return counts_df_coding, set(available_coding_genes)
        
    except Exception as e:
        print(f"Error processing {rds_path}: {e}")
        return None

def process_batch(rds_files, coding_genes, batch_num, intermediate_dir):
    """Process a batch of RDS files and save intermediate result"""
    print(f"\n=== Processing Batch {batch_num} ({len(rds_files)} files) ===")
    
    batch_dfs = []
    all_genes_in_batch = set()
    
    for rds_path in rds_files:
        result = process_single_rds(rds_path, coding_genes)
        if result is not None:
            counts_df, genes_in_sample = result
            batch_dfs.append(counts_df)
            all_genes_in_batch.update(genes_in_sample)
            
            # Force garbage collection
            gc.collect()
    
    if not batch_dfs:
        print(f"No valid samples in batch {batch_num}")
        return None
    
    print(f"Concatenating {len(batch_dfs)} samples in batch {batch_num}...")
    
    # Concatenate samples in this batch
    batch_combined = pd.concat(batch_dfs, axis=0, sort=False)
    
    # Clear individual dataframes
    del batch_dfs
    gc.collect()
    
    print(f"Batch {batch_num} combined shape: {batch_combined.shape}")
    
    # Save intermediate result
    intermediate_path = os.path.join(intermediate_dir, f"batch_{batch_num:02d}.parquet")
    batch_combined.to_parquet(intermediate_path)
    print(f"Saved batch {batch_num} to {intermediate_path}")
    
    return batch_combined, all_genes_in_batch

def combine_batches(intermediate_dir, output_dir, output_name, save_csv=True):
    """Combine all batch files into final matrix"""
    print("\n=== Combining all batches ===")
    
    # Find all batch files
    batch_files = sorted(glob.glob(os.path.join(intermediate_dir, "batch_*.parquet")))
    print(f"Found {len(batch_files)} batch files")
    
    if not batch_files:
        print("No batch files found!")
        return
    
    # Load first batch to get all gene names
    print("Loading batches to determine all genes...")
    all_genes = set()
    batch_shapes = []
    
    for batch_file in batch_files:
        print(f"  Checking {os.path.basename(batch_file)}...")
        batch_df = pd.read_parquet(batch_file)
        all_genes.update(batch_df.columns)
        batch_shapes.append(batch_df.shape)
        del batch_df
        gc.collect()
    
    all_genes = sorted(list(all_genes))
    print(f"Total unique genes across all batches: {len(all_genes)}")
    print(f"Batch shapes: {batch_shapes}")
    
    # Combine batches one by one, ensuring all have the same columns
    combined_dfs = []
    
    for i, batch_file in enumerate(batch_files):
        print(f"Loading batch {i+1}/{len(batch_files)}: {os.path.basename(batch_file)}")
        batch_df = pd.read_parquet(batch_file)
        
        # Reindex to ensure all batches have the same columns
        batch_df = batch_df.reindex(columns=all_genes, fill_value=0.0)
        combined_dfs.append(batch_df)
        
        # Monitor memory usage
        if i % 5 == 0:
            gc.collect()
    
    # Final concatenation
    print("Final concatenation...")
    final_df = pd.concat(combined_dfs, axis=0, sort=False)
    
    # Clear intermediate data
    del combined_dfs
    gc.collect()
    
    print(f"Final matrix shape: {final_df.shape}")
    
    # Save final result
    output_path = os.path.join(output_dir, f"{output_name}.parquet")
    final_df.to_parquet(output_path)
    print(f"Final matrix saved to: {output_path}")
    
    # Also save as CSV if requested
    if save_csv:
        csv_path = os.path.join(output_dir, f"{output_name}.csv")
        print(f"Saving CSV version to: {csv_path}")
        final_df.to_csv(csv_path)
    
    return final_df

def main():
    print("=== Memory-Efficient RDS Processing ===")
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Create output directories
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(intermediate_dir, exist_ok=True)
    
    print(f"RDS directory: {args.rds_dir}")
    print(f"GTF file: {args.gtf_path}")
    print(f"Output directory: {args.output_dir}")
    print(f"Batch size: {args.batch_size}")
    print(f"Output name: {args.output_name}")
    print(f"Save CSV: {args.save_csv}")
    print(f"Keep intermediates: {args.keep_intermediates}")
    
    # Load annotation
    coding_genes = load_annotation(args.gtf_path)
    
    # Find all RDS files
    rds_pattern = os.path.join(args.rds_dir, "*.RDS")
    rds_files = sorted(glob.glob(rds_pattern))
    print(f"Found {len(rds_files)} RDS files")
    
    if not rds_files:
        print(f"No RDS files found in {args.rds_dir}")
        return
    
    # Process in batches
    all_genes_seen = set()
    
    for i in range(0, len(rds_files), args.batch_size):
        batch_files = rds_files[i:i+args.batch_size]
        batch_num = (i // args.batch_size) + 1
        
        result = process_batch(batch_files, coding_genes, batch_num, intermediate_dir)
        if result is not None:
            batch_df, genes_in_batch = result
            all_genes_seen.update(genes_in_batch)
            # Clear batch data from memory
            del batch_df
            gc.collect()
    
    print(f"\nTotal unique coding genes seen: {len(all_genes_seen)}")
    
    # Combine all batches
    final_df = combine_batches(intermediate_dir, args.output_dir, args.output_name, args.save_csv)
    
    # Clean up intermediate files unless requested to keep them
    if not args.keep_intermediates:
        print("\nCleaning up intermediate files...")
        for batch_file in glob.glob(os.path.join(intermediate_dir, "batch_*.parquet")):
            os.remove(batch_file)
            print(f"Removed {os.path.basename(batch_file)}")
        
        os.rmdir(intermediate_dir)
        print("Cleanup complete!")
    else:
        print(f"\nIntermediate files kept in: {intermediate_dir}")
    
    print(f"\n=== Processing Complete ===")
    print(f"Final matrix shape: {final_df.shape}")
    print(f"Output saved to: {args.output_dir}")

if __name__ == "__main__":
    main()