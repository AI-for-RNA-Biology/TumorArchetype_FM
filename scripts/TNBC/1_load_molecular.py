#!/usr/bin/env python3
"""
Load molecular data from parquet files with tissue annotation filtering.

This script provides functions to efficiently load large parquet files containing
molecular data, with the ability to filter by specific tissue types using chunked processing.
"""

import pandas as pd
import pyarrow.parquet as pq
import pyarrow as pa
import argparse
import os
import sys


def load_parquet_with_spots(file_path, n_rows=None):
    """
    Load parquet file and properly set spot names as index
    
    Parameters:
    - file_path: Path to parquet file
    - n_rows: Number of rows to read (None for all rows)
    
    Returns:
    - DataFrame with spot_id as index
    """
    if n_rows:
        # Read limited rows
        parquet_file = pq.ParquetFile(file_path)
        first_batch = next(parquet_file.iter_batches(batch_size=n_rows))
        df = first_batch.to_pandas().head(n_rows)
    else:
        # Read full file (use with caution!)
        df = pd.read_parquet(file_path)
    
    # Set spot_id as index if it exists
    if 'spot_id' in df.columns:
        df = df.set_index('spot_id')
        print(f"✅ Set spot_id as index. Sample spot names: {df.index[:3].tolist()}")
    else:
        print(f"❌ No spot_id column found. Available columns: {df.columns.tolist()}")
    
    return df


def load_parquet_by_tissue_chunks(file_path, annotation, target_tissues, chunk_size=1000):
    """
    Load only spots from specific tissue types using chunked processing.
    
    Parameters:
    - file_path: Path to parquet file
    - annotation: DataFrame with spot annotations (index = spot names, 'label' column = tissue type)
    - target_tissues: List of tissue types to load (e.g., ['adipose tissue'])
    - chunk_size: Number of rows to process per chunk
    
    Returns:
    - DataFrame with only spots from target tissues
    """
    # Get target spot names
    target_spots = set(annotation[annotation['label'].isin(target_tissues)].index)
    print(f"Target tissues: {target_tissues}")
    print(f"Found {len(target_spots)} spots in target tissues")
    
    if len(target_spots) == 0:
        print("❌ No spots found for target tissues!")
        return pd.DataFrame()
    
    # Process file in chunks
    parquet_file = pq.ParquetFile(file_path)
    matching_chunks = []
    total_processed = 0
    
    print(f"Processing parquet file in chunks of {chunk_size}...")
    
    for i, batch in enumerate(parquet_file.iter_batches(batch_size=chunk_size)):
        chunk_df = batch.to_pandas()
        chunk_df['spot_id'] = chunk_df['spot_id'].apply(lambda x: x.split('parquet.')[1].replace('X', 'spot'))
        total_processed += len(chunk_df)
        
        # Check if spot_id column exists
        if 'spot_id' in chunk_df.columns:
            # Filter for target spots
            matching_spots = chunk_df[chunk_df['spot_id'].isin(target_spots)]
            if len(matching_spots) > 0:
                matching_spots = matching_spots.set_index('spot_id')
                
                matching_chunks.append(matching_spots)
                print(f"  Chunk {i+1}: Found {len(matching_spots)} matching spots")
        
        # Progress update every 10 chunks
        if (i + 1) % 10 == 0:
            print(f"  Processed {total_processed} rows...")
    
    # Combine all matching chunks
    if matching_chunks:
        result = pd.concat(matching_chunks, axis=0)
        print(f"✅ Successfully loaded {len(result)} spots from {len(matching_chunks)} chunks")
        
        # Show tissue distribution
        result_annotation = annotation.loc[result.index]
        tissue_counts = result_annotation['label'].value_counts()
        print("Tissue distribution in result:")
        for tissue, count in tissue_counts.items():
            print(f"  {tissue}: {count} spots")
        
        return result
    else:
        print("❌ No matching spots found in any chunk!")
        return pd.DataFrame()


def load_adipose_tissue_spots(file_path, annotation, chunk_size=1000):
    """
    Convenience function to load only adipose tissue spots.
    """
    return load_parquet_by_tissue_chunks(file_path, annotation, ['adipose tissue'], chunk_size)


def load_multiple_tissues(file_path, annotation, tissue_list, chunk_size=1000):
    """
    Convenience function to load spots from multiple tissue types.
    
    Parameters:
    - tissue_list: List of tissue types (e.g., ['adipose tissue', 'connective tissue'])
    """
    return load_parquet_by_tissue_chunks(file_path, annotation, tissue_list, chunk_size)


def main():
    parser = argparse.ArgumentParser(description='Load and filter molecular data by tissue type')
    parser.add_argument('input_parquet', help='Path to input parquet file')
    parser.add_argument('annotation_csv', help='Path to annotation CSV file with spot labels')
    parser.add_argument('output_parquet', help='Path to output filtered parquet file')
    parser.add_argument('--tissues', nargs='+', 
                       help='List of tissue types to include (default: all except undetermined)',
                       default=None)
    parser.add_argument('--chunk-size', type=int, default=1000,
                       help='Chunk size for processing (default: 1000)')
    parser.add_argument('--exclude-undetermined', action='store_true',
                       help='Exclude undetermined tissue spots (default: True)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.input_parquet):
        print(f"❌ Input parquet file not found: {args.input_parquet}")
        sys.exit(1)
    
    if not os.path.exists(args.annotation_csv):
        print(f"❌ Annotation CSV file not found: {args.annotation_csv}")
        sys.exit(1)
    
    # Load annotation
    print("Loading annotation data...")
    annotation = pd.read_csv(args.annotation_csv, index_col=0)
    
    # Show available tissue types
    print("\nAvailable tissue types:")
    tissue_counts = annotation['label'].value_counts()
    print(tissue_counts)
    
    # Determine tissues to load
    if args.tissues is None:
        if args.exclude_undetermined:
            tissue_list = list(annotation[annotation['label'] != 'undetermined']['label'].unique())
        else:
            tissue_list = list(annotation['label'].unique())
    else:
        tissue_list = args.tissues
    
    print(f"\nTissues to load: {tissue_list}")
    
    # Load filtered data
    print("\nLoading filtered molecular data...")
    filtered_data = load_multiple_tissues(args.input_parquet, annotation, tissue_list, args.chunk_size)
    
    if len(filtered_data) == 0:
        print("❌ No data loaded. Exiting.")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_parquet)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save filtered data
    print(f"\nSaving filtered data to: {args.output_parquet}")
    filtered_data.to_parquet(args.output_parquet)
    
    print(f"✅ Successfully saved {len(filtered_data)} spots with {len(filtered_data.columns)} features")
    print(f"Final data shape: {filtered_data.shape}")


if __name__ == "__main__":
    main()