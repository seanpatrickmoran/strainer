#!/usr/bin/env python3
"""
data2json.py - Enhanced database to JSON converter with feature source tracking
Prepares Hi-C image data for ML pipeline with complete metadata preservation
"""

import sqlite3
from pathlib import Path
import sys
import json
import numpy as np
from PIL import Image
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import datetime
import argparse
from typing import List, Dict, Optional, Tuple

def colormap():
    """Default colormap for Hi-C visualization"""
    colors = [(1, 1, 1), (1, 0, 0)]  # White to red
    return colors

def call(PATH: str, TIMEOUT: int = 10) -> Tuple[sqlite3.Connection, sqlite3.Cursor]:
    """Establish database connection"""
    connection = sqlite3.connect(PATH, timeout=TIMEOUT)
    cursor = connection.cursor()
    return connection, cursor

def parse_feature_source(meta_str: str) -> Tuple[Optional[str], Optional[int]]:
    """
    Parse feature source from metadata JSON
    Returns: (feature_file_path, row_number) or (None, None)
    """
    try:
        meta = json.loads(meta_str) if meta_str else {}
        feature_source = meta.get('feature_source', '')
        if feature_source and ':' in feature_source:
            parts = feature_source.split(':')
            if len(parts) == 2:
                return parts[0], int(parts[1])
    except (json.JSONDecodeError, ValueError):
        pass
    return None, None

def collate_data(connection_s: sqlite3.Connection, 
                 cursor_s: sqlite3.Cursor, 
                 limit: int, 
                 offset: int, 
                 image_path: str, 
                 dataset_store: List[Dict], 
                 running_key_id: int, 
                 resolution_filter: str = '5000',
                 train_mode: bool = False,
                 include_coordinates: bool = True) -> Tuple[List[Dict], int]:
    """
    Enhanced data collation with feature source tracking and coordinate preservation
    
    Args:
        connection_s: SQLite connection
        cursor_s: SQLite cursor
        limit: Number of records to fetch
        offset: Starting position
        image_path: Directory to save images
        dataset_store: List to append records to
        running_key_id: Current running key counter
        resolution_filter: Resolution to filter by
        train_mode: If True, only select unlabeled entries for training
        include_coordinates: If True, include coordinate information
    
    Returns:
        Updated dataset_store and running_key_id
    """
    cmap = LinearSegmentedColormap.from_list("mycmap", colormap())
    
    try:
        cursor_s.row_factory = sqlite3.Row
        
        # Build query based on mode
        if train_mode:
            # Training mode: only get unlabeled entries
            query = """
                SELECT key_id, resolution, viewing_vmax, numpyarr, seqA, seqB, 
                       dimensions, labels, coordinates, meta
                FROM imag_with_seqs 
                WHERE (labels IS NULL OR labels = '' OR labels = 'NONE') 
                  AND resolution = ?
                LIMIT ? OFFSET ?
            """
        else:
            # Inference mode: get labeled entries
            query = """
                SELECT key_id, resolution, viewing_vmax, numpyarr, seqA, seqB, 
                       dimensions, labels, coordinates, meta
                FROM imag_with_seqs 
                WHERE labels IS NOT NULL AND labels != '' AND labels != 'NONE'
                  AND resolution = ?
                LIMIT ? OFFSET ?
            """
        
        cursor_s.execute(query, (resolution_filter, limit, offset))
        
        for en in cursor_s.fetchall():
            try:
                # Verify data integrity
                dimensions = int(en['dimensions'])
                assert len(en['numpyarr']) == 4 * dimensions * dimensions
                
                running_key_id += 1
                
                # Extract numpy array and create image
                arr = np.frombuffer(en['numpyarr'], dtype=np.float32).reshape(dimensions, dimensions)
                arr_norm = (arr - arr.min()) / (en['viewing_vmax'] - arr.min()) if en['viewing_vmax'] > arr.min() else arr
                rgba = cmap(arr_norm)
                rgba = (rgba * 255).astype(np.uint8)
                img = Image.fromarray(rgba, mode="RGBA")
                
                # Save image
                img_name = f"{image_path}/{running_key_id}_image.png"
                img.save(img_name)
                
                # Parse coordinates if available
                coords_dict = {}
                if include_coordinates and en['coordinates']:
                    coords = en['coordinates'].split(',')
                    if len(coords) == 6:
                        coords_dict = {
                            'chr1': coords[0],
                            'start1': int(coords[1]),
                            'end1': int(coords[2]),
                            'chr2': coords[3],
                            'start2': int(coords[4]),
                            'end2': int(coords[5])
                        }
                
                # Parse feature source from metadata
                feature_path, row_num = parse_feature_source(en['meta'])
                
                # Build record
                record = {
                    "key_id": running_key_id,
                    "original_key_id": en['key_id'],
                    "resolution": en['resolution'],
                    "image": img_name,
                    "sequenceA": en['seqA'] if en['seqA'] else "",
                    "sequenceB": en['seqB'] if en['seqB'] else "",
                    "labels": en['labels'] if en['labels'] else "NONE",
                    "feature_source": {
                        "file": feature_path,
                        "row": row_num,
                        "full_path": f"{feature_path}:{row_num}" if feature_path and row_num else None
                    }
                }
                
                # Add coordinates if requested
                if include_coordinates and coords_dict:
                    record["coordinates"] = coords_dict
                
                dataset_store.append(record)
                
            except Exception as error:
                print(f"Error processing record {en['key_id']}: {error}")
                continue
                
    except Exception as error:
        print(f"Error in collate_data: {error}")
    
    return dataset_store, running_key_id

def process_databases(db_sources: List[Tuple[str, int]], 
                     output_json: str,
                     image_dir: str,
                     resolution_filter: str = '5000',
                     batch_size: int = 1024,
                     train_mode: bool = False,
                     include_coordinates: bool = True) -> None:
    """
    Process multiple databases and create JSON output
    
    Args:
        db_sources: List of (database_path, limit) tuples
        output_json: Path to output JSON file
        image_dir: Directory to save images
        resolution_filter: Resolution to filter by
        batch_size: Number of records per batch
        train_mode: If True, only select unlabeled entries
        include_coordinates: If True, include coordinate information
    """
    # Create image directory if needed
    Path(image_dir).mkdir(parents=True, exist_ok=True)
    
    dataset_store = []
    running_key_id = 0
    
    # Process each database
    for db_index, (db_path, db_limit) in enumerate(db_sources):
        try:
            print(f"\n{'='*50}")
            print(f"Processing database {db_index + 1}/{len(db_sources)}")
            print(f"Path: {db_path}")
            print(f"Limit: {db_limit}")
            print(f"Mode: {'TRAINING (unlabeled only)' if train_mode else 'INFERENCE (labeled only)'}")
            print(f"Starting running_key_id: {running_key_id}")
            print(f"{'='*50}")
            
            connection_s, cursor_s = call(db_path, 5)
            offset = 0
            
            while offset < db_limit:
                tnow = datetime.datetime.now()
                
                # Process batch
                dataset_store, running_key_id = collate_data(
                    connection_s=connection_s,
                    cursor_s=cursor_s,
                    limit=batch_size,
                    offset=offset,
                    image_path=image_dir,
                    dataset_store=dataset_store,
                    running_key_id=running_key_id,
                    resolution_filter=resolution_filter,
                    train_mode=train_mode,
                    include_coordinates=include_coordinates
                )
                
                offset += batch_size
                
                print(f"Time: {datetime.datetime.now() - tnow}, Current time: {datetime.datetime.now()}")
                print(f"Total records processed: {len(dataset_store)}")
                print(f"Current offset: {offset}/{db_limit}")
                print(f"Current running_key_id: {running_key_id}")
                print()
                
        except Exception as error:
            print(f"Error processing database {db_path}: {error}")
        finally:
            if 'cursor_s' in locals():
                cursor_s.close()
            if 'connection_s' in locals():
                connection_s.close()
            print(f"Closed database: {db_path}")
    
    # Save to JSON
    print(f"\n{'='*50}")
    print(f"Saving {len(dataset_store)} total records to JSON...")
    with open(output_json, "w") as file:
        json.dump(dataset_store, file, indent=2)
    print(f"Data saved to: {output_json}")
    print(f"Final running_key_id: {running_key_id}")
    print(f"{'='*50}")

def main():
    parser = argparse.ArgumentParser(description='Convert SQLite database to JSON for ML pipeline')
    
    # Required arguments
    parser.add_argument('databases', nargs='+', help='Database paths (format: path:limit)')
    parser.add_argument('--output', '-o', required=True, help='Output JSON path')
    parser.add_argument('--image-dir', '-i', required=True, help='Directory for images')
    
    # Optional arguments
    parser.add_argument('--resolution', '-r', default='5000', help='Resolution filter (default: 5000)')
    parser.add_argument('--batch-size', '-b', type=int, default=1024, help='Batch size (default: 1024)')
    parser.add_argument('--train', action='store_true', help='Training mode: only select unlabeled entries')
    parser.add_argument('--no-coords', action='store_true', help='Exclude coordinate information')
    
    args = parser.parse_args()
    
    # Parse database sources
    db_sources = []
    for db_spec in args.databases:
        if ':' in db_spec:
            path, limit = db_spec.rsplit(':', 1)
            db_sources.append((path, int(limit)))
        else:
            # If no limit specified, use a large number
            db_sources.append((db_spec, 999999))
    
    # Process databases
    process_databases(
        db_sources=db_sources,
        output_json=args.output,
        image_dir=args.image_dir,
        resolution_filter=args.resolution,
        batch_size=args.batch_size,
        train_mode=args.train,
        include_coordinates=not args.no_coords
    )

if __name__ == "__main__":
    now = datetime.datetime.now()
    main()
    print(f"\nTotal execution time: {datetime.datetime.now() - now}")