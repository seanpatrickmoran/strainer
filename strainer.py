#!/usr/bin/env python3
"""
Unified Hi-C Image Pipeline
"""

import json
import yaml
import sqlite3
import numpy as np
import hicstraw
import cooler
import cooltools
import py2bit
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from skimage.filters import threshold_otsu
from skimage.transform import hough_line, hough_line_peaks
from matplotlib.colors import LinearSegmentedColormap

# Constants
DEFAULT_RESOLUTIONS = [2000, 5000, 10000]
DIMENSION_MAP = {
    2000: 162,
    5000: 64,
    10000: 32
}
REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])

class UnifiedHiCPipeline:
    def __init__(self, config_path: str):
        """Initialize pipeline with configuration file"""
        self.config_path = Path(config_path)
        self.config = self.load_config()
        self.feature_mapping = {}  # Maps key_id to original file:row
        self.current_key_id = 1
        
    def load_config(self) -> dict:
        """Load configuration from JSON or YAML file"""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")
        
        with open(self.config_path, 'r') as f:
            if self.config_path.suffix in ['.yaml', '.yml']:
                return yaml.safe_load(f)
            else:  # Default to JSON
                return json.load(f)
    
    def validate_dimensions(self, resolutions: List[int]) -> Dict[int, int]:
        """Validate and calculate dimensions for given resolutions"""
        dimensions = {}
        
        # Use predefined dimensions or calculate
        for res in resolutions:
            if res in DIMENSION_MAP:
                dimensions[res] = DIMENSION_MAP[res]
            else:
                # Scale proportionally - smaller resolution needs larger dimension
                # 2000bp -> 163px, 5000bp -> 65px, 10000bp -> 33px
                # Approximately: dimension = 326000 / resolution
                dim = int(326000 / res)
                dimensions[res] = dim
            
            # Warn if dimension exceeds 200
            if dimensions[res] > 200:
                response = input(f"Warning: Resolution {res} will create {dimensions[res]}x{dimensions[res]} images. Continue? [Y/n]: ")
                if response.lower() == 'n':
                    raise ValueError("User cancelled due to large image dimensions")
        
        return dimensions
    
    def choose_vmax(self, nimage: np.ndarray) -> Tuple[np.ndarray, float]:
        """Calculate viewing maximum for Hi-C matrix visualization"""
        thresh = threshold_otsu(nimage)
        binary = nimage > thresh
        tested_angles = np.linspace(-np.pi / 2, np.pi / 2, 360, endpoint=False)
        h, theta, d = hough_line(binary, theta=tested_angles)
        
        c_angle, c_dist = 0, 0
        chosen = 10
        
        for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
            x0, y0 = dist * np.array([np.cos(angle), np.sin(angle)])
            slope = np.tan(angle + np.pi / 2)
            
            if abs(1 - abs(slope)) < abs(1 - abs(chosen)):
                c_angle = angle
                c_dist = dist
                chosen = slope
        
        x0, y0 = c_dist * np.array([np.cos(c_angle), np.sin(c_angle)])
        slope = np.tan(c_angle + np.pi / 2)
        n = nimage.shape[0]
        
        if slope > 10:
            return nimage, np.max(nimage)
        
        mimage = nimage.copy()
        if int(np.round(y0 - x0)) <= n // 10:
            mimage = np.triu(mimage, k=n//10+2)
        else:
            mimage = np.triu(mimage, k=-n//10)
        
        return mimage, np.max(mimage) if np.max(mimage) > 0 else 1.0
    
    def windowing(self, x1: int, x2: int, y1: int, y2: int, res: int, width: int) -> Tuple[int, int, int, int]:
        """Calculate window coordinates for feature extraction"""
        x1, x2, y1, y2, res = int(x1), int(x2), int(y1), int(y2), int(res)
        target = res * width
        
        if x2 - x1 >= target and y2 - y1 >= target:
            return x1, x2, y1, y2
        
        cx = x1 + (x2 - x1) // 2
        cy = y1 + (y2 - y1) // 2
        
        r1 = cx - target // 2
        r2 = cx + target // 2
        r3 = cy - target // 2
        r4 = cy + target // 2
        
        # Adjust if out of bounds
        if r1 < 0:
            shift = -r1
            r1, r2 = 0, r2 + shift
        if r3 < 0:
            shift = -r3
            r3, r4 = 0, r4 + shift
        
        return r1, r2, r3, r4
    
    def process_hic_file(self, hic_path: str, feature_path: str, resolution: int, 
                        dimension: int, norm: str = "NONE") -> Dict:
        """Process a single Hi-C file with features"""
        hold_np = {}
        

        # Check files exist
        if not os.path.exists(hic_path):
            print(f"ERROR: Hi-C file not found: {hic_path}")
            return {}
        
        if not os.path.exists(feature_path):
            print(f"ERROR: Feature file not found: {feature_path}")
            return {}
        
        print(f"    Loading {hic_path}")
        print(f"    With features from {feature_path}")

        if hic_path.endswith('.hic'):
            hic = hicstraw.HiCFile(hic_path)
            last_chr = None
            matrix_obj = None
            
            with open(feature_path, 'r') as f:
                # Skip header if exists
                first_line = f.readline()
                if first_line.startswith('#') or 'chr' in first_line.lower():
                    pass  # Header line, skip it
                else:
                    f.seek(0)  # No header, go back to start
                
                for row_num, line in enumerate(f, start=1):
                    try:
                        parts = line.strip().split('\t')
                        if len(parts) < 6:
                            continue
                        
                        c1, x1, x2, c2, y1, y2 = parts[:6]
                        x1, x2, y1, y2 = int(x1), int(x2), int(y1), int(y2)
                        
                        # Strip 'chr' prefix if present for hicstraw
                        c1_clean = c1.lstrip('chr')
                        c2_clean = c2.lstrip('chr')
                        
                        # Load matrix if chromosome changed
                        if c1_clean != last_chr:
                            matrix_obj = hic.getMatrixZoomData(c1_clean, c2_clean, 
                                                              "observed", norm, "BP", resolution)
                            last_chr = c1_clean
                        
                        # Calculate window
                        r1, r2, r3, r4 = self.windowing(x1, x2, y1, y2, resolution, dimension)
                        
                        if r1 < 0 or r3 < 0:
                            continue
                        
                        # Extract matrix
                        np_mat = matrix_obj.getRecordsAsMatrix(r1, r2, r3, r4)
                        
                        if np_mat.shape != (dimension+1, dimension+1):
                            print(f"      Matrix shape: {np_mat.shape}, Expected: ({dimension+1}, {dimension+1})")
                            print(f"      Skipping - shape mismatch")
                            continue
                        

                        # Calculate histograms
                        mum_mat, vmax = self.choose_vmax(np_mat)
                        true_max = np.max(np_mat)
                        
                        hist_rel = np.zeros(256, dtype=np.int32)
                        hist_true = np.zeros(256, dtype=np.int32)
                        
                        for i in range(np_mat.shape[0]):
                            for j in range(np_mat.shape[1]):
                                # rel_val = int(np.clip(np_mat[i,j] / vmax * 255, 0, 255))
                                # true_val = int(np.clip(np_mat[i,j] / true_max * 255, 0, 255))
                                # hist_rel[rel_val] += 1
                                # hist_true[true_val] += 1
                                _n = round(np_mat[i][j]/vmax*255)
                                rel_idx = _n if _n < 255 else 255 if _n < 0 else 0
                                _n = round(np_mat[i][j]/true_max *255)
                                true_idx = _n if _n < 255 else 255 if _n < 0 else 0
                                hist_rel[rel_idx] += 1
                                hist_true[true_idx] += 1





                        
                        # Store result with feature mapping
                        feature_key = f"{feature_path}:{row_num}"
                        self.feature_mapping[self.current_key_id] = feature_key
                        
                        hold_np[self.current_key_id] = {
                            'coordinates': [c1, x1, x2, c2, y1, y2],
                            'numpy_window': np_mat.astype(np.float32),
                            'viewing_vmax': vmax,
                            'true_max': true_max,
                            'hist_rel': hist_rel,
                            'hist_true': hist_true,
                            'feature_source': feature_key
                        }
                        
                        self.current_key_id += 1
                        
                    except Exception as e:
                        print(f"Error processing line {row_num}: {e}")
                        continue
        
        elif hic_path.endswith('.cool') or hic_path.endswith('.mcool'):
            # Handle cooler format files
            if hic_path.endswith('.mcool'):
                # For mcool files, need to specify resolution
                cooler_path = f"{hic_path}::/resolutions/{resolution}"
            else:
                cooler_path = hic_path
            
            try:
                clr = cooler.Cooler(cooler_path)
            except Exception as e:
                print(f"Error loading cooler file: {e}")
                return hold_np
            
            # Get matrix balancing weights if using normalization
            if norm != "NONE" and norm in ["KR", "VC", "VC_SQRT"]:
                balance = True
            else:
                balance = False
            
            last_chr = None
            matrix_data = None
            
            with open(feature_path, 'r') as f:
                # Skip header if exists
                first_line = f.readline()
                if first_line.startswith('#') or 'chr' in first_line.lower():
                    pass  # Header line, skip it
                else:
                    f.seek(0)  # No header, go back to start
                
                for row_num, line in enumerate(f, start=1):
                    try:
                        parts = line.strip().split('\t')
                        if len(parts) < 6:
                            # Handle 3-column format (chr, start, end)
                            if len(parts) >= 3:
                                c1 = parts[0]
                                x1 = int(parts[1])
                                y1 = int(parts[2])
                                x2 = x1 + resolution
                                c2 = c1
                                y2 = y1 + resolution
                            else:
                                continue
                        else:
                            c1, x1, x2, c2, y1, y2 = parts[:6]
                            x1, x2, y1, y2 = int(x1), int(x2), int(y1), int(y2)
                        
                        # Ensure chromosome names match cooler format
                        if not c1.startswith('chr'):
                            c1 = 'chr' + c1
                            c2 = 'chr' + c2
                        
                        # Check if chromosomes exist in cooler file
                        if c1 not in clr.chromnames or c2 not in clr.chromnames:
                            print(f"Warning: Chromosome {c1} or {c2} not found in cooler file")
                            continue
                        
                        # Calculate window coordinates
                        r1, r2, r3, r4 = self.windowing(x1, x2, y1, y2, resolution, dimension)
                        
                        if r1 < 0 or r3 < 0:
                            continue
                        
                        # Check bounds against chromosome sizes
                        chr1_size = clr.chromsizes[c1]
                        chr2_size = clr.chromsizes[c2]
                        
                        if r2 > chr1_size:
                            r2 = chr1_size
                        if r4 > chr2_size:
                            r4 = chr2_size
                        
                        # Fetch matrix data
                        # Need to load new matrix if chromosome changed
                        if c1 != last_chr:
                            if balance:
                                matrix_data = clr.matrix(balance=True).fetch(c1, c2)
                            else:
                                matrix_data = clr.matrix(balance=False).fetch(c1, c2)
                            last_chr = c1
                        
                        # Extract submatrix using bin coordinates
                        bin_size = resolution
                        bin_r1 = r1 // bin_size
                        bin_r2 = r2 // bin_size
                        bin_r3 = r3 // bin_size
                        bin_r4 = r4 // bin_size
                        
                        # Ensure we have the right dimensions
                        expected_bins = dimension+1
                        if (bin_r2 - bin_r1) != expected_bins or (bin_r4 - bin_r3) != expected_bins:
                            # Adjust to get exact dimension
                            bin_r2 = bin_r1 + expected_bins
                            bin_r4 = bin_r3 + expected_bins
                        
                        # Extract submatrix
                        np_mat = matrix_data[bin_r1:bin_r2, bin_r3:bin_r4]
                        
                        # Handle NaN values (common in cooler files)
                        np_mat = np.nan_to_num(np_mat, nan=0.0, posinf=0.0, neginf=0.0)
                        
                        # Ensure correct shape
                        if np_mat.shape != (dimension+1, dimension+1):
                            # Pad if necessary
                            if np_mat.shape[0] < dimension+1 or np_mat.shape[1] < dimension+1:
                                padded = np.zeros((dimension+1, dimension+1), dtype=np.float32)
                                padded[:np_mat.shape[0], :np_mat.shape[1]] = np_mat
                                np_mat = padded
                            else:
                                # Trim if too large
                                np_mat = np_mat[:dimension+1, :dimension+1]
                        
                        # Convert to float32
                        np_mat = np_mat.astype(np.float32)
                        
                        # Calculate viewing parameters
                        mum_mat, vmax = self.choose_vmax(np_mat)
                        true_max = np.max(np_mat) if np.max(np_mat) > 0 else 1.0
                        
                        # Calculate histograms
                        hist_rel = np.zeros(256, dtype=np.int32)
                        hist_true = np.zeros(256, dtype=np.int32)
                        
                        for i in range(np_mat.shape[0]):
                            for j in range(np_mat.shape[1]):
                                if vmax > 0:
                                    rel_val = int(np.clip(np_mat[i,j] / vmax * 255, 0, 255))
                                else:
                                    rel_val = 0
                                
                                if true_max > 0:
                                    true_val = int(np.clip(np_mat[i,j] / true_max * 255, 0, 255))
                                else:
                                    true_val = 0
                                
                                hist_rel[rel_val] += 1
                                hist_true[true_val] += 1
                        
                        # Store result with feature mapping
                        feature_key = f"{feature_path}:{row_num}"
                        self.feature_mapping[self.current_key_id] = feature_key
                        
                        hold_np[self.current_key_id] = {
                            'coordinates': [c1, x1, x2, c2, y1, y2],
                            'numpy_window': np_mat,
                            'viewing_vmax': vmax,
                            'true_max': true_max,
                            'hist_rel': hist_rel,
                            'hist_true': hist_true,
                            'feature_source': feature_key
                        }
                        
                        self.current_key_id += 1
                        
                    except Exception as e:
                        print(f"Error processing line {row_num}: {e}")
                        continue
        
        return hold_np
    
    def extract_sequences(self, coordinates: List, genome: str) -> Tuple[str, str]:
        """Extract DNA sequences for given coordinates"""
        genome_path = self.config.get(f'{genome.upper()}_PATH')
        if not genome_path or not os.path.exists(genome_path):
            return "", ""
        
        tb = py2bit.open(genome_path)
        
        try:
            chr_a, x1, x2, chr_b, y1, y2 = coordinates
            x1, x2, y1, y2 = int(x1), int(x2), int(y1), int(y2)
            
            seq_a = tb.sequence(chr_a, x1, x2)
            seq_b = tb.sequence(chr_b, y1, y2)
            
            return seq_a, seq_b
        except Exception as e:
            print(f"Error extracting sequences: {e}")
            return "", ""
        finally:
            tb.close()
    
    def create_database(self, output_path: str):
        """Create output SQLite database with simplified schema"""
        conn = sqlite3.connect(output_path)
        cursor = conn.cursor()
        
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS imag_with_seqs (
                key_id INTEGER PRIMARY KEY,
                name TEXT,
                dataset TEXT,
                condition TEXT,
                coordinates TEXT,
                numpyarr BLOB,
                viewing_vmax REAL,
                true_max REAL,
                hist_rel BLOB,
                hist_true BLOB,
                dimensions INTEGER,
                hic_path TEXT,
                resolution INTEGER,
                norm TEXT,
                seqA TEXT,
                seqB TEXT,
                toolsource TEXT,
                featuretype TEXT,
                labels TEXT DEFAULT '',
                meta TEXT
            )
        """)
        
        # Create feature mapping table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS feature_mapping (
                key_id INTEGER PRIMARY KEY,
                source_file TEXT,
                source_row INTEGER,
                FOREIGN KEY (key_id) REFERENCES imag_with_seqs(key_id)
            )
        """)
        
        conn.commit()
        conn.close()

    def run(self, output_db: str):
        """Execute the full pipeline with batch writing"""
        print("Starting unified Hi-C image pipeline...")
        print(f"Config has {len(self.config.get('datasets', []))} datasets")
        
        for dataset in self.config.get('datasets', []):
            print(f"  Dataset: {dataset['name']}")
            print(f"    Hi-C: {dataset.get('hic_path', 'MISSING')}")
            print(f"    Features: {dataset.get('feature_path', 'MISSING')}")
        
        # Create output database
        self.create_database(output_db)
        conn = sqlite3.connect(output_db)
        cursor = conn.cursor()
        
        # Batch writing parameters
        batch_size = 1000
        batch_data = []
        feature_mappings = []
        
        try:
            # Process each dataset
            for dataset in self.config.get('datasets', []):
                print(f"\nProcessing dataset: {dataset['name']}")
                
                hic_path = dataset['hic_path']
                feature_path = dataset['feature_path']
                
                # Get resolutions and options
                resolutions = dataset.get('resolutions', DEFAULT_RESOLUTIONS)
                genome = dataset.get('genome', 'hg38')
                norm = dataset.get('options', {}).get('norm', 'NONE')
                toolsource = dataset.get('options', {}).get('toolsource', 'unknown')
                featuretype = dataset.get('options', {}).get('featuretype', 'unknown')
                
                # Validate dimensions
                dimensions = self.validate_dimensions(resolutions)
                
                # Process each resolution
                for resolution in resolutions:
                    print(f"  Processing resolution: {resolution}bp")
                    dimension = dimensions[resolution]
                    
                    # Process Hi-C file
                    image_dict = self.process_hic_file(
                        hic_path, feature_path, resolution, dimension, norm
                    )
                    
                    # Process each image
                    for key_id, data in image_dict.items():
                        coords = data['coordinates']
                        coord_str = f"{coords[0]},{coords[1]},{coords[2]},{coords[3]},{coords[4]},{coords[5]}"
                        
                        # Extract sequences if genome file available
                        seq_a, seq_b = self.extract_sequences(coords, genome)
                        
                        # Prepare metadata
                        meta = json.dumps({
                            'genome': genome,
                            'feature_source': data['feature_source'],
                            'dataset_name': dataset['name']
                        })
                        
                        # Add to batch
                        batch_data.append((
                            key_id,
                            f"{dataset['name']}_{resolution}_{key_id}",
                            dataset.get('dataset', dataset['name']),
                            dataset.get('condition', ''),
                            coord_str,
                            data['numpy_window'].tobytes(),
                            data['viewing_vmax'],
                            data['true_max'],
                            data['hist_rel'].tobytes(),
                            data['hist_true'].tobytes(),
                            dimension,
                            hic_path,
                            resolution,
                            norm,
                            seq_a,
                            seq_b,
                            toolsource,
                            featuretype,
                            '',  # Empty labels column
                            meta
                        ))
                        
                        # Add feature mapping
                        source_parts = data['feature_source'].split(':')
                        feature_mappings.append((
                            key_id, 
                            source_parts[0], 
                            int(source_parts[1])
                        ))
                        
                        # Write batch if full
                        if len(batch_data) >= batch_size:
                            try:
                                cursor.executemany("""
                                    INSERT INTO imag_with_seqs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                """, batch_data)
                                
                                cursor.executemany("""
                                    INSERT INTO feature_mapping VALUES (?, ?, ?)
                                """, feature_mappings)
                                
                                conn.commit()
                                print(f"    Wrote batch of {len(batch_data)} records")
                                
                                batch_data = []
                                feature_mappings = []
                                
                            except Exception as e:
                                print(f"    Error writing batch: {e}")
                                # Try to write individually
                                for record in batch_data:
                                    try:
                                        cursor.execute("""
                                            INSERT INTO imag_with_seqs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                                        """, record)
                                    except Exception as e2:
                                        print(f"      Failed to write record {record[0]}: {e2}")
                                conn.commit()
                                batch_data = []
                                feature_mappings = []
            
            # Write remaining batch
            if batch_data:
                try:
                    cursor.executemany("""
                        INSERT INTO imag_with_seqs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, batch_data)
                    
                    cursor.executemany("""
                        INSERT INTO feature_mapping VALUES (?, ?, ?)
                    """, feature_mappings)
                    
                    conn.commit()
                    print(f"  Wrote final batch of {len(batch_data)} records")
                    
                except Exception as e:
                    print(f"  Error writing final batch: {e}")
                    # Try individual writes as fallback
                    for record in batch_data:
                        try:
                            cursor.execute("""
                                INSERT INTO imag_with_seqs VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            """, record)
                        except Exception as e2:
                            print(f"    Failed to write record {record[0]}: {e2}")
                    conn.commit()
        
        finally:
            conn.close()
        
        # Save feature mapping to separate file
        mapping_path = Path(output_db).with_suffix('.mapping.json')
        with open(mapping_path, 'w') as f:
            json.dump(self.feature_mapping, f, indent=2)
        
        print(f"\nPipeline complete!")
        print(f"Database saved to: {output_db}")
        print(f"Feature mapping saved to: {mapping_path}")
        print(f"Total images processed: {self.current_key_id - 1}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python strainer.py <config.json/yaml> [output.db]")
        sys.exit(1)
    
    config_path = sys.argv[1]
    output_db = sys.argv[2] if len(sys.argv) > 2 else "output.db"
    
    pipeline = UnifiedHiCPipeline(config_path)
    pipeline.run(output_db)


if __name__ == "__main__":
    main()
