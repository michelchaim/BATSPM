#!/usr/bin/env python3
"""
Convert Bruker MRI images to NIfTI (.nii) format
Enhanced for T2 fruitbat scan analysis pipeline with proper SPM compatibility
"""

import os
import sys
import subprocess
from pathlib import Path
import re
import numpy as np
import nibabel as nib

# --------------------------------------------------------------------------
# PARSING FUNCTIONS
# --------------------------------------------------------------------------

def read_bruker_visu_pars(visu_pars_file):
    params = {}
    try:
        with open(visu_pars_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        size_match = re.search(r'##\$VisuCoreSize=\(\s*\d+\s*\)\s*([\d\s]+)', content)
        if size_match:
            sizes = [int(x) for x in size_match.group(1).strip().split()]
            params['VisuCoreSize'] = sizes
        
        frame_match = re.search(r'##\$VisuCoreFrameCount=(\d+)', content)
        if frame_match:
            params['VisuCoreFrameCount'] = int(frame_match.group(1))
        
        orient_match = re.search(r'##\$VisuCoreOrientation=\(\s*(\d+),\s*(\d+)\s*\)\s*([\d.\s\-]+)', content)
        if orient_match:
            try:
                vals_str = orient_match.group(3).strip()
                vals = [float(x) for x in vals_str.split()]
                if len(vals) >= 9:
                    params['VisuCoreOrientation'] = np.array(vals[:9]).reshape(3, 3)
            except: pass
        
        pos_match = re.search(r'##\$VisuCorePosition=\(\s*(\d+),\s*(\d+)\s*\)\s*([\d.\s\-]+)', content)
        if pos_match:
            try:
                vals_str = pos_match.group(3).strip()
                vals = [float(x) for x in vals_str.split()]
                if len(vals) >= 3:
                    params['VisuCorePosition'] = np.array(vals).reshape(-1, 3)
            except: pass
        
        thick_match = re.search(r'##\$VisuCoreFrameThickness=\(\s*\d+\s*\)\s*([\d.\s]+)', content)
        if thick_match:
            thick = [float(x) for x in thick_match.group(1).strip().split()]
            params['VisuCoreFrameThickness'] = thick[0] if thick else 0.8
            
    except Exception as e:
        print(f"  Warning: Could not fully parse visu_pars: {e}")
    return params

def read_bruker_reco(reco_file):
    params = {}
    try:
        with open(reco_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        size_match = re.search(r'##\$RECO_size=\(\s*\d+\s*\)\s*([\d\s]+)', content)
        if size_match:
            sizes = [int(x) for x in size_match.group(1).strip().split()]
            params['RECO_size'] = sizes
        
        wordtype_match = re.search(r'##\$RECO_wordtype=(.+)', content)
        if wordtype_match:
            params['RECO_wordtype'] = wordtype_match.group(1).strip()
            
    except Exception as e:
        print(f"  Warning: Could not parse reco: {e}")
    return params

def read_bruker_method(parent_dir):
    params = {}
    try:
        method_file = Path(parent_dir) / "method"
        if not method_file.exists(): return params  
        with open(method_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        fov_match = re.search(r'##\$PVM_Fov=\(\s*\d+\s*\)\s*([\d.\s]+)', content)
        if fov_match:
            params['PVM_Fov'] = [float(x) for x in fov_match.group(1).strip().split()]
        
        thick_match = re.search(r'##\$PVM_SliceThickness=(.+)', content)
        if thick_match:
            params['PVM_SliceThickness'] = float(thick_match.group(1).strip())
            
    except Exception as e:
        print(f"  Warning: Could not parse method: {e}")
    return params

def build_affine_from_bruker(visu_params, reco_params, method_params):
    affine = np.eye(4)
    reco_size = reco_params.get('RECO_size', [192, 160])
    
    voxel_x, voxel_y = 1.0, 1.0
    voxel_z = 0.8
    
    if 'PVM_Fov' in method_params and len(method_params['PVM_Fov']) >= 2:
        fov = method_params['PVM_Fov']
        voxel_x = fov[0] / reco_size[1]
        voxel_y = fov[1] / reco_size[0] if len(fov) > 1 else voxel_x
    
    if 'PVM_SliceThickness' in method_params:
        voxel_z = method_params['PVM_SliceThickness']
    elif 'VisuCoreFrameThickness' in visu_params:
        voxel_z = visu_params['VisuCoreFrameThickness']
    
    affine[0, 0] = voxel_x
    affine[1, 1] = voxel_y
    affine[2, 2] = voxel_z
    
    if 'VisuCoreOrientation' in visu_params:
        rotation = visu_params['VisuCoreOrientation']
        affine[:3, :3] = rotation @ np.diag([voxel_x, voxel_y, voxel_z])
    
    if 'VisuCorePosition' in visu_params and visu_params['VisuCorePosition'].shape[0] > 0:
        affine[:3, 3] = visu_params['VisuCorePosition'][0]
    
    return affine

# --------------------------------------------------------------------------
# CORE CONVERSION LOGIC
# --------------------------------------------------------------------------

def read_bruker_2dseq(pdata_dir):
    pdata_dir = Path(pdata_dir)
    visu_params = read_bruker_visu_pars(pdata_dir / "visu_pars")
    reco_params = read_bruker_reco(pdata_dir / "reco")
    seqfile = pdata_dir / "2dseq"
    
    if not seqfile.exists(): raise FileNotFoundError(f"2dseq file not found in {pdata_dir}")
    
    try:
        dtype = np.int16
        if 'RECO_wordtype' in reco_params:
            if '32BIT' in reco_params['RECO_wordtype']: dtype = np.int32
            elif '_8BIT' in reco_params['RECO_wordtype']: dtype = np.int8
        
        with open(seqfile, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=dtype)
        
        print(f"  Data type: {dtype.__name__}, Total elements: {len(data)}")
        
        # --- PRIORITY CHECK: KNOWN "WEIRD" SIZES FIRST ---
        # Bypasses header for files with known slice counts based on file size
        
        if len(data) == 552960:
            print(f"  -> MATCH: 18 slices detected (Header Correction)")
            return data.reshape((18, 192, 160))

        if len(data) == 276480:
            print(f"  -> MATCH: 9 slices detected")
            return data.reshape((9, 192, 160))
            
        if len(data) == 1013760:
            print(f"  -> MATCH: 33 slices detected")
            return data.reshape((33, 192, 160))
            
        # --- STANDARD CHECK ---
        # Uses header info if size matches expectation
        if 'RECO_size' in reco_params:
            reco_size = reco_params['RECO_size']
            frame_count = visu_params.get('VisuCoreFrameCount', 1)
            expected_size = int(np.prod(reco_size)) * frame_count
            
            if len(data) == expected_size:
                return data.reshape((frame_count, reco_size[0], reco_size[1]))
        
        print(f"  Warning: Could not determine exact shape. Data size: {len(data)}")
        return data
        
    except Exception as e:
        print(f"  Error reading 2dseq: {e}")
        raise

def convert_bruker_to_nii(bruker_scan_dir, output_dir=None):
    bruker_scan_dir = Path(bruker_scan_dir)
    if not bruker_scan_dir.exists(): return None
    
    if output_dir is None: output_dir = bruker_scan_dir.parent
    else: output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Converting Bruker scan: {bruker_scan_dir.name}")
    
    # Try dcm2niix (optional)
    try:
        subprocess.run(["dcm2niix", "-f", bruker_scan_dir.name, "-o", str(output_dir), str(bruker_scan_dir)], 
                      capture_output=True, timeout=5)
        nii_files = list(output_dir.glob(f"{bruker_scan_dir.name}*.nii*"))
        if nii_files: 
            print("✓ Converted with dcm2niix")
            return nii_files[0]
    except: pass
    
    # Python Conversion
    try:
        print("Using Python-based conversion...")
        pdata_1 = bruker_scan_dir / "pdata" / "1"
        if not pdata_1.exists(): return None
        
        img_data = read_bruker_2dseq(pdata_1)
        reco_params = read_bruker_reco(pdata_1 / "reco")
        visu_params = read_bruker_visu_pars(pdata_1 / "visu_pars")
        method_params = read_bruker_method(bruker_scan_dir)
        
        affine = build_affine_from_bruker(visu_params, reco_params, method_params)
        nii_img = nib.Nifti1Image(img_data, affine)
        nii_img.header['descrip'] = b'Bruker T2 MRI converted to NIfTI'
        
        output_file = output_dir / f"{bruker_scan_dir.name}.nii.gz"
        nib.save(nii_img, str(output_file))
        print(f"✓ Successfully converted to: {output_file}\n")
        return output_file
            
    except Exception as e:
        print(f"Error: {e}\n")
        return None

def batch_convert(parent_dir, output_dir=None):
    parent_dir = Path(parent_dir)
    if not parent_dir.exists():
        print(f"Error: Directory not found: {parent_dir}")
        return
        
    output_dir = Path(output_dir) if output_dir else parent_dir / "nii_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    scans = [x for x in parent_dir.iterdir() if x.is_dir() and (x/"pdata").exists()]
    print(f"Found {len(scans)} scans. Saving to: {output_dir}\n")
    
    for scan in sorted(scans):
        convert_bruker_to_nii(scan, output_dir)
    print("Batch conversion complete.")

if __name__ == "__main__":
    # Default Paths
    in_folder = "/Users/micheljanzen/Desktop/BATSPM/T2 scans"
    out_folder = "/Users/micheljanzen/Desktop/BATSPM/nii_output"
    
    if len(sys.argv) > 1: convert_bruker_to_nii(sys.argv[1], out_folder)
    else: batch_convert(in_folder, out_folder)
