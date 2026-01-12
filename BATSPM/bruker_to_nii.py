#!/usr/bin/env python3
"""
Convert Bruker MRI images to NIfTI (.nii) format
Enhanced for T2 fruitbat scan analysis pipeline with proper SPM compatibility

Features:
- Reads Bruker raw 2dseq binary data
- Parses visu_pars, reco, and method files for proper spatial information
- Builds correct affine transformation matrices with orientation and position
- Creates SPM-compatible NIfTI headers
- Supports batch conversion of multiple scans
"""

import os
import sys
import subprocess
from pathlib import Path
import re
import numpy as np
import nibabel as nib


def read_bruker_visu_pars(visu_pars_file):
    """
    Parse Bruker visu_pars file to extract image parameters
    
    Returns:
        dict with VisuCoreSize, VisuCoreFrameCount, orientation matrices, etc.
    """
    params = {}
    try:
        with open(visu_pars_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        # Extract VisuCoreSize
        size_match = re.search(r'##\$VisuCoreSize=\(\s*\d+\s*\)\s*([\d\s]+)', content)
        if size_match:
            sizes = [int(x) for x in size_match.group(1).strip().split()]
            params['VisuCoreSize'] = sizes
        
        # Extract VisuCoreFrameCount
        frame_match = re.search(r'##\$VisuCoreFrameCount=(\d+)', content)
        if frame_match:
            params['VisuCoreFrameCount'] = int(frame_match.group(1))
        
        # Extract VisuCoreOrientation (rotation matrix)
        orient_match = re.search(r'##\$VisuCoreOrientation=\(\s*(\d+),\s*(\d+)\s*\)\s*([\d.\s\-]+)', content)
        if orient_match:
            try:
                vals_str = orient_match.group(3).strip()
                vals = [float(x) for x in vals_str.split()]
                if len(vals) >= 9:
                    params['VisuCoreOrientation'] = np.array(vals[:9]).reshape(3, 3)
            except:
                pass
        
        # Extract VisuCorePosition (position of slices)
        pos_match = re.search(r'##\$VisuCorePosition=\(\s*(\d+),\s*(\d+)\s*\)\s*([\d.\s\-]+)', content)
        if pos_match:
            try:
                vals_str = pos_match.group(3).strip()
                vals = [float(x) for x in vals_str.split()]
                if len(vals) >= 3:
                    params['VisuCorePosition'] = np.array(vals).reshape(-1, 3)
            except:
                pass
        
        # Extract VisuCoreExtent (physical extent in mm)
        extent_match = re.search(r'##\$VisuCoreExtent=\(\s*\d+\s*\)\s*([\d.\s]+)', content)
        if extent_match:
            extents = [float(x) for x in extent_match.group(1).strip().split()]
            params['VisuCoreExtent'] = extents
        
        # Extract frame thickness (slice thickness)
        thick_match = re.search(r'##\$VisuCoreFrameThickness=\(\s*\d+\s*\)\s*([\d.\s]+)', content)
        if thick_match:
            thick = [float(x) for x in thick_match.group(1).strip().split()]
            params['VisuCoreFrameThickness'] = thick[0] if thick else 0.8
        
    except Exception as e:
        print(f"  Warning: Could not fully parse visu_pars: {e}")
    
    return params


def read_bruker_reco(reco_file):
    """
    Parse Bruker reco file to extract image parameters
    
    Returns:
        dict with RECO_size and RECO_wordtype
    """
    params = {}
    try:
        with open(reco_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        # Extract RECO_size
        size_match = re.search(r'##\$RECO_size=\(\s*\d+\s*\)\s*([\d\s]+)', content)
        if size_match:
            sizes = [int(x) for x in size_match.group(1).strip().split()]
            params['RECO_size'] = sizes
        
        # Extract RECO_wordtype (data type: _16BIT_SGN_INT, _32BIT_SGN_INT, etc.)
        wordtype_match = re.search(r'##\$RECO_wordtype=(.+)', content)
        if wordtype_match:
            wordtype = wordtype_match.group(1).strip()
            params['RECO_wordtype'] = wordtype
        
    except Exception as e:
        print(f"  Warning: Could not parse reco: {e}")
    
    return params


def read_bruker_method(parent_dir):
    """
    Parse Bruker method file to extract acquisition parameters
    
    Returns:
        dict with FOV, slice thickness, and other parameters
    """
    params = {}
    try:
        method_file = Path(parent_dir) / "method"
        if not method_file.exists():
            return params
            
        with open(method_file, 'r', encoding='latin1', errors='ignore') as f:
            content = f.read()
        
        # Extract PVM_Fov (Field of View in mm)
        fov_match = re.search(r'##\$PVM_Fov=\(\s*\d+\s*\)\s*([\d.\s]+)', content)
        if fov_match:
            fov = [float(x) for x in fov_match.group(1).strip().split()]
            params['PVM_Fov'] = fov
        
        # Extract PVM_SliceThickness
        thick_match = re.search(r'##\$PVM_SliceThickness=(.+)', content)
        if thick_match:
            params['PVM_SliceThickness'] = float(thick_match.group(1).strip())
        
        # Extract PVM_NSlices
        nslices_match = re.search(r'##\$PVM_NSlices=(\d+)', content)
        if nslices_match:
            params['PVM_NSlices'] = int(nslices_match.group(1))
        
    except Exception as e:
        print(f"  Warning: Could not parse method: {e}")
    
    return params


def build_affine_from_bruker(visu_params, reco_params, method_params):
    """
    Build a proper 4x4 affine transformation matrix from Bruker parameters
    
    This matrix maps voxel indices to physical coordinates (mm).
    Essential for SPM and other neuroimaging tools.
    
    Args:
        visu_params: dict from read_bruker_visu_pars
        reco_params: dict from read_bruker_reco
        method_params: dict from read_bruker_method
    
    Returns:
        4x4 numpy affine matrix
    """
    affine = np.eye(4)
    
    # Get image dimensions from reco
    reco_size = reco_params.get('RECO_size', [192, 160])
    
    # Initialize voxel sizes with defaults
    voxel_x = 1.0  # mm
    voxel_y = 1.0  # mm
    voxel_z = 0.8  # mm (default slice thickness)
    
    # Calculate in-plane voxel sizes from FOV and matrix size
    if 'PVM_Fov' in method_params and len(method_params['PVM_Fov']) >= 2:
        fov = method_params['PVM_Fov']
        # voxel_size = FOV / matrix_size
        voxel_x = fov[0] / reco_size[1]
        voxel_y = fov[1] / reco_size[0] if len(fov) > 1 else voxel_x
    
    # Get slice thickness (z-direction voxel size)
    if 'PVM_SliceThickness' in method_params:
        voxel_z = method_params['PVM_SliceThickness']
    elif 'VisuCoreFrameThickness' in visu_params:
        voxel_z = visu_params['VisuCoreFrameThickness']
    
    # Set diagonal scaling (voxel sizes)
    affine[0, 0] = voxel_x
    affine[1, 1] = voxel_y
    affine[2, 2] = voxel_z
    
    # Apply rotation if available from visu_pars
    if 'VisuCoreOrientation' in visu_params:
        rotation = visu_params['VisuCoreOrientation']
        # Combine rotation with voxel scaling
        affine[:3, :3] = rotation @ np.diag([voxel_x, voxel_y, voxel_z])
    
    # Apply position offset (origin in physical space)
    if 'VisuCorePosition' in visu_params and visu_params['VisuCorePosition'].shape[0] > 0:
        position = visu_params['VisuCorePosition'][0]
        affine[:3, 3] = position
    
    return affine


def read_bruker_2dseq(pdata_dir):
    """
    Read Bruker 2dseq binary file and return image data
    
    Args:
        pdata_dir: Path to pdata/1 directory containing 2dseq file
    
    Returns:
        numpy array with image data (shape: frames x height x width)
    """
    pdata_dir = Path(pdata_dir)
    
    # Read visu_pars to get frame count and dimensions
    visu_pars_file = pdata_dir / "visu_pars"
    visu_params = read_bruker_visu_pars(visu_pars_file)
    
    # Read reco file to get image dimensions and data type
    reco_file = pdata_dir / "reco"
    reco_params = read_bruker_reco(reco_file)
    
    # Read 2dseq binary file
    seqfile = pdata_dir / "2dseq"
    if not seqfile.exists():
        raise FileNotFoundError(f"2dseq file not found in {pdata_dir}")
    
    try:
        # Determine data type from RECO_wordtype
        dtype = np.int16  # default: signed 16-bit
        if 'RECO_wordtype' in reco_params:
            wordtype = reco_params['RECO_wordtype']
            if '32BIT' in wordtype:
                dtype = np.int32 if 'SGN' in wordtype else np.uint32
            elif '_16BIT' in wordtype:
                dtype = np.int16 if 'SGN' in wordtype else np.uint16
            elif '_8BIT' in wordtype:
                dtype = np.int8 if 'SGN' in wordtype else np.uint8
        
        # Read binary data
        with open(seqfile, 'rb') as f:
            data = np.frombuffer(f.read(), dtype=dtype)
        
        print(f"  Data type: {dtype.__name__}, Total elements: {len(data)}")
        
        # Get dimensions and reshape
        if 'RECO_size' in reco_params:
            reco_size = reco_params['RECO_size']
            frame_count = visu_params.get('VisuCoreFrameCount', 1)
            
            expected_size = int(np.prod(reco_size)) * frame_count
            
            if len(data) == expected_size:
                # Reshape: (frames, height, width)
                data = data.reshape((frame_count, reco_size[0], reco_size[1]))
                print(f"  Final shape: {data.shape}")
                return data
        
        # Fallback for common shapes
        if len(data) == 276480:
            data = data.reshape((9, 192, 160))
            print(f"  Using fallback shape: {data.shape}")
            return data
        
        print(f"  Warning: Could not determine exact shape. Data size: {len(data)}")
        return data
        
    except Exception as e:
        print(f"  Error reading 2dseq: {e}")
        raise


def convert_bruker_to_nii(bruker_scan_dir, output_dir=None):
    """
    Convert a single Bruker MRI scan to NIfTI format with proper headers
    
    Args:
        bruker_scan_dir: Path to the Bruker scan directory (parent folder containing pdata)
        output_dir: Output directory for NII files (default: same as input)
    
    Returns:
        Path to the created NII file
    """
    bruker_scan_dir = Path(bruker_scan_dir)
    
    if not bruker_scan_dir.exists():
        print(f"Error: Path does not exist: {bruker_scan_dir}")
        return None
    
    # Set output directory
    if output_dir is None:
        output_dir = bruker_scan_dir.parent
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        print(f"Converting Bruker scan: {bruker_scan_dir.name}")
        
        # Try dcm2niix first if available
        try:
            result = subprocess.run(
                ["which", "dcm2niix"],
                capture_output=True,
                timeout=5
            )
            if result.returncode == 0:
                cmd = [
                    "dcm2niix",
                    "-f", bruker_scan_dir.name,
                    "-o", str(output_dir),
                    str(bruker_scan_dir)
                ]
                print(f"Using dcm2niix: {' '.join(cmd)}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                
                if result.returncode == 0:
                    print(f"✓ Successfully converted with dcm2niix")
                    nii_files = list(output_dir.glob(f"{bruker_scan_dir.name}*.nii*"))
                    return nii_files[0] if nii_files else None
        except:
            pass
        
        # Fallback: Pure Python conversion
        print("Using Python-based conversion...")
        pdata_1 = bruker_scan_dir / "pdata" / "1"
        if not pdata_1.exists():
            print(f"Error: pdata/1 not found in {bruker_scan_dir}")
            return None
        
        # Read the 2dseq data
        img_data = read_bruker_2dseq(pdata_1)
        
        # Read acquisition parameters for proper affine
        reco_file = pdata_1 / "reco"
        reco_params = read_bruker_reco(reco_file)
        visu_pars_file = pdata_1 / "visu_pars"
        visu_params = read_bruker_visu_pars(visu_pars_file)
        method_params = read_bruker_method(bruker_scan_dir)
        
        # Build enhanced affine matrix with orientation and position
        affine = build_affine_from_bruker(visu_params, reco_params, method_params)
        
        # Create NIfTI image with proper affine
        nii_img = nib.Nifti1Image(img_data, affine)
        
        # Set header information for SPM compatibility
        nii_img.header.set_data_dtype(img_data.dtype)
        nii_img.header.set_slope_inter(1.0, 0.0)  # Linear scaling: y = 1.0*x + 0.0
        nii_img.header['descrip'] = b'Bruker T2 MRI converted to NIfTI'
        
        # Print affine information for debugging
        print(f"  Affine matrix (voxel->physical mapping):")
        print(f"    {affine[0,:]}")
        print(f"    {affine[1,:]}")
        print(f"    {affine[2,:]}")
        print(f"    {affine[3,:]}")
        
        # Save NIfTI file
        output_file = output_dir / f"{bruker_scan_dir.name}.nii.gz"
        nib.save(nii_img, str(output_file))
        
        print(f"✓ Successfully converted to: {output_file}")
        print(f"  Shape: {img_data.shape}, Data type: {img_data.dtype}")
        
        return output_file
            
    except Exception as e:
        print(f"Error converting Bruker file: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


def batch_convert(parent_dir, output_dir=None):
    """
    Convert all Bruker scans in a directory to NIfTI format
    
    Args:
        parent_dir: Parent directory containing Bruker scan folders
        output_dir: Output directory for NII files (default: creates nii_output folder)
    """
    parent_dir = Path(parent_dir)
    
    if not parent_dir.exists():
        print(f"Error: Directory does not exist: {parent_dir}")
        return
    
    if output_dir is None:
        output_dir = parent_dir / "nii_output"
    else:
        output_dir = Path(output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all Bruker scan directories (contain pdata subfolder)
    bruker_scans = []
    for item in parent_dir.iterdir():
        if item.is_dir() and (item / "pdata").exists():
            bruker_scans.append(item)
    
    if not bruker_scans:
        print(f"No Bruker scan directories found in: {parent_dir}")
        return
    
    print(f"Found {len(bruker_scans)} Bruker scan(s)")
    print(f"Output directory: {output_dir}\n")
    
    successful = 0
    for scan_path in sorted(bruker_scans):
        if convert_bruker_to_nii(scan_path, output_dir):
            successful += 1
        print()  # Empty line between scans
    
    print(f"\n{'='*60}")
    print(f"✓ Conversion complete: {successful}/{len(bruker_scans)} successful")
    print(f"{'='*60}")


if __name__ == "__main__":
    import sys
    
    # Configuration
    output_folder = "/Users/micheljanzen/Desktop/BATSPM/nii_output"
    
    if len(sys.argv) > 1:
        # Convert specific scan if path provided as argument
        scan_path = sys.argv[1]
        print(f"Converting single Bruker scan to NIfTI...\n")
        convert_bruker_to_nii(scan_path, output_folder)
    else:
        # Default: Batch convert all scans in T2 scans directory
        parent_directory = "/Users/micheljanzen/Desktop/BATSPM/T2 scans"
        print("="*60)
        print("Batch converting all Bruker scans in T2 scans directory...")
        print("="*60)
        print()
        batch_convert(parent_directory, output_folder)
