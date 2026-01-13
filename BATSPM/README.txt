# **BATSPM**: Small Mammal MRI Analysis Toolbox

**Version:** 1.1.1
**Default Template:** Fruitbat
**Adaptability:** Compatible with other small mammals (Mice, Rats, other Bats) if reference templates are provided.
**Dependencies:** MATLAB, SPM25 (should be compatible with spm12 aswell)
**Input Format:** NIfTI (`.nii`)

---

**## 0. Introduction & Capabilities**

**### The Challenge**
Processing MRI data from small animal scanners presents unique challenges compared to human neuroimaging:
* **Orientation Issues:**  Animal scans often suffer from "Pencil Effects" (swapped dimensions) or incorrect origin points due to varying bed positioning.
* **Noise & Ghosting:**  High-field scanners frequently produce significant background noise and artifacts that interfere with standard segmentation algorithms.
* **Manual Labor:**  Standard SPM pipelines usually require manual re-orientation, cropping, and origin-setting for every single subject.

**### The Solution**
**BATSPM** is a fully automated preprocessing and morphometry pipeline designed to standardize these images before running a complete Voxel-Based Morphometry (VBM) workflow. While optimized for bat brains, the underlying algorithms for geometry repair and noise removal should be applicable to most small mammal neuroimaging datasets.

**Key Capabilities:**
* **Input Handling:**  Processes standard `.nii` files (raw Bruker/DICOM data must be converted to NIfTI first).
* **Auto-Repair & Reslicing:**  Automatically detects swapped header dimensions and reslices images to isotropic resolution (standardized voxel sizes) to ensure correct anatomical orientation.
* **Auto-Crop:**  Uses cluster analysis to isolate the brain and remove background "ghosts" without manual intervention.
* **Dual Analysis:**  Supports both Whole-Brain VBM (Statistical Mapping) and ROI-based Volumetry (Atlas Statistics).

---

**## 1. Installation**

**### Setup**
**1.  Download:** Unzip the `BATSPM` folder.
**2.  MATLAB Path:**
    * Open MATLAB.
    * Right-click the `BATSPM` folder -> **Add to Path** -> **Selected Folders and Subfolders**.
    * Run `savepath` in the Command Window to make this permanent


**### Verification:** Is my installation correct?
Before running, ensure your **BATSPM** folder contains exactly these files:


BATSPM/
 ├── batspm_run.m                    (Main Pipeline)
 ├── batspm_batch_run.m         (Batch Processor)
 ├── batspm_atlas_stats.m        (ROI Statistics Tool)
 ├── batspm_save_snapshot.m (QC Generator)
 ├── batspm_reslice_iso.m        (Geometry Fixer)
 ├── Batlas.nii                            (The Reference Atlas)
 ├── Bat_tpm.nii                         (Tissue Probability Map)
 └── README.md                      (This file)
---

**## 2. How to Use**

**### Preparation**
**1.  Create a Results Folder:**  You must manually create an empty folder where the output will be stored.
**2.  Prepare Input Folder:**  Place all your subject `.nii` files into a single Input Folder.
**3.  Check File Names:**  MATLAB struggles with spaces and special symbols.
    * **[X] Bad:** `Bat 1 (final).nii`, `bat-scan+date.nii`
    * **[V] Good:** `Bat01_final.nii`, `T2_Bat05.nii`

 **### A. Batch Processing (Recommended)**
Use this mode to process an entire folder of subjects automatically.

**/matlab**
% 1. Define your folders (Drag and drop folder into MATLAB to get path)
% Ensure these folders exist before running!
input_folder  = '/Users/lab/Data/Raw_Cohort_A';
output_folder = '/Users/lab/Data/Results_Cohort_A';

% 2. Run Batch
% The 'true' flag calculates Atlas Statistics for every subject.
batspm_batch_run(input_folder, output_folder, true);
  **### B. Single Subject Processing**
Use this mode for testing or debugging a specific file 
**/matlab**
raw_file = '/Users/lab/Data/Raw/Bat01.nii';
out_dir  = '/Users/lab/Data/Results';

batspm_run(raw_file, out_dir, true);
 **## 3. Understanding the Output** For every subject processed, a folder is created containing these files:  **### Primary Results** (For Analysis)


----------------------------------------------------------------------------------------
File Name             | Description              | Use Case
----------------------------------------------------------------------------------------

Stats_SubjectID.csv   | Regional Statistics      | Excel table listing volume change (%)
                      |                          | for every Atlas region.

j_y_SubjectID.nii     | Jacobian Determinant     | The “Growth Map”. Use this for group
                      |                          | VBM statistics (T-tests in SPM).

wT2SubjectID.nii      | Normalized Anatomy       | The clean, standard brain image. Use
                      |                          | for creating group-average templates.

QC_SubjectID.jpg     | Quality Control Snapshot | CHECK FIRST.
                      |                          | Left = Gray Anatomy
                      |                          | Right = Colorful Growth Map
----------------------------------------------------------------------------------------

Color Scale Notes:
- Blue  = −30% volume
- Red   = +30% volume
- Gray/Green ≈ <5% difference (healthy / expected)


**### Secondary Files** (Group Quality checks)

---------------------------------------------------
File Name          | Description
---------------------------------------------------

c1SubjectID.nii    | Gray Matter Mask
                   | Cortical segmentation

c2SubjectID.nii    | White Matter Mask
                   | Deep brain segmentation

c3SubjectID.nii    | CSF Mask
                   | Ventricular / CSF segmentation
---------------------------------------------------




**## 4. The Atlas**

Regions 1–25 = Left side, 26–50 = Right side

**### Left Hemisphere**
--------------------------------------------------------
Region | Structure
--------------------------------------------------------

  1    | Cingulate cortex
  2    | Entorhinal cortex
  3    | Temporal cortex
  4    | Auditory cortex
  5    | Pyriform cortex
  6    | Insular cortex
  7    | Posterior parietal cortex
  8    | Primary somatosensory cortex
  9    | Primary visual cortex
 10    | Lateral visual cortex
 11    | Medial visual cortex
 12    | Primary motor cortex
 13    | Limbic frontal cortex
 14    | Orbitofrontal cortex
 15    | Frontal cortex
 16    | Hippocampus
 17    | Amygdala
 18    | Striatum
 19    | Hypothalamus
 20    | Thalamus
 21    | Cerebellum
 22    | Brain stem
 23    | White matter
 24    | Chambers and fluids   (CSF)
 25    | Chambers and fluids   (CSF)
--------------------------------------------------------

**### Right hemisphere**
Regions 26–50 are mirrored versions of 1–25.



**## 5. TROUBLESHOOTING

------------------------------------------------------------

ISSUE:
"Undefined function 'spm_vol'..."

CAUSE:
MATLAB has lost the path to the SPM software.

FIX:
Run the following in MATLAB:

    addpath('/path/to/your/spm')
    savepath


------------------------------------------------------------

ISSUE:
"File not found" errors

CAUSE:
Usually caused by spaces in filenames or missing folders.

FIX:
Rename files to use underscores (e.g., Bat_01.nii)
and ensure your Output Folder exists before starting.

------------------------------------------------------------

### Bonus Tool: Bruker to NIfTI Converter (bruker_to_nii.py)

**Description:**
This is a standalone Python utility included in the toolbox to convert raw Bruker MRI data (2dseq) into SPM-compatible NIfTI files (.nii.gz). It is specifically engineered to handle "bad headers" where the metadata (visu_pars) conflicts with the actual binary data size—a common issue in animal scanners.

**Features:**
* **Automatic Affine Generation:** Calculates the correct orientation matrix (s-form) from the method/reco files, preventing the "swapped dimensions" issue common in raw exports.
* **Header Correction:** Includes a fallback logic that mathematically deduces the correct slice count if the Bruker header file is corrupted or incorrect.
* **Batch Processing:** Converts an entire directory of scans in one go.

**Requirements:**
* Python 3.x
* Libraries: `nibabel`, `numpy` (Install via: `pip install nibabel numpy`)
* Optional: `dcm2niix` (The script tries to use this first, then falls back to its own internal Python converter if dcm2niix fails).

**How to Use:**
1.  Place the script `bruker_to_nii.py` in your main project folder.
2.  Create a folder named exactly `T2 scans` in the same directory.
3.  Drag your raw Bruker subject folders (containing `pdata`, `method`, etc.) into `T2 scans`.
4.  Run the script from your terminal:
    python3 bruker_to_nii.py
5.  Converted files will appear in a new folder called `nii_output`.


