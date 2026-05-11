# 3D Nucleus Morphology Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A Python tool to analyze nuclear shapes from 3D microscopy images.**

Takes segmented images (each nucleus = a number) and calculates:
- **Volume** (picoliters)
- **Mean curvature** (how irregular the surface is)
- **Positive/negative surface percentage** (bulging vs indentation)
- **Comparison to perfect sphere**

Outputs Excel files with statistics and individual measurements.


## How It Works

**Step 1:** Extract each nucleus from the segmentation mask  
**Step 2:** Create a 3D mesh surface  
**Step 3:** Calculate curvature at every point using the Meyer et al. (2003) algorithm  
**Step 4:** Remove top/bottom 7% (avoids segmentation artifacts)  
**Step 5:** Calculate area-weighted mean curvature (larger surface areas count more)  
**Step 6:** Export everything to Excel

**Key details:**
- Curvature radius: 1.0 µm (adjustable for larger/smaller nuclei)
- Area-weighting prevents vertex density bias


## Quick Start (5 minutes)

**1. Install packages:**
```bash
pip install numpy scipy imageio vedo trimesh pandas openpyxl
```
**2.Edit `src/config.py`: (just change the paths and group names)**
```bash
INPUT_FOLDER = r"C:\your\images\folder"
OUTPUT_EXCEL = r"C:\your\results.xlsx"
GROUPS_ORDERED = ["Control", "Treated"]
GROUP_PATTERNS = {
    "Control": ["control", "DMSO"],
    "Treated": ["treated", "drug"],
}
```
**3. Run**
``` bash 
python run_analysis.py
```


## Configuration

Edit `src/config.py`:

| Setting | What it does | Default |
|---------|--------------|---------|
| `ORIGINAL_SPACING` | Your microscope's pixel size (microns) | (0.258, 0.258, 1.0) |
| `CURVATURE_RADIUS` | Smoothing level - larger = smoother | 1.0 µm |
| `Z_EXCLUSION` | How much of top/bottom to remove | 0.07 (7%) |
| `USE_PARALLEL` | Use multiple CPU cores | True |

**Quick tips:**
- Bigger nuclei? Increase `CURVATURE_RADIUS` to 2.0
- Windows crashing? Set `USE_PARALLEL = False`
- Want different pole exclusion? Change `Z_EXCLUSION` (0.05 = 5%)


## Output Excel File

The Excel file contains:

**1. Summary sheet** - Statistics per group
- Number of nuclei, mean volume, mean curvature, positive/negative surface percentages

**2. [Group]_data sheet** - Individual nuclei
- Volume, mean curvature, positive/negative percentages for each nucleus



# Installation

**Requirements:** Python 3.7+

**From GitHub:**
```bash
git clone https://github.com/YOUR_USERNAME/3D_curvature_analysis.git
cd 3D_curvature_analysis
pip install -r requirements.txt
```

