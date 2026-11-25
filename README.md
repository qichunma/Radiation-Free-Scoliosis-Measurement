# Radiation-Free Scoliosis Measurement via k-d-Tree Filtering and Three-Stage Screening

This repository contains the official MATLAB implementation of the core algorithms described in the paper:
**"Radiation-Free Scoliosis Measurement via k-d-Tree Filtering and Three-Stage Screening"**
*(Manuscript under review at Biomedical Optics Express)*

## 📖 Overview
This project proposes a radiation-free automated method for Adolescent Idiopathic Scoliosis (AIS) screening using structured light 3D reconstruction. The provided script (`main_algorithm.m`) integrates the complete pipeline:
- **Iterative Bilateral Filtering (k-d-tree accelerated)**: Suppresses surface noise while preserving the geometric feature of the spinal sulcus.
- **Slice-based Contour Extraction**: Histogram binning and breakpoint segmentation to extract valid back contours.
- **First-derivative Zero-crossings**: Robustly identifies spinous process candidates, even for shallow or flat spinal depressions.
- **Three-Stage Spatial Screening**: Y-range filtering, nearest-point merging, and adjacent-slice distance comparison.
- **B-spline Fitting & Cobb Angle Calculation**: Generates smooth coronal curves and computes the Cobb angle using tangent angle extrema with a distance threshold.

## ⚙️ Requirements
- **MATLAB** (R2021a or higher recommended).
- **Toolboxes**:
  - Curve Fitting Toolbox (required for `spapi`, `fnval`).
  - Statistics and Machine Learning Toolbox (required for `createns`, `knnsearch`).
  - Image Processing Toolbox (optional, for `smoothdata`).

## 🚀 How to Run
1. Download or clone this repository.
2. Open `main_algorithm.m` in MATLAB.
3. Ensure the input file path is correct. A de-identified sample point cloud (`sample_data.txt`) is provided in the repository.
4. Run the script. 
   - The script will output the filtered point cloud, candidate spinous processes, B-spline fitted curve, and the top calculated Cobb angles in the command window.
   - Figures will automatically pop up to visualize each processing step (filtering comparison, candidate points, 3D/2D spinal midline fitting, and Cobb angle annotations).

## 📊 Data Availability
Due to strict patient privacy and institutional ethical restrictions, the full raw point cloud dataset (containing identifiable body surface information) cannot be publicly released. However, we have provided a **de-identified sample point cloud** (`sample_data.txt`) to ensure the reproducibility of the core algorithmic results presented in the paper.

## ✉️ Contact
- **Qichun Ma** (First Author)
- **Jinzhuang Xiao** (Corresponding Author, E-mail: robot@hbu.edu.cn)
Hebei University, Baoding, China.
