# pyXRDatDCS - X-ray Scattering GUI

## Overview
**pyXRDatDCS** is a tool for generating powder X-ray Diffraction (XRD) patterns based on spectra available at the Dynamic Compression Sector (DCS). It allows users to simulate XRD images for ideal powder samples, considering various experimental parameters, detector configurations, and material properties.

The goal of this application is to provide a planning tool for users to decide on appropriate sample thicknesses, x-ray energies, detector distances, etc. 

## Key Features

### 1. Simulation & Analysis
   ** Powder XRD Simulation: Generates 2D XRD patterns for ideal powders.
   ** Calculated Spectra: Using srwpy.srwlib or oasys_srw.srwlib
   ** Polar and azimuthal-dependent attenuation from the sample and from Post-Sample Filters
   ** Debye-Waller Factors: Uses tabulated B-factors for elemental crystals (Peng et al.). If not elemental crystal, it can approximate them for CIF files from the materials project (i.e. with elasticity data) or the user can specify the factors for each element. Temperature is assumed to be 300K
   

### 2. Experimental Setup
*   **Detector Selection**: Choose from pre-configured detectors:
    *   Rayonix SX165 (CsI, GOS, Ideal)
    *   Keck PAD (Si, CdTe)
    *   DCS 4-Frame Detectors
    *   Ideal Square Detector
*   **Geometry Configuration**: User adjustable parameters:
    *   Sample-to-detector distance
    *   Beam center and angle
    *   Pixel size and resolution
*   **Post-Sample Filters**: Add and configure multiple filters (material, thickness, density) to simulate attenuation.

### 3. Material & Spectrum Management
*   **Crystal Structure Loading**:
    *   Load `.cif` files directly.
    *   Fetch structures from **Materials Project** (requires API key).
    *   Search and download from the **Crystallography Open Database (COD)**.
*   **Spectrum Handling**:
    *   **Load Spectrum**: Import energy spectra from `.csv` or `.h5` files.
    *   **Calculate Spectrum**: Built-in tools to generate spectra based on experimental parameters.

### 4. Data Management
*   **Save Results**: Export simulated images as **TIFF** files.
*   **H5 Metadata**: Automatically saves a companion `.h5` file containing all simulation parameters (detector, beam, target, filters, spectrum, and full CIF content).
*   **Load Simulation**: Restore the entire simulation state (including crystal structure and spectrum) from a previously saved `.h5` file via `File > Load Simulation`.

## 5. Some Additional Tools
XtalDiffractionViewer.py is a GUI for visualizing single crystal diffraction peaks in a particular crystal orientation, allowing you to select energy ranges to include. 
SingleCrystalPeaks_Manual_WithStrain.py is similarly a GUI for visualizing single xtal peaks. It has a nicer interface for the Laser/Gun geometry. It currently has a feature where if you click on a reflection it will generate a bunch of slightly shifted orientations (not so useful, just a curiosity). You can edit the strain tensor from the menu bar. 


## Installation & Usage

### Prerequisites
*   Conda installation (e.g. miniconda or anaconda)

### Running the Application
To run the application, simply execute the provided batch file:

```batch
pyXRDatDCS.bat
```

This script will:
1.  Map the necessary network drive (`Z:` to `\\DCS100\Internal`) if not already connected.
2.  Set up the environment variables.
3.  Launch the `XrayScatteringGUI.py` application using the configured Python environment (`env_pyXRD1`).

### Running without the .bat (Linux/macOS)
Alternatively, you can just make the environment manually

conda env create --file environment.yml      (create environment from yml)
conda activate env_pyXRD1      				(activate it)
python .\XrayScatteringGUI.py   			(run it)

## Dependencies
*   **GUI**: PyQt5
*   **Computation**: NumPy, SciPy, XrayLib
*   **Plotting**: Matplotlib
*   **Data Handling**: h5py, pandas, PIL
*   **Crystallography**: ASE, PyFAI, Spglib

## Author
**Nick Sinclair** (Washington State University)
*   Email: nicholas.sinclair@wsu.edu
*   Date: Dec 11, 2024
