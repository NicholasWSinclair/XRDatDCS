# XRDatDCS

## Overview
**XRDatDCS** is a tool for generating powder X-ray Diffraction (XRD) patterns based on spectra available at the Dynamic Compression Sector (DCS). The goal of this application is to provide a planning tool for users to decide on appropriate sample thicknesses, x-ray energies, detector distances, etc. It allows users to simulate XRD images for ideal powder samples, considering various experimental parameters, detector configurations, and material properties. While most real materials are not ideal powders, the predicted intensities are reasonable estimates for experiment planning purposes, unless your target is composed of very large crystals or very highly oriented crystals. For example, the intensity from rolled metal foils is very textured, but the azimuthally integrated intensity from this program is usually within a factor of 2-3 of each peak (some will be higher than predicted, some lower). 


## Key Features

### 1. Simulation & Analysis
   ** Powder XRD Simulation: Generates 2D XRD patterns for ideal powders.
   
   ** Calculated Spectra: Using srwpy.srwlib 
   
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
    *   **Calculate Spectrum**: Built-in tool to generate spectra based on experimental parameters.

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

### Installation
Create the environment:

> conda env create --file environment.yml      (create environment from yml)

> conda activate env_pyXRD0      				(activate it)

> python .\XrayScatteringGUI.py   			(run the image simulation GUI)


to run this later: 
On Windows: Open anaconda powershell. Otherwise, open a terminal.
> conda activate env_pyXRD0
> cd <install folder>
> python .\XrayScatteringGUI.py

On windows, if you want to make a shortcut on your desktop for this: 
1) find your environment's exe file:
   >> python -c "import sys; print(sys.executable)"
2) make a shortcut to that exe file with the .py file passed as an argument.
   e.g. make a shortcut to: C:\Users\username\AppData\Local\miniconda3\envs\env_pyXRD0\python.exe "C:\YourDocuments\code\XRDatDCS\XrayScatteringGUI.py"
There's an .ico file for the program in the \ui folder.


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
