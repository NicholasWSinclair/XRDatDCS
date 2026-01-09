"""
Tool for Generating powder XRD patterns from spectra available at the Dynamic 
Compression Sector (DCS)

Assumptions: The sample is an ideal powder, with crystal sizes large enough 
to not broaden the spectrum substantially, considering the detector PSF
 ~(>100nm)

Debye-Waller Factors are derived either from Peng et al, for elemental 
materials, or are approximated using a Debye temperature, calculated with the 
formula for cubic crystals (from Warren 'X-ray Diffraction') 
-- For non-elemental, non-cubic crystals, type in a more appropriate value in 
the edit box. 

Peng et al: Acta Cryst. (1996). A52, 456-470
Debye-Waller Factors and Absorptive Scattering Factors of Elemental Crystals
L.M. PENG, G. REN, S. L. DUDAREV, AND M. J. WHELAN


Nick Sinclair 2024 Dec 11
Washington State University
nicholas.sinclair@wsu.edu
"""

import sys
import numpy as np
import csv
from matplotlib import use
use("QtAgg")

from collections import defaultdict
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QVBoxLayout, QLineEdit, QPushButton, QLabel, QWidget, QMenuBar, QAction,
    QGroupBox, QHBoxLayout, QSizePolicy, QSlider, QCheckBox, QTableWidget, QTableWidgetItem, QMessageBox, QHeaderView,
    QTabWidget, QStyle, QDockWidget, QSystemTrayIcon, QMenu
)
import re
from collections import Counter
from math import gcd
from functools import reduce


from emmet.core.summary import HasProps
from PyQt5.QtGui import QPalette, QColor, QIcon
from PyQt5.QtCore import Qt, QTimer, QSettings
try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
except ImportError:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from scipy.integrate import quad
import math
from PIL import Image
import h5py
import concurrent.futures

import matplotlib.pyplot as plt

from scipy.fft import fft2, ifft2
from scipy.signal import find_peaks

from matplotlib.figure import Figure

from ase.io import read
from xoppylib.scattering_functions.xoppy_calc_f0 import xoppy_calc_f0 
import spglib
from PyQt5.QtWidgets import QProgressBar
# from ase.spacegroup import crystal
from xoppylib.scattering_functions.f1f2_calc import f1f2_calc
# from xoppylib.scattering_functions.material_constants_library import MaterialConstantsLibrary
# from dabax.dabax_xraylib import DabaxXraylib
import xraylib
from scipy.interpolate import interp1d
import os 
# from pathlib import Path
dir_path = os.path.dirname(os.path.realpath(__file__))
if dir_path not in sys.path:
    sys.path.append(dir_path)
lib_path = os.path.join(dir_path,"lib")
if lib_path not in sys.path:
    sys.path.append(lib_path)
import pandas as pd
# import multiprocessing
from mp_api.client import MPRester
# API_KEY = '' #For materials project API (specific to user)
API_KEY_LOCATION = os.path.join(dir_path,"API_KEY.txt")
Debye_Waller_elementalcrystalpath=os.path.join(dir_path,"DebyeWallerFactorsB_FromPeng_ElementalCrystals.csv")

#pyFAI
import pyFAI.gui.matplotlib
from pyFAI.calibrant import get_calibrant
from pyFAI.integrator.azimuthal import AzimuthalIntegrator

#for searching COD
import requests
from bs4 import BeautifulSoup
import wget
import urllib.parse


#local resources
import importlib
from lib import DCSSpectrum as DCSSpectrum
from lib import PostSampleFilterGUI_tableversion as PostSampleFilterGUI
from lib import GeometryVisualizer

# importlib.reload(PostSampleFilterGUI)
# importlib.reload(GeometryVisualizer)
from time import time, perf_counter

hc = 12398.4198 # electron volts * angstrom
r0 = 2.81794032E-5 # angstroms
JouleIneV = 6.242e+18 #eV per Joule
seed = 328375183878393377431811483858591401665
# material_constants_library = DabaxXraylib()


class XRayScatteringApp(QMainWindow):
    def __init__(self,ciffile=None):
        super().__init__()
        self.subwindows = []
        self.DWBox = None
        self.ciffiletoload = ciffile
        self.init_settings()
        self.init_ui()
        

    def init_settings(self):
        self.settings = QSettings("WSU", "XrayScatteringGUI")
        self.use_cod_sql = self.settings.value("CIF_Sources/COD_UseSQL", True, type=bool)
        self.api_key_path = self.settings.value("CIF_Sources/MP_API_Key_File_Path", os.path.join(dir_path,"API_KEY.txt"), type=str)
        self.MP_APIkeyValid = self.validate_mp_api_key(self.api_key_path)

    def init_ui(self):
                
        
        #Application Icon        
        app_icon = QIcon(os.path.join(lib_path,'AppIcon.png'))
        self.setWindowIcon(app_icon)
        
        # System Tray Icon
        self.tray_icon = QSystemTrayIcon(self)
        self.tray_icon.setIcon(app_icon)
        
        # Tray Icon Menu
        tray_menu = QMenu()
        
        show_action = QAction("Show", self)
        show_action.triggered.connect(self.show)
        tray_menu.addAction(show_action)
        
        hide_action = QAction("Hide", self)
        hide_action.triggered.connect(self.hide)
        tray_menu.addAction(hide_action)
        
        quit_action = QAction("Exit", self)
        quit_action.triggered.connect(QApplication.instance().quit)
        tray_menu.addAction(quit_action)
        
        self.tray_icon.setContextMenu(tray_menu)
        self.tray_icon.show()
        
        # Connect activation (double click to show)
        self.tray_icon.activated.connect(self.tray_icon_activated)
        
        
        #### File Menu
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu('File')
        select_file_action = QAction('Select CIF File', self)        
        select_file_action.triggered.connect(self.select_file)
        file_menu.addAction(select_file_action)
        
        save_image_action = QAction('Save Image as TIF', self)        
        save_image_action.triggered.connect(self.save_image)
        file_menu.addAction(save_image_action)

        load_sim_action = QAction('Load Simulation From H5', self)
        load_sim_action.triggered.connect(self.LoadSimulation)
        file_menu.addAction(load_sim_action)
        
        
        spec_menu = menu_bar.addMenu('Spectrum')
        
        calcspec_action = QAction('Calc Spectrum', self)
        calcspec_action.triggered.connect(self.CalcSpectrum)
        spec_menu.addAction(calcspec_action)
        
        loadspec_action = QAction('Load Spectrum', self)
        loadspec_action.triggered.connect(self.LoadSpectrum)
        spec_menu.addAction(loadspec_action)
        
        #### CIF Sources Menu
        cif_sources_menu = menu_bar.addMenu('CIF Sources')
        
        # COD Submenu
        cod_menu = cif_sources_menu.addMenu('COD')
        self.use_cod_sql_action = QAction('Use SQL interface (Recommended)', self, checkable=True)
        self.use_cod_sql_action.setChecked(self.use_cod_sql)
        self.use_cod_sql_action.triggered.connect(self.set_cod_sql_preference)
        cod_menu.addAction(self.use_cod_sql_action)
        
        # Materials Project Submenu
        mp_menu = cif_sources_menu.addMenu('Materials Project')
        
        specify_api_key_action = QAction('Specify API Key File', self)
        specify_api_key_action.triggered.connect(self.select_mp_api_key_file)
        mp_menu.addAction(specify_api_key_action)
        
        check_api_key_action = QAction('Check API Key', self)
        check_api_key_action.triggered.connect(self.check_mp_api_key)
        mp_menu.addAction(check_api_key_action)
        
        
        #### Detector Menu
        detector_menu = menu_bar.addMenu('Select Detector')        
  
        def adddetectormenu(detname,detectorstaticmethod):
            make_det_action = QAction(detname, self)
            make_det_action.triggered.connect(lambda: self.makedetectorobj(detectorstaticmethod))
            detector_menu.addAction(make_det_action)
            
            
        
        
        adddetectormenu('CsI Rayonix', 'makeCsIRayonix')
        adddetectormenu('GOS Rayonix', 'makeGOSRayonix')
        adddetectormenu('DCS 4-Frame 75mm', 'makeDCSFourFrameDetector_LSO_75mm')
        adddetectormenu('DCS 4-Frame 120mm', 'makeDCSFourFrameDetector_LSO_120mm')
        adddetectormenu('DCS 4-Frame 150mm', 'makeDCSFourFrameDetector_LSO_150mm')
        adddetectormenu('Ideal Square Detector', 'makeGenericSquareDetector')
        adddetectormenu('Ideal Rayonix', 'makeIdealRayonix')
        adddetectormenu('Si KeckPad', 'makeSiKeckPAD')
        adddetectormenu('CdTe KeckPad', 'makeCdTeKeckPAD')
        
        GeometryVis_Menu = menu_bar.addMenu('Expt. Geometry')
        
        DriverGeometry_menu=GeometryVis_Menu.addMenu('Set Geometry')
        
        def addexpttype(ExptTypeName,ExptStaticMethod):
            make_expt_action = QAction(ExptTypeName, self)
            make_expt_action.triggered.connect(lambda: self.makeexptobj(ExptStaticMethod))
            DriverGeometry_menu.addAction(make_expt_action)
        
        addexpttype('E-station Gun', 'makeImpactExperiment_Estation')
        addexpttype('D-station Gun', 'makeImpactExperiment_Dstation')
        addexpttype('Laser Shock', 'makeLaserShockExperiment')
        # setGeo_GunEstation_action = QAction('Show Expt with Lowest E Peaks', self)
        # setGeo_GunEstation_action.triggered.connect()
        # DriverGeometry_menu.addAction(setGeo_GunEstation_action)
        
        
        showgeoWithPeaks_action = QAction('Show Expt with Lowest E Peaks', self)
        showgeoWithPeaks_action.triggered.connect(lambda: self.showgeometry(1))
        GeometryVis_Menu.addAction(showgeoWithPeaks_action)
        
        showgeo_action = QAction('Show Expt without XRD Peaks', self)
        showgeo_action.triggered.connect(lambda: self.showgeometry(0))
        GeometryVis_Menu.addAction(showgeo_action)
        
        View_Menu = menu_bar.addMenu('View')
        # ShowParams_action= QAction('Show Parameter Window',self)
        # ShowParams_action.triggered.connect(lambda: self.showparamwindow())
        # View_Menu.addAction(ShowParams_action)
        
        # ShowParams_action= QAction('Dock Parameter Window',self)
        # ShowParams_action.triggered.connect(lambda: self.showparamwindow())
        
        
        # Dock_targetProps_action= QAction('Show Target Properties',self)
        # Dock_targetProps_action.triggered.connect(lambda: self.dockGroup(self.TargetGroupBox))
        # View_Menu.addAction(Dock_targetProps_action)
        
        # Dock_DetProps_action= QAction('Show Detector Properties',self)
        # Dock_DetProps_action.triggered.connect(lambda: self.dockGroup(self.DetectorGroupBox))
        # View_Menu.addAction(Dock_DetProps_action)
        
        Tools_Menu = menu_bar.addMenu('Tools')
        
        suggest_thickness_action = QAction('Suggest sample thickness', self)
        suggest_thickness_action.triggered.connect(self.suggest_sample_thickness)
        Tools_Menu.addAction(suggest_thickness_action)
        
        suggest_energy_action = QAction('Suggest xray energy', self)
        suggest_energy_action.triggered.connect(self.suggest_xray_energy)
        Tools_Menu.addAction(suggest_energy_action)
        # View_Menu.addAction(ShowParams_action)
        
        # ShowParams_action= QAction('Dock Parameter Window',self)
        # ShowParams_action.triggered.connect(lambda: self.showparamwindow())
        
        
        # Dock_targetProps_action= QAction('Show Target Properties',self)
        # Dock_targetProps_action.triggered.connect(lambda: self.dockGroup(self.TargetGroupBox))
        # View_Menu.addAction(Dock_targetProps_action)
        
        # Dock_DetProps_action= QAction('Show Detector Properties',self)
        # Dock_DetProps_action.triggered.connect(lambda: self.dockGroup(self.DetectorGroupBox))
        # View_Menu.addAction(Dock_DetProps_action)
        
        
        def addedit(layout=None,label=None,placeholdertext=None,defaulttext=None,unitlabel=None,checkfunction=None):
            hboxlayout = QHBoxLayout()
            layout.addLayout(hboxlayout)
            hboxlayout.addWidget(QLabel(label))
            editbox = QLineEdit(self)
            editbox.setPlaceholderText(placeholdertext)
            editbox.setToolTip(placeholdertext)
            editbox.setText(defaulttext)
            hboxlayout.addWidget(editbox)
            if unitlabel is not None:
                hboxlayout.addWidget(QLabel(unitlabel))
            if checkfunction is not None:
                editbox.editingFinished.connect(lambda: checkfunction(editbox))
            editbox.setAlignment((Qt.AlignVCenter | Qt.AlignRight))
            
            return editbox
        
        def addtextline(layout,label,defaulttext,unitlabel):
            hboxlayout = QHBoxLayout()
            layout.addLayout(hboxlayout)
            hboxlayout.addWidget(QLabel(label))
            outlabel = QLabel(self)
            outlabel.setText(defaulttext)
            hboxlayout.addWidget(outlabel)
            outlabel.setAlignment(Qt.AlignVCenter | Qt.AlignRight)
            if unitlabel is not None:
                hboxlayout.addWidget(QLabel(unitlabel))
            return outlabel
        
        def addcheckbox(layout,label,ischecked,donotshow=False):
            
            checkbox = QCheckBox(self)
            if not donotshow:
                hboxlayout = QHBoxLayout()
                layout.addLayout(hboxlayout)
                hboxlayout.addWidget(QLabel(label))
                if ischecked:
                    checkbox.setCheckState(2)
                else:
                    checkbox.setCheckState(0)
                hboxlayout.addWidget(checkbox)
            else:
                checkbox.hide()
            return checkbox
        
        def addbutton(layout,label,callback):
            hboxlayout = QHBoxLayout()
            layout.addLayout(hboxlayout)
            # hboxlayout.addWidget(QLabel(label))
            button = QPushButton(label)
            button.clicked.connect(callback)

            hboxlayout.addWidget(button)
            return button
        
        # Progress bar
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setAlignment(Qt.AlignCenter)
        self.progress_bar.setValue(0)
        
        self.setWindowTitle("XRD Image Simulation")
        self.setGeometry(100, 50, 800, 1000)
        

        # Central widget
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        layout = QVBoxLayout()
        self.central_widget.setLayout(layout)
        
        # self.file_label = QLabel("CIF File: Not selected", self)
        # layout.addWidget(self.file_label)
        
        
        
                
        self.DetectorGroupBox = QGroupBox('Detector Params')
        self.BeamGroupBox = QGroupBox('Beam Params')
        self.TargetGroupBox = QGroupBox('Target Params')
        
        #extra containers for the groupboxes to dock to
        self.detgroupboxcontainer=QVBoxLayout()
        self.beamgroupboxcontainer=QVBoxLayout()
        self.targetgroupboxcontainer=QVBoxLayout()
        
        self.detgroupboxcontainer.addWidget(self.DetectorGroupBox)
        self.beamgroupboxcontainer.addWidget(self.BeamGroupBox)
        self.targetgroupboxcontainer.addWidget(self.TargetGroupBox)
        
        self.DetectorGroupBox.dockingcontainer = self.detgroupboxcontainer
        self.BeamGroupBox.dockingcontainer = self.beamgroupboxcontainer
        self.TargetGroupBox.dockingcontainer = self.targetgroupboxcontainer
        
        
        SampleGroupBox = QGroupBox()
        self.PostSampleFilters_GroupBox = PostSampleFilterGUI.PostSampleFilterManager(self)
        self.PostSampleFilters_GroupBox.add_filter("LiF", 1, 2.64, 52,0.0)
        self.PostSampleFilters_GroupBox.setMinimumSize(400,50)
        paramboxeslayout = QHBoxLayout()
        # paramboxeslayout.addWidget(self.BeamGroupBox)
        paramboxeslayout.addLayout(self.beamgroupboxcontainer)
        
        # paramboxeslayout.addWidget(self.TargetGroupBox)
        paramboxeslayout.addLayout(self.targetgroupboxcontainer)
        targetvlayout = QVBoxLayout()
        self.TargetGroupBox.setLayout(targetvlayout)
        # self.TargetGroupBox.undockgroupbutton = addbutton(targetvlayout, '', lambda: self.move_groupbox_to_new_window(self.TargetGroupBox))
        # self.TargetGroupBox.undockgroupbutton.setIcon(self.style().standardIcon(QStyle.SP_TitleBarShadeButton))
        hsamplespecifylayout=QHBoxLayout()
        targetvlayout.addLayout(hsamplespecifylayout)
        self.SetSampleButton=addbutton(hsamplespecifylayout,'Select CIF file',self.select_file)
        self.SpecifyChemForm_edit = addedit(hsamplespecifylayout,'Or Type:','Chemical Formula','',None)
        self.SpecifyChemForm_edit.returnPressed.connect(self.findchemical)
        self.matproject_checkbox = QCheckBox('Mat Project')
        self.matproject_checkbox.setChecked(False)
        if not self.MP_APIkeyValid:
            self.matproject_checkbox.setEnabled(False)
            self.matproject_checkbox.setToolTip("Valid API Key not found. Please configure in 'CIF Sources' menu IF you want to search Materials Project CIFs.")
        else:
             self.matproject_checkbox.setToolTip("Fetch from Materials Project IF you want to search Materials Project CIFs.")
        hsamplespecifylayout.addWidget(self.matproject_checkbox)
        self.cod_checkbox = QCheckBox('COD')
        self.cod_checkbox.setChecked(True)
        self.cod_checkbox.setToolTip("Fetch from COD IF you want to search COD CIFs.")
        hsamplespecifylayout.addWidget(self.cod_checkbox)
        targetvlayout.addWidget(SampleGroupBox)
        targetvlayout.addWidget(self.PostSampleFilters_GroupBox)
        self.SampleLayout = QVBoxLayout()
        SampleGroupBox.setLayout(self.SampleLayout)
        
        # paramboxeslayout.addWidget(self.DetectorGroupBox)
        paramboxeslayout.addLayout(self.detgroupboxcontainer)
        BeamLayout = QVBoxLayout()
        self.BeamGroupBox.setLayout(BeamLayout)

        DetectorLayout = QVBoxLayout()
        self.DetectorGroupBox.setLayout(DetectorLayout)
        
        paramwidget = QWidget()
        paramwidget.setLayout(paramboxeslayout)
        
        # layout.addLayout(paramboxeslayout)
        self.paramdock = QDockWidget('Simulation Parameters',self)
        self.paramdock.setWidget(paramwidget)
        View_Menu.addAction(self.paramdock.toggleViewAction())
        # layout.addWidget(self.paramdock)
        self.addDockWidget(Qt.TopDockWidgetArea, self.paramdock, Qt.Horizontal)
        self.paramdock.setAllowedAreas(Qt.AllDockWidgetAreas)
        self.paramdock.setFloating(False)
        
        
        for groupbox in (self.DetectorGroupBox, self.BeamGroupBox, SampleGroupBox, self.TargetGroupBox, self.PostSampleFilters_GroupBox):
            groupbox.setSizePolicy(
                QSizePolicy.Expanding,  # Horizontal expansion
                QSizePolicy.Preferred   # Respect preferred height
            )
        
        

        self.DetectorGroupBox.windowlabel='Detector Parameters'
        self.TargetGroupBox.windowlabel='Target Parameters'
        self.BeamGroupBox.windowlabel='Beam Parameters'
        

        # self.file_button = QPushButton("Select CIF File", self)
        # self.file_button.clicked.connect(self.select_file)
        # layout.addWidget(self.file_button)
        
        
        # show_atom_3d_viewer(atoms, cell)
        
    

        #### Beam Params
        # self.BeamGroupBox.undockgroupbutton = addbutton(BeamLayout, '', lambda: self.move_groupbox_to_new_window(self.BeamGroupBox))
        # self.BeamGroupBox.undockgroupbutton.setIcon(self.style().standardIcon(QStyle.SP_TitleBarShadeButton))
        self.SetSpectrumButton=addbutton(BeamLayout,'Set Spectrum',self.CalcSpectrum)
        BeamLayout.addWidget(QLabel(" "))
        self.polarizedcheck = addcheckbox(BeamLayout,'Polarized Beam',0)
        self.enableFluorescence = addcheckbox(BeamLayout,'EnableFluorescence',0,donotshow=True)
        self.onlyfluorescence = addcheckbox(BeamLayout,'Only Fluorescence',0,donotshow=True)
        
        self.xcenterpix_input = addedit(BeamLayout,'Beam Center Hpix#','H-center pixel # (0 is left-most)','1025','    ')
        self.ycenterpix_input = addedit(BeamLayout,'Beam Center Vpix#','V-center pixel # (0 is left-most)','1025','    ')
        
        BeamLayout.addWidget(QLabel(" "))
        DWGroupbox = QGroupBox('Debye-Waller Config')
        self.DWlayout = QVBoxLayout()
        DWGroupbox.setLayout(self.DWlayout)
        BeamLayout.addWidget(DWGroupbox)
        
        self.DebyeWallerCheckbox =  addcheckbox(self.DWlayout,'Enable Debye Waller',0)
        self.CalcDebyeFactorFromMaterialsProjectButton=addbutton(self.DWlayout,'Get DW From Mat. Proj.',self.CalcDebyeFactorFromMaterialsProject)
        
        BeamLayout.setAlignment(Qt.AlignTop)
        
        #### Detector Params  
        # self.DetectorGroupBox.undockgroupbutton = addbutton(DetectorLayout, '', lambda: self.move_groupbox_to_new_window(self.DetectorGroupBox))
        # self.DetectorGroupBox.undockgroupbutton.setIcon(self.style().standardIcon(QStyle.SP_TitleBarShadeButton))
        self.pixelsize_input = addedit(DetectorLayout,'Pixel Size (µm)','Pixel Size','79','µm',checkfunction=self.isvalidvalue)
        self.SampleToDetDistance_input = addedit(DetectorLayout,'Sample To Det Distance','Det Distance','200','mm',checkfunction=self.isvalidvalue)
        self.MHorizontal_input = addedit(DetectorLayout,'# Horiz. Pix (Half), 2*','M','1024','+1',checkfunction=self.isvalidvalue)
        self.NVertical_input = addedit(DetectorLayout,'# Vert. Pix (Half), 2*','N','1024','+1',checkfunction=self.isvalidvalue)
        self.InstrFWHM_input = addedit(DetectorLayout,'Instr FWHM','Instr. FWHM','0.2','mm',checkfunction=self.isvalidvalue)
        
        self.ScintThickness_input = addedit(DetectorLayout,'Scint. Thickness','Thickness (mm(','0.08','mm',checkfunction=self.isvalidvalue)
        self.readnoiseoffset_input = addedit(DetectorLayout,'Read Noise offset','Read Noise Offset (cts)','11','counts',checkfunction=self.isvalidvalue)
        self.readnoiseSTD_input = addedit(DetectorLayout,'Read Noise STD Dev','Read Noise STD (cts)','2','counts',checkfunction=self.isvalidvalue)
        self.gainADUperkeV = addedit(DetectorLayout,'Gain (ADU/keV)','Gain (ADU/keV)','0.076','ADU/keV',checkfunction=self.isvalidvalue)
        self.SimNoiseCheckbox =  addcheckbox(DetectorLayout,'Simulate Noise',1)
        DetectorLayout.setAlignment(Qt.AlignTop)
        
        
        
        #### Sample Params

        self.chemform_out = addtextline(self.SampleLayout,'Formula:','','')
        self.thickness_input =addedit(self.SampleLayout,'Thickness','Thickness','0.0075','mm',checkfunction=self.isvalidvalue)
        self.sampleangle_input =addedit(self.SampleLayout,'Beam Angle w.r.t. Normal','Beam Angle w.r.t. Normal','52','deg',checkfunction=self.isvalidangle)
        self.density_out =  addtextline(self.SampleLayout,'Density (From CIF):','','g/cc')
        # self.crystalsize_input =addedit(self.SampleLayout,'Typical Crystal Size','Typical Xtal Size','0.01','µm')
        
        self.SampleLayout.setAlignment(Qt.AlignTop)
        
        

        #### Action Buttons
        ButtonGroupBox = QGroupBox()
        layout.addWidget(ButtonGroupBox)
        ButtonLayout = QHBoxLayout()
        ButtonGroupBox.setLayout(ButtonLayout)
        
        # Show Spectrum
        self.ShowSpectrum_button = QPushButton("Show Spectrum", self)
        self.ShowSpectrum_button.clicked.connect(self.showspectrum)
        ButtonLayout.addWidget(self.ShowSpectrum_button)
        
        
        # Calculate button
        self.calc_button = QPushButton("Calculate 1D Pattern (No Atten.)", self)
        self.calc_button.clicked.connect(self.calculate_and_plot_1D)
        ButtonLayout.addWidget(self.calc_button)
        
        # Calculate Image button
        self.calcsingle_button = QPushButton("Calculate Image", self)
        self.calcsingle_button.clicked.connect(self.calculate_and_plot_2D)
        ButtonLayout.addWidget(self.calcsingle_button)
        
        # Display Atom Positions button
        self.displayatoms_button = QPushButton("Display Atom Positions", self)
        self.displayatoms_button.clicked.connect(self.displayatompositions)
        ButtonLayout.addWidget(self.displayatoms_button)
        
        # Display Calc Fluor button
        self.calcfluor_button = QPushButton("Calc Fluorescence", self)
        self.calcfluor_button.clicked.connect(self.calctotalfluorescence)
        ButtonLayout.addWidget(self.calcfluor_button)
        
        # Display Calc QE button
        self.calcQE_button = QPushButton("Plot Detector QE", self)
        self.calcQE_button.clicked.connect(self.plotQE)
        ButtonLayout.addWidget(self.calcQE_button)

        # Output and Plot
        self.result_label = QLabel("", self)
        layout.addWidget(self.result_label)
        
        

        
        # Sliders (XRD Image)
        self.vmin_slider = QSlider(Qt.Horizontal)
        self.vmin_slider.setMinimum(0)
        self.vmin_slider.setMaximum(110)
        self.vmin_slider.setValue(0)
        self.vmin_slider.valueChanged.connect(self.update_image_display)        
        self.vmax_slider = QSlider(Qt.Horizontal)
        self.vmax_slider.setMinimum(0)
        self.vmax_slider.setMaximum(110)        
        self.vmax_slider.setValue(100)
        self.vmax_slider.valueChanged.connect(self.update_image_display)
        
        # Sliders (Cake))
        self.vmin_slider_cake = QSlider(Qt.Horizontal)
        self.vmin_slider_cake.setMinimum(0)
        self.vmin_slider_cake.setMaximum(110)
        self.vmin_slider_cake.setValue(0)
        self.vmin_slider_cake.valueChanged.connect(self.update_image_display_cake)        
        self.vmax_slider_cake = QSlider(Qt.Horizontal)
        self.vmax_slider_cake.setMinimum(0)
        self.vmax_slider_cake.setMaximum(110)        
        self.vmax_slider_cake.setValue(100)
        self.vmax_slider_cake.valueChanged.connect(self.update_image_display_cake)

        self.plottabs = QTabWidget()
        layout.addWidget(self.plottabs)
        self.XRDImageTab_GroupBox = QGroupBox()
        self.AzInt1D_Tab_GroupBox = QGroupBox()
        self.AzInt2D_Tab_GroupBox = QGroupBox()
        self.OtherPlots_Tab_GroupBox = QGroupBox()
        
        self.plottabs.addTab(self.XRDImageTab_GroupBox, 'XRD Image')
        self.plottabs.addTab(self.AzInt1D_Tab_GroupBox, 'XRD Azim. Integration 1D')
        self.plottabs.addTab(self.AzInt2D_Tab_GroupBox, 'XRD Azim. Integration 2D (Cake)')
        self.plottabs.addTab(self.OtherPlots_Tab_GroupBox, 'Other Plots')

        #### XRDImage_Tab
        XRDimagelayout = QVBoxLayout()
        self.XRDImageTab_GroupBox.setLayout(XRDimagelayout)
        self.fig = Figure(figsize=(8, 6))        
        self.plot_canvas = FigureCanvas(self.fig)
        self.fig.tight_layout()
        hlayouttoolbar = QHBoxLayout()
        XRDimagelayout.addLayout(hlayouttoolbar)
        self.toolbar = NavigationToolbar_sub(self.plot_canvas, self)
        hlayouttoolbar.addWidget(self.toolbar)
        hlayouttoolbar.addWidget(QLabel("ColorScale:"))
        sliderlayout = QVBoxLayout()
        hlayouttoolbar.addLayout(sliderlayout)
        sliderlayout.addWidget(self.vmin_slider)
        sliderlayout.addWidget(self.vmax_slider)
        hlayout_imageandButtons = QHBoxLayout()
        XRDimagelayout.addLayout(hlayout_imageandButtons)
        hlayout_imageandButtons.addWidget(self.plot_canvas)
        self.plot_canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        vlayoutbuttons = QVBoxLayout()
        hlayout_imageandButtons.addLayout(vlayoutbuttons)
        vlayoutbuttons.setAlignment(Qt.AlignTop)
        addbutton(vlayoutbuttons,'Azim. Integration',self.makeAzInt)
        addbutton(vlayoutbuttons,'Save Image (.tif)',self.save_image)
        # addbutton(vlayoutbuttons,'Azim. Integration (Cake)',self.makecake)

        self.ax = self.plot_canvas.figure.subplots()
        
        #### AzInt1D_Tab
        AzInt1Dlayout = QVBoxLayout()
        self.AzInt1D_Tab_GroupBox.setLayout(AzInt1Dlayout)
        self.fig_1d = Figure(figsize=(8, 6))        
        self.plot_canvas_1d = FigureCanvas(self.fig_1d)
        self.fig_1d.tight_layout()
        
        hlayouttoolbar_1d = QHBoxLayout()
        
        AzInt1Dlayout.addLayout(hlayouttoolbar_1d)
        self.toolbar_1d = NavigationToolbar_sub(self.plot_canvas_1d, self)
        hlayouttoolbar_1d.addWidget(self.toolbar_1d)
        hlayout_imageandButtons_1d = QHBoxLayout()
        
        AzInt1Dlayout.addLayout(hlayout_imageandButtons_1d)
        hlayout_imageandButtons_1d.addWidget(self.plot_canvas_1d)
        vlayoutbuttons_1d = QVBoxLayout()
        hlayout_imageandButtons_1d.addLayout(vlayoutbuttons_1d)
        vlayoutbuttons_1d.setAlignment(Qt.AlignTop)
        addbutton(vlayoutbuttons_1d,'Save Curve as ASCII',self.SaveAzInt1DData)
        self.ax_1d = self.plot_canvas_1d.figure.subplots()
        
        #### AzInt2D_Tab
        AzInt2Dlayout = QVBoxLayout()
        self.AzInt2D_Tab_GroupBox.setLayout(AzInt2Dlayout)
        self.fig_2d = Figure(figsize=(8, 6))        
        self.plot_canvas_2d = FigureCanvas(self.fig_2d)
        self.fig_2d.tight_layout()
        hlayouttoolbar_2d = QHBoxLayout()

        self.toolbar_2d = NavigationToolbar_sub(self.plot_canvas_2d, self)
        hlayouttoolbar_2d.addWidget(self.toolbar_2d)
        sliderlayout_2d = QVBoxLayout()
        hlayouttoolbar_2d.addLayout(sliderlayout_2d)
        sliderlayout_2d.addWidget(self.vmin_slider_cake)
        sliderlayout_2d.addWidget(self.vmax_slider_cake)
        AzInt2Dlayout.addLayout(hlayouttoolbar_2d)
        
        AzInt2Dlayout.addWidget(self.plot_canvas_2d)
        self.ax_2d = self.plot_canvas_2d.figure.subplots()

        #### OtherPlots Tab
        OtherPlotslayout = QVBoxLayout()
        self.OtherPlots_Tab_GroupBox.setLayout(OtherPlotslayout)
        self.fig_other = Figure(figsize=(8, 6))        
        self.plot_canvas_other = FigureCanvas(self.fig_other)
        self.fig_other.tight_layout()
        hlayouttoolbar_other = QHBoxLayout()
        OtherPlotslayout.addLayout(hlayouttoolbar_other)
        self.toolbar_other = NavigationToolbar_sub(self.plot_canvas_other, self)
        hlayouttoolbar_other.addWidget(self.toolbar_other)
        OtherPlotslayout.addWidget(self.plot_canvas_other)
        self.ax_other = self.plot_canvas_other.figure.subplots()
        
        
        # layout.addWidget(self.XRDImageTab_GroupBox)
        
        layout.addWidget(self.progress_bar)
        
        

        #### Placeholders
        self.SpectrumManager = None
        self.bunchfreq = None
        self.spec_EeV = None
        self.spec_Wpere = None
        self.cif_file = None
        self.structure_data = None
        self.im = None
        self.imgshown = None
        self.selectedmat = None  # This will hold the selected material
        self.API_KEY = None
        self.noiseimage = None
        self.maskfun = None
        self.cakeimgshown = None
        self.im_2d = None
        self.twotheta = None
        self.IntvsTh = None
        self.radial = None
        self.cakeintensity = None
        
        if self.ciffiletoload is not None:
            self.load_cif(self.ciffiletoload)
            
        #### Final init:
        self.makedetectorobj('makeCsIRayonix')
        
    # def move_groupbox_to_new_window(self,groupbox):
    #     new_window = NewWindow(groupbox)
    #     self.subwindows.append(new_window)
    #     groupbox.window = new_window
    #     new_window.setWindowTitle(groupbox.windowlabel)
    #     groupbox.undockgroupbutton.setText('')
    #     groupbox.undockgroupbutton.setIcon(self.style().standardIcon(QStyle.SP_TitleBarUnshadeButton))
    #     groupbox.undockgroupbutton.clicked.disconnect()
    #     groupbox.undockgroupbutton.clicked.connect(lambda: self.dockGroup(groupbox))
    #     new_window.show()
        
    # def dockGroup(self,groupbox):
    #     groupbox.dockingcontainer.addWidget(groupbox)
    #     groupbox.undockgroupbutton.setText('')
    #     groupbox.undockgroupbutton.setIcon(self.style().standardIcon(QStyle.SP_TitleBarShadeButton))
    #     groupbox.window.close()
    #     groupbox.undockgroupbutton.clicked.disconnect()
    #     groupbox.undockgroupbutton.clicked.connect(lambda: self.move_groupbox_to_new_window(groupbox))
    # def showparamwindow(self):
    #     self.paramdock.setVisible(True)
        
    
        
    
    def tray_icon_activated(self, reason):
        if reason == QSystemTrayIcon.DoubleClick:
            if self.isVisible():
                self.hide()
            else:
                self.show()
                self.activateWindow()

    def plotQE(self):
        if self.XRDDetector is None:
            self.reportmessage('Define a Detector First')
            return
        Edomain = np.linspace(5000,100000,200)
        
        try:
            thickness = float(self.ScintThickness_input.text())
        except ValueError:
            print('Bad Thickness Value')
            return
        
        QElist = []
        DQElist = []
        for energy in Edomain:
            QE = self.XRDDetector.QEfunc(
                    mat_thickness_inmm=thickness,
                    mat_density=self.XRDDetector.QE_materialdensity,
                    mat_formula=self.XRDDetector.QE_mat,energyineV=energy)
            DQE = self.XRDDetector.DQEfunc(
                    mat_thickness_inmm=thickness,
                    mat_density=self.XRDDetector.QE_materialdensity,
                    mat_formula=self.XRDDetector.QE_mat,energyineV=energy)   
            QElist.append(QE)
            DQElist.append(DQE)
        ax = self.ax_other
        ax.clear()
        
        ax.plot(Edomain,QElist)
        ax.plot(Edomain,DQElist)
        self.plot_canvas_other.draw_idle()
        self.plottabs.setCurrentWidget(self.OtherPlots_Tab_GroupBox)
        
        
    def findchemical(self):
        
        if self.matproject_checkbox.isChecked():
            self.findchemical_fromMatProject()
        
        if self.cod_checkbox.isChecked():
            self.findchemical_fromCOD()

    def findchemical_fromCOD(self):
        
        chemform = self.SpecifyChemForm_edit.text()
        
        # avog =6.0221408e23
        # unitcellmass = 0
        # for el,natom in zip(comps.Elements,comps.nAtoms):
        #     atomicweight = xraylib.AtomicWeight(el)/avog
        
        formattedformula = to_hill_notation(chemform)
        print(f'Formula sent to COD:{formattedformula}')

        
        if formattedformula=='':
            print('Invalid Formula')
            return
        
        try: 
            table_rows,table_header = search_cod_fortabledata(formattedformula)
            if table_rows==[]:
                print('No Data Retrieved')
                return
            dictlist=readtabledata(table_rows,table_header)

            if dictlist==[]:
                print('No COD records found.')
                return
            for el in dictlist:
                print('-------------')
                for key in el.keys():
                    print(f'{key}: {el[key]}')
        except: 
            print('Connection Error. Perhaps COD is down?')
            return    
                
                
                
        self.matwidget_COD = MaterialTableWidget_COD(self,dictlist)
        self.matwidget_COD.show()
        
        
    def findchemical_fromMatProject(self):
        chemform = self.SpecifyChemForm_edit.text()
        
        
        self.API_KEY = self.get_mp_api_key()
        if not self.API_KEY:
            print('No API Key')
            self.result_label.setText('No API Key')
            return
            
        with MPRester(self.API_KEY) as mpr: 
            
            try:
                materials = mpr.materials.summary.search(formula=chemform, fields=["material_id", "formula_pretty", "band_gap","density",'symmetry']) 
            except (ConnectionError,OSError):
                print('Connection Error. Perhaps Mat Project is down?')
                return
            except:
                self.reportmessage('Query Failed: Bad Formula?')
                return
            
            self.materials=materials
            if materials: 
                print(f"Materials found for {chemform}:") 
                for material in materials: 
                    # print(f"ID: {material['material_id']}, Formula: {material['pretty_formula']}, " 
                    #       f"Band Gap: {material['band_gap']} eV, Space Group: {material['spacegroup']['symbol']}") 
                    print(f"ID: {material.material_id}, Formula: {material.formula_pretty}, Density: {material.density:.03f}, {material.symmetry.crystal_system.value} Lattice, PtGroup: {material.symmetry.point_group}") 
                    # self.sym = material.symmetry
            else: 
                print(f"No materials found for {chemform}.")
                self.result_label.setText(f"No materials found for {chemform}.")
                return
            
        attributes = ['material_id','formula_pretty','density','symmetry.crystal_system.value','symmetry.point_group']
        attribute_labels = ['material_id','formula_pretty','density','System','Point Group']
        self.matwidget = MaterialTableWidget_MatProject(self, materials, attributes, attribute_labels,selectmethod_str = 'select_material_forCIF')
        self.matwidget.show()
        
        
        
    def showgeometry(self,showXRDcones = 1):

        #make filterlist
        filterlist = []

        for chemformula,thickness,Density,Angle_deg,zpos in zip(
                self.PostSampleFilters_GroupBox.activeCompounds,
                self.PostSampleFilters_GroupBox.activeThicknesses, 
                self.PostSampleFilters_GroupBox.activeDensities, 
                self.PostSampleFilters_GroupBox.activeAngles,
                self.PostSampleFilters_GroupBox.activezpos):
            filterlist.append({'thickness':thickness,'angle':Angle_deg,'mat':chemformula,'zpos':zpos})
            
         
        
        try: 
            Mx= 2*int(self.MHorizontal_input.text())+1
            Ny = 2*int(self.NVertical_input.text())+1
            pixelsizeinmm = float(self.pixelsize_input.text())*0.001
            ycenterpix = int(self.ycenterpix_input.text())
            xcenterpix = int(self.xcenterpix_input.text())
            detdistinmm = float(self.SampleToDetDistance_input.text())
            # energy = float(self.energy_input.text())
            thicknessinmm = float(self.thickness_input.text())
            Alpharadians =  float(self.sampleangle_input.text())*np.pi/180
            halfx = (Mx-1)/2 +1
            halfy = (Ny-1)/2 +1
        except ValueError: 
            self.reportmessage('Bad Values')
            return

        if showXRDcones==1:
            #make list of XRD cones
            
            corners_x=halfx*np.array([-1,-1,1,1])
            corners_y=halfy*np.array([-1, 1,1,-1])
            distcorners = pixelsizeinmm*np.sqrt(np.square(xcenterpix-halfx-corners_x)+np.square(ycenterpix-halfy-corners_y))
            theta_max = 180/np.pi*np.arctan(np.max(distcorners)/detdistinmm)
            
            if self.SpectrumManager is None:
                print('No Spectrum Specified')
                self.result_label.setText('No Spectrum Specified') 
                return
            if not self.cif_file or not self.structure_data:
                self.result_label.setText("Please select a valid CIF file.")
                return
                
            self.getspectrum()
            
            if self.spec_EeV is None:
                print('No Spectrum Specified')
                self.result_label.setText('No Spectrum Specified') 
                return
            self.ax.clear()
            
            #find spectral peaks
            x = self.spec_EeV
            y = self.spec_WpereV
            # Callback to update progress bar

            dx = x[1]-x[0]
            self.num_peaks, self.peak_locs = find_peaks_in_data(x, y, height=0.001, distance=int(4000/dx))
            lowestEpeak = np.min(self.peak_locs)

            def dont_update_progress(current, total):
                pass            
            # Compute structure factors
            hklvals, hkl_angles, structure_factors,multiplicities,elements,f0_values = compute_structure_factors(self.structure_data, theta_max, lowestEpeak, dont_update_progress)
            conelist = []  
            for twotheta,strfactor in zip(hkl_angles,structure_factors):
                if strfactor>1:
                    conelist.append(twotheta)
        else:
            conelist = None
                    
        geoobj = GeometryVisualizer.GeometryVisualizer(
            sample_thickness=thicknessinmm, 
            sample_angle=np.degrees(Alpharadians), 
            sample_to_detector_distance=detdistinmm, 
            pixel_size=pixelsizeinmm, 
            num_pixels=Mx,
            filterlist = filterlist, 
            conelist = conelist, 
            detectorxycenter = (pixelsizeinmm*(xcenterpix-halfx),
                                pixelsizeinmm*(ycenterpix-halfx)))
        self.subwindows.append(geoobj)
        geoobj.show()
        
    def makeAzInt(self):
        
        if self.noiseimage is None:
            return
        
        try:
            pixelsizeinmm = float(self.pixelsize_input.text())*0.001
            ycenterpix = int(self.ycenterpix_input.text())
            xcenterpix = int(self.xcenterpix_input.text())
            detdistinmm = float(self.SampleToDetDistance_input.text())
        except ValueError:
            self.result_label.setText("Please enter a valid energy value.")
            return
        
        energy = self.peak_locs[0] #
        #we're going to integrate into 2th anyway, just use the lowest energy peak. 
        
        self.twotheta,self.IntvsTh = AziIntegrate1D_single(img = self.noiseimage,d = detdistinmm,x0 = xcenterpix,y0 = ycenterpix,energy = energy,pixsize = pixelsizeinmm)
        self.radial, self.cakeintensity= AziIntegrate2D(img = self.noiseimage,d = detdistinmm,x0 = xcenterpix,y0 = ycenterpix,energy = energy,pixsize = pixelsizeinmm)
        
        self.plot1DIntegration(self.twotheta,self.IntvsTh)
        self.plotcake(self.radial,self.cakeintensity)
        
        self.plottabs.setCurrentWidget(self.AzInt1D_Tab_GroupBox)
        
        pass
    
    # def makecake(self):

    #     pass
        
    def CalcDebyeFactorFromMaterialsProject(self):
        self.API_KEY = self.get_mp_api_key()
        if not self.API_KEY:
            print('No API Key')
            self.result_label.setText('No API Key')
            return
            
        structure = self.structure_data
        with MPRester(self.API_KEY) as mpr: 
            # formula = structure.get_chemical_formula('reduce') # Query for all data on materials matching the formula 
            # 
            elements = self.structure_data.get_chemical_symbols()
            formula = chemformula(''.join(elements))                
            
            try:
                materials = mpr.materials.summary.search(formula=formula, has_props = [HasProps.elasticity],fields=["material_id", "formula_pretty", "band_gap","density",'symmetry']) 
            except (ConnectionError,OSError):
                print('Connection Error. Perhaps Mat Project is down?')
                return
            self.materials=materials
            if materials: 
                print(f"Materials found for {formula}:") 
                for material in materials: 
                    # print(f"ID: {material['material_id']}, Formula: {material['pretty_formula']}, " 
                    #       f"Band Gap: {material['band_gap']} eV, Space Group: {material['spacegroup']['symbol']}") 
                    print(f"ID: {material.material_id}, Formula: {material.formula_pretty}, Density: {material.density:.03f}, {material.symmetry.crystal_system.value} Lattice, PtGroup: {material.symmetry.point_group}") 
                    # self.sym = material.symmetry
            else: 
                print(f"No materials found for {formula}.")
                self.result_label.setText(f"No materials found for {formula}.")
                return
            
        attributes = ['material_id','formula_pretty','density','symmetry.crystal_system.value','symmetry.point_group']
        attribute_labels = ['material_id','formula_pretty','density','System','Point Group']
        self.matwidget = MaterialTableWidget_MatProject(self, materials, attributes, attribute_labels)
        self.matwidget.show()
        
    def save_poni_file(self, poni_file_name):
        """
        Save a .poni calibration file with simulation parameters.
        
        Parameters:
            poni_file_name: Full path to the .poni file to create
        """
        try:
            # Get parameters from GUI
            pixelsizeinmm = float(self.pixelsize_input.text()) * 0.001
            ycenterpix = int(self.ycenterpix_input.text())
            xcenterpix = int(self.xcenterpix_input.text())
            detdistinmm = float(self.SampleToDetDistance_input.text())
            
            # Use the lowest energy peak for wavelength calculation
            if self.peak_locs is None or len(self.peak_locs) == 0:
                print("Warning: No peak locations available, cannot create .poni file")
                return False
                
            energy = self.peak_locs[0]  # Use lowest energy peak
            
            # Setup azimuthal integrator parameters
            d = detdistinmm
            x0 = xcenterpix
            y0 = ycenterpix
            pixsize = pixelsizeinmm
            
            tiltval = 0
            tiltplanerotation = 0
            lam = hc / energy * 1e-10  # Convert to meters
            pyFAI_orientation = 3
            
            # Create detector
            detcustom = pyFAI.detectors.Detector(
                pixel1=pixsize/1e3, 
                pixel2=pixsize/1e3, 
                orientation=pyFAI_orientation
            )
            
            # Create azimuthal integrator and set geometry
            ai = AzimuthalIntegrator(detector=detcustom, wavelength=lam)
            #Handling some misorientation between Dioptas and this code
            orientation_val_for_Dioptas = 2
            max_shape_dim_X = 2 * int(self.MHorizontal_input.text()) + 1
            max_shape_dim_Y = 2 * int(self.NVertical_input.text()) + 1
            x0=x0+0.5
            y0=max_shape_dim_Y-y0-0.5

            ai.setFit2D(d, x0, y0, tilt=tiltval, tiltPlanRotation=tiltplanerotation)
            
            # Get pyFAI parameters
            poniparams = ai.getPyFAI()
            
            # Build detector name
            detector_name = "Simulated"
            if self.XRDDetector is not None and hasattr(self.XRDDetector, 'name'):
                detector_name = "Simulated_" + self.XRDDetector.name
            
            # Calculate pixel size in meters and max shape
            pixelsizeinm = float(self.pixelsize_input.text()) * 1e-6

            
            # Get orientation value (handle enum)
            #orientation_val = poniparams['orientation'].value if hasattr(poniparams['orientation'], 'value') else poniparams['orientation']
            

            
            # Write .poni file
            import datetime
            with open(poni_file_name, 'w') as f:
                f.write("# Nota: C-Order, 1 refers to the Y axis, 2 to the X axis\n")
                f.write(f"# Simulation done on {datetime.datetime.now().strftime('%a %b %d %H:%M:%S %Y')}\n")
                f.write("poni_version: 2.1\n")
                f.write(f"Detector: Detector\n")
                f.write(f'Detector_config: {{"pixel1": {pixelsizeinm}, "pixel2": {pixelsizeinm}, "max_shape": [{max_shape_dim_Y}, {max_shape_dim_X}], "orientation": {orientation_val_for_Dioptas}}}\n')
                f.write(f"Distance: {poniparams['dist']}\n")
                f.write(f"Poni1: {poniparams['poni1']}\n")
                f.write(f"Poni2: {poniparams['poni2']}\n")
                f.write(f"Rot1: {poniparams['rot1']}\n")
                f.write(f"Rot2: {poniparams['rot2']}\n")
                f.write(f"Rot3: {poniparams['rot3']}\n")
                f.write(f"Wavelength: {poniparams['wavelength']}\n")
            
            print(f'PONI file saved: {poni_file_name}')
            return True
            
        except ValueError as e:
            print(f"Error: Invalid parameter values for .poni file creation: {e}")
            return False
        except Exception as e:
            print(f"Error creating .poni file: {e}")
            return False
    
    def save_image(self):
        # self.recordsimparametersindict()
        if self.noiseimage is not None:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "TIFF Files (*.tif);;All Files (*)", options=options)
        
            if file_name:
                # Save TIFF image
                image = Image.fromarray(self.noiseimage)
                image.save(file_name, format='TIFF')
                
                # Save H5 file with simulation parameters
                h5_file_name = file_name.rsplit('.', 1)[0] + '.h5'
                try:
                    with h5py.File(h5_file_name, 'w') as h5f:
                        # Create groups for different parameter categories
                        detector_grp = h5f.create_group('detector_parameters')
                        beam_grp = h5f.create_group('beam_parameters')
                        target_grp = h5f.create_group('target_parameters')
                        filter_grp = h5f.create_group('post_sample_filters')
                        spectrum_grp = h5f.create_group('energy_spectrum')
                        
                        # Detector Parameters
                        try:
                            detector_grp.attrs['pixel_size_um'] = float(self.pixelsize_input.text())
                        except:
                            detector_grp.attrs['pixel_size_um'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['sample_to_det_distance_mm'] = float(self.SampleToDetDistance_input.text())
                        except:
                            detector_grp.attrs['sample_to_det_distance_mm'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['num_horizontal_pixels_half'] = int(self.MHorizontal_input.text())
                        except:
                            detector_grp.attrs['num_horizontal_pixels_half'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['num_vertical_pixels_half'] = int(self.NVertical_input.text())
                        except:
                            detector_grp.attrs['num_vertical_pixels_half'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['instrument_FWHM_mm'] = float(self.InstrFWHM_input.text())
                        except:
                            detector_grp.attrs['instrument_FWHM_mm'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['scintillator_thickness_mm'] = float(self.ScintThickness_input.text())
                        except:
                            detector_grp.attrs['scintillator_thickness_mm'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['read_noise_offset_counts'] = float(self.readnoiseoffset_input.text())
                        except:
                            detector_grp.attrs['read_noise_offset_counts'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['read_noise_std_counts'] = float(self.readnoiseSTD_input.text())
                        except:
                            detector_grp.attrs['read_noise_std_counts'] = 'N/A'
                        
                        try:
                            detector_grp.attrs['gain_ADU_per_keV'] = float(self.gainADUperkeV.text())
                        except:
                            detector_grp.attrs['gain_ADU_per_keV'] = 'N/A'
                        
                        detector_grp.attrs['simulate_noise'] = bool(self.SimNoiseCheckbox.isChecked())
                        
                        if self.XRDDetector is not None:
                            detector_grp.attrs['detector_name'] = self.XRDDetector.name if hasattr(self.XRDDetector, 'name') else 'Unknown'
                        
                        # Beam Parameters
                        try:
                            beam_grp.attrs['beam_center_horizontal_pixel'] = int(self.xcenterpix_input.text())
                        except:
                            beam_grp.attrs['beam_center_horizontal_pixel'] = 'N/A'
                        
                        try:
                            beam_grp.attrs['beam_center_vertical_pixel'] = int(self.ycenterpix_input.text())
                        except:
                            beam_grp.attrs['beam_center_vertical_pixel'] = 'N/A'
                        
                        beam_grp.attrs['polarized_beam'] = bool(self.polarizedcheck.isChecked())
                        beam_grp.attrs['enable_fluorescence'] = bool(self.enableFluorescence.isChecked())
                        beam_grp.attrs['only_fluorescence'] = bool(self.onlyfluorescence.isChecked())
                        beam_grp.attrs['enable_debye_waller'] = bool(self.DebyeWallerCheckbox.isChecked())
                        
                        if self.bunchfreq is not None:
                            beam_grp.attrs['bunch_frequency_Hz'] = self.bunchfreq
                        
                        # Target/Sample Parameters
                        target_grp.attrs['chemical_formula'] = self.chemform_out.text()
                        
                        try:
                            target_grp.attrs['thickness_mm'] = float(self.thickness_input.text())
                        except:
                            target_grp.attrs['thickness_mm'] = 'N/A'
                        
                        try:
                            target_grp.attrs['beam_angle_deg'] = float(self.sampleangle_input.text())
                        except:
                            target_grp.attrs['beam_angle_deg'] = 'N/A'
                        
                        target_grp.attrs['density_g_per_cc'] = self.density_out.text()
                        
                        if self.cif_file is not None:
                            target_grp.attrs['cif_file'] = self.cif_file
                        
                        if self.structure_data is not None:
                            try:
                                cell_params = self.structure_data.get_cell().cellpar()
                                target_grp.attrs['cell_a'] = cell_params[0]
                                target_grp.attrs['cell_b'] = cell_params[1]
                                target_grp.attrs['cell_c'] = cell_params[2]
                                target_grp.attrs['cell_alpha'] = cell_params[3]
                                target_grp.attrs['cell_beta'] = cell_params[4]
                                target_grp.attrs['cell_gamma'] = cell_params[5]
                                target_grp.attrs['cell_volume_angstrom3'] = self.structure_data.get_volume()
                            except:
                                pass
                        
                        # Debye-Waller Factors
                        if hasattr(self, 'atomDWstructure') and self.atomDWstructure is not None:
                            dw_grp = target_grp.create_group('debye_waller_factors')
                            dw_grp.attrs['description'] = "Debye-Waller B factors (in exp(-2B*s^2) form) for each element"
                            for atom_dw_data in self.atomDWstructure:
                                element = atom_dw_data['atom']
                                try:
                                    b_factor = float(atom_dw_data['editbox'].text())
                                    dw_grp.attrs[f'{element}_B_factor_angstrom2'] = b_factor
                                except:
                                    dw_grp.attrs[f'{element}_B_factor_angstrom2'] = 'N/A'
                        
                        # Post-Sample Filters
                        if hasattr(self.PostSampleFilters_GroupBox, 'activeCompounds'):
                            for i, (compound, thickness, density, angle, zpos) in enumerate(zip(
                                self.PostSampleFilters_GroupBox.activeCompounds,
                                self.PostSampleFilters_GroupBox.activeThicknesses,
                                self.PostSampleFilters_GroupBox.activeDensities,
                                self.PostSampleFilters_GroupBox.activeAngles,
                                self.PostSampleFilters_GroupBox.activezpos
                            )):
                                filter_subgrp = filter_grp.create_group(f'filter_{i}')
                                filter_subgrp.attrs['compound'] = compound
                                filter_subgrp.attrs['thickness_mm'] = thickness
                                filter_subgrp.attrs['density_g_per_cc'] = density
                                filter_subgrp.attrs['angle_deg'] = angle
                                filter_subgrp.attrs['z_position_mm'] = zpos
                        
                        
                        # Energy Spectrum
                        if self.spec_EeV is not None and self.spec_WpereV is not None:
                            spectrum_grp.create_dataset('energy_eV', data=self.spec_EeV, compression='gzip')
                            spectrum_grp.create_dataset('spectral_power_W_per_eV', data=self.spec_WpereV, compression='gzip')
                            spectrum_grp.attrs['description'] = 'Energy spectrum used in simulation'
                        
                        # CIF File Content
                        if self.cif_file is not None and os.path.exists(self.cif_file):
                            try:
                                with open(self.cif_file, 'r', encoding='utf-8') as cif:
                                    cif_content = cif.read()
                                    # Store CIF content as a string dataset
                                    cif_grp = h5f.create_group('cif_data')
                                    # Use special dtype for variable-length string
                                    dt = h5py.string_dtype(encoding='utf-8')
                                    cif_dataset = cif_grp.create_dataset('cif_file_content', data=cif_content, dtype=dt)
                                    cif_grp.attrs['cif_filename'] = os.path.basename(self.cif_file)
                                    cif_grp.attrs['cif_full_path'] = self.cif_file
                                    cif_grp.attrs['description'] = 'Complete CIF file content for crystal structure'
                                    print(f'CIF file content saved to H5: {os.path.basename(self.cif_file)}')
                            except Exception as e:
                                print(f'Warning: Could not save CIF file content: {e}')
                        
                        # Add timestamp
                        import datetime
                        h5f.attrs['creation_time'] = datetime.datetime.now().isoformat()
                        h5f.attrs['software'] = 'XRD Image Simulation'
                        
                        # Save SpectrumManager Configuration
                        if self.SpectrumManager is not None:
                            try:
                                spec_config_grp = h5f.create_group('SpectrumManager')
                                self.SpectrumManager.save_to_h5_group(spec_config_grp)
                            except Exception as e:
                                print(f"Warning: Could not save SpectrumManager configuration: {e}")
                        
                    print(f'H5 file saved: {h5_file_name}')
                except Exception as e:
                    print(f'Error saving H5 file: {e}')
                
                # Save .poni calibration file
                poni_file_name = file_name.rsplit('.', 1)[0] + '.poni'
                poni_success = self.save_poni_file(poni_file_name)
                
                # Update result label based on what was saved
                if poni_success:
                    self.result_label.setText(f'TIFF, H5, and PONI files saved successfully')
                else:
                    self.result_label.setText(f'TIFF and H5 files saved successfully (PONI creation skipped)')


    def LoadSimulation(self):
        """Load simulation parameters from an H5 file"""
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(
            self, 
            "Load Simulation", 
            "", 
            "H5 Files (*.h5);;All Files (*)", 
            options=options
        )
        
        if file_name:
            try:
                with h5py.File(file_name, 'r') as h5f:
                    if 'detector_parameters' in h5f:
                        det_grp = h5f['detector_parameters']
                        det_name = det_grp.attrs.get('detector_name', 'Unknown')
                        
                        # Find matching static method
                        det_method_name = None
                        for attr_name in dir(XRDDetector):
                            if attr_name.startswith('make'):
                                attr = getattr(XRDDetector, attr_name)
                                if callable(attr):
                                    try:
                                        # Create temp object to check name
                                        temp_det = attr()
                                        if hasattr(temp_det, 'name') and temp_det.name == det_name:
                                            det_method_name = attr_name
                                            break
                                    except:
                                        pass
                        
                        if det_method_name:
                            self.makedetectorobj(det_method_name)
                        else:
                            print(f"Warning: Could not find detector method for {det_name}")
                        
                        # Update detector parameters from H5
                        attrs = det_grp.attrs
                        if 'pixel_size_um' in attrs and attrs['pixel_size_um'] != 'N/A':
                            self.pixelsize_input.setText(str(attrs['pixel_size_um']))
                        if 'sample_to_det_distance_mm' in attrs and attrs['sample_to_det_distance_mm'] != 'N/A':
                            self.SampleToDetDistance_input.setText(str(attrs['sample_to_det_distance_mm']))
                        if 'num_horizontal_pixels_half' in attrs and attrs['num_horizontal_pixels_half'] != 'N/A':
                            self.MHorizontal_input.setText(str(attrs['num_horizontal_pixels_half']))
                        if 'num_vertical_pixels_half' in attrs and attrs['num_vertical_pixels_half'] != 'N/A':
                            self.NVertical_input.setText(str(attrs['num_vertical_pixels_half']))
                        if 'instrument_FWHM_mm' in attrs and attrs['instrument_FWHM_mm'] != 'N/A':
                            self.InstrFWHM_input.setText(str(attrs['instrument_FWHM_mm']))
                        if 'scintillator_thickness_mm' in attrs and attrs['scintillator_thickness_mm'] != 'N/A':
                            self.ScintThickness_input.setText(str(attrs['scintillator_thickness_mm']))
                        if 'read_noise_offset_counts' in attrs and attrs['read_noise_offset_counts'] != 'N/A':
                            self.readnoiseoffset_input.setText(str(attrs['read_noise_offset_counts']))
                        if 'read_noise_std_counts' in attrs and attrs['read_noise_std_counts'] != 'N/A':
                            self.readnoiseSTD_input.setText(str(attrs['read_noise_std_counts']))
                        if 'gain_ADU_per_keV' in attrs and attrs['gain_ADU_per_keV'] != 'N/A':
                            self.gainADUperkeV.setText(str(attrs['gain_ADU_per_keV']))
                        if 'simulate_noise' in attrs:
                            self.SimNoiseCheckbox.setChecked(bool(attrs['simulate_noise']))

                    if 'beam_parameters' in h5f:
                        beam_grp = h5f['beam_parameters']
                        attrs = beam_grp.attrs
                        if 'beam_center_horizontal_pixel' in attrs and attrs['beam_center_horizontal_pixel'] != 'N/A':
                            self.xcenterpix_input.setText(str(attrs['beam_center_horizontal_pixel']))
                        if 'beam_center_vertical_pixel' in attrs and attrs['beam_center_vertical_pixel'] != 'N/A':
                            self.ycenterpix_input.setText(str(attrs['beam_center_vertical_pixel']))
                        if 'polarized_beam' in attrs:
                            self.polarizedcheck.setChecked(bool(attrs['polarized_beam']))
                        if 'enable_fluorescence' in attrs:
                            self.enableFluorescence.setChecked(bool(attrs['enable_fluorescence']))
                        if 'only_fluorescence' in attrs:
                            self.onlyfluorescence.setChecked(bool(attrs['only_fluorescence']))
                        if 'enable_debye_waller' in attrs:
                            self.DebyeWallerCheckbox.setChecked(bool(attrs['enable_debye_waller']))
                        if 'bunch_frequency_Hz' in attrs:
                            self.bunchfreq = float(attrs['bunch_frequency_Hz'])

                    if 'target_parameters' in h5f:
                        target_grp = h5f['target_parameters']
                        attrs = target_grp.attrs
                        if 'thickness_mm' in attrs and attrs['thickness_mm'] != 'N/A':
                            self.thickness_input.setText(str(attrs['thickness_mm']))
                        if 'beam_angle_deg' in attrs and attrs['beam_angle_deg'] != 'N/A':
                            self.sampleangle_input.setText(str(attrs['beam_angle_deg']))
                        
                        # Load CIF
                        if 'cif_data' in h5f:
                            # Use existing method to load CIF from H5
                            self.load_cif_from_h5(file_name)
                        
                        # Load Debye-Waller Factors
                        if 'debye_waller_factors' in target_grp:
                            dw_grp = target_grp['debye_waller_factors']
                            if hasattr(self, 'atomDWstructure') and self.atomDWstructure:
                                for atom_data in self.atomDWstructure:
                                    element = atom_data['atom']
                                    key = f'{element}_B_factor_angstrom2'
                                    if key in dw_grp.attrs and dw_grp.attrs[key] != 'N/A':
                                        atom_data['editbox'].setText(str(dw_grp.attrs[key]))

                    if 'post_sample_filters' in h5f:
                        filter_grp = h5f['post_sample_filters']
                        # Clear existing filters
                        if hasattr(self.PostSampleFilters_GroupBox, 'table'):
                            for _ in range(self.PostSampleFilters_GroupBox.table.rowCount()):
                                self.PostSampleFilters_GroupBox.table.removeRow(0)
                            self.PostSampleFilters_GroupBox.activeCompounds = []
                            self.PostSampleFilters_GroupBox.activeThicknesses = []
                            self.PostSampleFilters_GroupBox.activeDensities = []
                            self.PostSampleFilters_GroupBox.activeAngles = []
                            self.PostSampleFilters_GroupBox.activezpos = []

                        # Add filters from H5
                        # Keys are filter_0, filter_1, etc. Sort them to maintain order.
                        filter_keys = sorted([k for k in filter_grp.keys() if k.startswith('filter_')], 
                                          key=lambda x: int(x.split('_')[1]))
                        
                        for key in filter_keys:
                            f_sub = filter_grp[key]
                            attrs = f_sub.attrs
                            self.PostSampleFilters_GroupBox.add_filter(
                                attrs['compound'],
                                float(attrs['thickness_mm']),
                                float(attrs['density_g_per_cc']),
                                float(attrs['angle_deg']),
                                float(attrs['z_position_mm']),
                                checked=1
                            )

                    if 'energy_spectrum' in h5f:
                        spec_grp = h5f['energy_spectrum']
                        if 'energy_eV' in spec_grp and 'spectral_power_W_per_eV' in spec_grp:
                            energy_dat = np.array(spec_grp['energy_eV'])
                            spectralpower_dat = np.array(spec_grp['spectral_power_W_per_eV'])
                            
                            if self.SpectrumManager is None:
                                self.SpectrumManager = DCSSpectrum.DCSSpectrum()
                                self.SpectrumManager.Emax_edit.setText(str(30000))
                                self.subwindows.append(self.SpectrumManager)
                            
                            self.SpectrumManager.energy_dat = energy_dat
                            self.SpectrumManager.spectralpower_dat = spectralpower_dat
                            
                            self.result_label.setText(f"Simulation loaded from {os.path.basename(file_name)}")
                            
                    # Load SpectrumManager Configuration
                    if 'SpectrumManager' in h5f:
                        try:
                            if self.SpectrumManager is None:
                                self.SpectrumManager = DCSSpectrum.DCSSpectrum()
                                self.subwindows.append(self.SpectrumManager)
                            
                            self.SpectrumManager.load_from_h5_group(h5f['SpectrumManager'])
                        except Exception as e:
                            print(f"Warning: Could not load SpectrumManager configuration: {e}")

            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load simulation: {e}")
                print(f"Error loading simulation: {e}")

    def suggest_sample_thickness(self):
        """ Suggest sample thickness based on attenuation length of spectral peaks."""
        if self.SpectrumManager is None:
            self.result_label.setText("No Spectrum Data!")
            QMessageBox.warning(self, "No Spectrum", "Please load or calculate a spectrum first.")
            return

        # Ensure spectrum is loaded/calculated
        self.getspectrum()
        
        if self.spec_EeV is None:
             self.result_label.setText("No Spectrum Data!")
             QMessageBox.warning(self, "No Spectrum", "Please load or calculate a spectrum first.")
             return
             
        # Check if structure/material is defined
        if not hasattr(self, 'structure_data') or self.structure_data is None:
             QMessageBox.warning(self, "No Material", "Please load a CIF file or define a material structure.")
             return

        x = self.spec_EeV
        y = self.spec_WpereV
        
        # Find peaks
        dx = x[1] - x[0] if len(x) > 1 else 1.0
        num_peaks, peak_locs = find_peaks_in_data(x, y, height=0.001, distance=int(4000/dx))
        
        if num_peaks == 0:
            QMessageBox.information(self, "No Peaks", "No spectral peaks found to analyze.")
            return
            
        density = calculate_density(self.structure_data)
        chem_formula = self.structure_data.get_chemical_formula()
        
        # Get Sample Angle
        try:
             sample_angle_deg = float(self.sampleangle_input.text())
             cos_theta = np.cos(np.deg2rad(sample_angle_deg))
        except:
             sample_angle_deg = 0
             cos_theta = 1.0
        
        # Calculate attenuation lengths
        atten_lengths = get_attenlengthinmm(chem_formula, density, peak_locs)
        
        # Prepare message
        msg = f"Material: {chem_formula}\nDensity: {density:.3f} g/cm³\n"
        msg += f"Beam Angle: {sample_angle_deg:.1f}°\n\n"
        msg += "Suggestions:\n"
        for energy, atten_len in zip(peak_locs, atten_lengths):
            recommended_thickness = atten_len * np.abs(cos_theta)
            msg += f"Peak: {energy:.1f} eV:\n"
            msg += f"  Rec. Path Length (Atten. Length): {atten_len:.4f} mm\n"
            msg += f"  Sample Thickness for {sample_angle_deg:.1f}°: {recommended_thickness:.4f} mm\n\n"
            
        QMessageBox.information(self, "Suggested Sample Thickness", msg)

    def suggest_xray_energy(self):
        """ Suggest X-ray energy matching the current sample thickness."""
        try:
            sample_thickness = float(self.thickness_input.text())
        except ValueError:
             QMessageBox.warning(self, "Invalid Thickness", "Please enter a valid numeric Sample Thickness.")
             return

        # Check if structure/material is defined
        if not hasattr(self, 'structure_data') or self.structure_data is None:
             QMessageBox.warning(self, "No Material", "Please load a CIF file or define a material structure.")
             return
             
        # Get Sample Angle
        try:
             sample_angle_deg = float(self.sampleangle_input.text())
             cos_theta = np.cos(np.deg2rad(sample_angle_deg))
        except:
             sample_angle_deg = 0
             cos_theta = 1.0
             
        # Effective path length through the sample
        target_path_length = sample_thickness / np.abs(cos_theta) if abs(cos_theta) > 1e-6 else sample_thickness * 1000 # Avoid divide by zero
             
        density = calculate_density(self.structure_data)
        chem_formula = self.structure_data.get_chemical_formula()
        
        # Generate energy range (100 eV to 100 keV)
        # Using more points for better resolution 
        energies = np.linspace(100, 100000, 2000)
        
        # Calculate attenuation lengths
        atten_lengths = get_attenlengthinmm(chem_formula, density, energies)
        
        # Find intersections where Attenuation Length == Target Path Length
        diff = atten_lengths - target_path_length
        # Detect indices where sign changes
        crossing_indices = np.where(np.diff(np.sign(diff)))[0]
        
        suggested_energies = []
        
        if len(crossing_indices) > 0:
            for idx in crossing_indices:
                # Linear interpolation for better precision
                e1, e2 = energies[idx], energies[idx+1]
                a1, a2 = atten_lengths[idx], atten_lengths[idx+1]
                
                if a2 != a1:
                    energy_crossing = e1 + (target_path_length - a1) * (e2 - e1) / (a2 - a1)
                    suggested_energies.append(energy_crossing)
                else:
                    suggested_energies.append(e1)
        else:
            # No exact match found, find closest
            idx = (np.abs(atten_lengths - target_path_length)).argmin()
            suggested_energies.append(energies[idx])

        msg = f"Material: {chem_formula}\n"
        msg += f"Sample Thickness: {sample_thickness} mm\n"
        msg += f"Beam Angle: {sample_angle_deg:.1f}°\n"
        msg += f"Target Path Length: {target_path_length:.4f} mm\n\n"
        msg += f"Suggested Energies (Atten. Length ~ Path Length):\n"
        
        for energy in suggested_energies:
            msg += f"  - {energy:.1f} eV ({energy/1000:.3f} keV)\n"
        
        # Add warning if no crossing was found
        if len(crossing_indices) == 0:
            msg += "\n(Note: No exact match found. Closest value shown.)"

        QMessageBox.information(self, "Suggested X-ray Energy", msg)
        
        # Plot Attenuation Length vs Energy
        ax = self.ax_other
        ax.clear()
        ax.plot(energies, atten_lengths, label='Attenuation Length')
        ax.axhline(y=target_path_length, color='r', linestyle='--', label=f'Target Path Length: {target_path_length:.3f} mm')
        
        # Mark all suggested energies
        for i, energy in enumerate(suggested_energies):
            label = f'{energy/1000:.1f} keV' if i == 0 else None # Label only first to avoid clutter? Or check count
            ax.axvline(x=energy, color='g', linestyle='--', label=label)
            # Add text annotation
            # ax.text(energy, target_path_length, f'{energy/1000:.1f} keV',rotation=90, verticalalignment='bottom')

        ax.set_title(f'Attenuation Length vs Energy: {chem_formula}')
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("Attenuation Length (mm)")
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="-")
        
        # Handle legend duplicate labels if multiple lines added
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        
        self.plot_canvas_other.draw()
        self.plottabs.setCurrentWidget(self.OtherPlots_Tab_GroupBox)

    def makedetectorobj(self,detectorstaticmethodstring):
        detectorstaticmethod = getattr(XRDDetector, detectorstaticmethodstring)
        self.XRDDetector = detectorstaticmethod()
        self.loaddetector(self.XRDDetector)

    def makeexptobj(self,ExptStaticMethod):
        expt_gen_method = getattr(ExperimentGeometry, ExptStaticMethod)
        self.Expt_Geometry = expt_gen_method()
        self.loadgeometry(self.Expt_Geometry)
        
    def loadgeometry(self,exptgeoobj):

        #clear filter list
        for el in reversed(range(self.PostSampleFilters_GroupBox.table.rowCount())):
            self.PostSampleFilters_GroupBox.table.removeRow(el)

        
        for num,filterobj in enumerate(exptgeoobj.filterlist,start=1):
            self.PostSampleFilters_GroupBox.add_filter(filterobj['mat'], filterobj['thickness'], filterobj['density'], filterobj['angle'],filterobj['zpos'],checked=1)
            
        self.sampleangle_input.setText(f'{exptgeoobj.sampleangle_deg:.3f}')
        self.SampleToDetDistance_input.setText(f'{exptgeoobj.detdistance:.3f}')
            # self.PostSampleFilters_GroupBox.table[]
            
        # exptgeo.filterlist=[
        #     {'mat':'LiF','density':2.64 ,'thickness':2.0,'angle':62,'zpos':0},
        #     {'mat':'C16H14O3','density':1.26 ,'thickness':8.0,'angle':-28,'zpos':80},
        #     ]
        # exptgeo.sampleangle_deg = 62
        # exptgeo.detdistance = 140
        
    def showimageinnewwindow(self,img):
        self.subwindows.append(ExtraFigureWindow(self))
        newfig = self.subwindows[-1]
        newfig.ax.imshow(img)
        # newfig.ax.set_aspect(self.aspectratio)
        newfig.show()
        
    def calcthetaphi_eachPixel(self,Mx,Ny,pixelsizeinmm, xcenterpix, ycenterpix, detdistinmm):
        imagecolnumbers_x = np.tile(np.arange(Mx).reshape(1, -1), (Ny, 1))
        imagerownumbers_y = np.tile(np.arange(Ny).reshape(-1, 1), (1, Mx))
    
        
        # distancefromcentermat = pixelsizeinmm*np.sqrt(np.square(imagecolnumbers_x-xcenterpix) + np.square(imagerownumbers_y-ycenterpix))
        # distancefromcentermat[ycenterpix,xcenterpix] = 1e-9
        
        '''Make 3x3 matrix of 2D theta,phi arrays
        with displacements from pixel center:
            [[(+h,+h),(+h,0),(+h,-h)],
             [(0,+h),(0,0),(0,-h)],
             [(-h,+h),(-h,0),(-h,-h)]]
        '''
        h = 1/2 #offset distance in pixels for simpsons rule integration
        offsets = np.array([[(+h,+h),(+h,0),(+h,-h)],
         [(0,+h),(0,0),(0,-h)],
         [(-h,+h),(-h,0),(-h,-h)]])
        offsetval = np.zeros((9,2))
        indnum = 0
        distancefromcentermat = pixelsizeinmm*np.sqrt(np.square(imagecolnumbers_x-xcenterpix) + np.square(imagerownumbers_y-ycenterpix))
        distancefromcentermat[distancefromcentermat<1e-9] = 1e-9
        rmat = distancefromcentermat
        Twothetas_rad = np.arctan(distancefromcentermat/detdistinmm); #array of 2*theta values at each pixel center
        Phi_rad = np.arccos(pixelsizeinmm*(imagecolnumbers_x-xcenterpix)/distancefromcentermat) # array of phi values at each pixel center
        return Twothetas_rad, Phi_rad, rmat
    
    def calcthetaphi_eachPixel_1x1(self,Mx,Ny,pixelsizeinmm, xcenterpix, ycenterpix, detdistinmm):
        """
        Calculate geometry for the center of each pixel only (optimization over 9x1 version).
        
        Returns:
            Twothetas_rad
            Phi_rad
            rmat
        """
        imagecolnumbers_x = np.tile(np.arange(Mx).reshape(1, -1), (Ny, 1))
        imagerownumbers_y = np.tile(np.arange(Ny).reshape(-1, 1), (1, Mx))
        
        distancefromcentermat = pixelsizeinmm*np.sqrt(np.square(imagecolnumbers_x-xcenterpix) + np.square(imagerownumbers_y-ycenterpix))
        distancefromcentermat[distancefromcentermat<1e-9] = 1e-9
        
        # In the 9x1 version, rmat is set only if dx==0 and dy==0, which is this case.
        rmat = distancefromcentermat

        Twothetas_rad = np.arctan(distancefromcentermat/detdistinmm)
        Phi_rad = np.arccos(pixelsizeinmm*(imagecolnumbers_x-xcenterpix)/distancefromcentermat) # array of phi values at each pixel center

        return Twothetas_rad, Phi_rad, rmat


    def calcthetaphi_eachPixel_9x1(self,Mx,Ny,pixelsizeinmm, xcenterpix, ycenterpix, detdistinmm):
        imagecolnumbers_x = np.tile(np.arange(Mx).reshape(1, -1), (Ny, 1))
        imagerownumbers_y = np.tile(np.arange(Ny).reshape(-1, 1), (1, Mx))
    
        
        # distancefromcentermat = pixelsizeinmm*np.sqrt(np.square(imagecolnumbers_x-xcenterpix) + np.square(imagerownumbers_y-ycenterpix))
        # distancefromcentermat[ycenterpix,xcenterpix] = 1e-9
        
        '''Make 3x3 matrix of 2D theta,phi arrays
        with displacements from pixel center:
            [[(+h,+h),(+h,0),(+h,-h)],
             [(0,+h),(0,0),(0,-h)],
             [(-h,+h),(-h,0),(-h,-h)]]
        '''
        h = 1/2 #offset distance in pixels for simpsons rule integration
        offsets = np.array([[(+h,+h),(+h,0),(+h,-h)],
         [(0,+h),(0,0),(0,-h)],
         [(-h,+h),(-h,0),(-h,-h)]])
        
        Twothetas_rad_9x1 = np.zeros((9,Ny,Mx))
        Phi_rad_9x1 = np.zeros((9,Ny,Mx))
        weights_simpsonOneThird_9x1 = [1,4,1,4,16,4,1,4,1] #weights for simpson's 1/3 rule
        offsetval = np.zeros((9,2))
        indnum = 0
        for ij, row in enumerate(offsets):
            for jk, el in enumerate(row):
                dy = el[0]
                dx = el[1]
                # offsetval[indnum] = [dy,dx]
                distancefromcentermat = pixelsizeinmm*np.sqrt(np.square(imagecolnumbers_x+dx-xcenterpix) + np.square(imagerownumbers_y+dy-ycenterpix))
                distancefromcentermat[distancefromcentermat<1e-9] = 1e-9
                if dx==0 and dy==0:
                    rmat = distancefromcentermat
                Twothetas_rad_9x1[indnum] = np.arctan(distancefromcentermat/detdistinmm); #array of 2*theta values at each pixel center
                Phi_rad_9x1[indnum] = np.arccos(pixelsizeinmm*(imagecolnumbers_x+dx-xcenterpix)/distancefromcentermat) # array of phi values at each pixel center
                
                indnum+=1
        # dPhi = Phi_rad_9x1[]
                    
        return Twothetas_rad_9x1, Phi_rad_9x1, weights_simpsonOneThird_9x1, rmat    
    
    
    def showspectrum(self):
        
        if self.SpectrumManager is None:
            self.result_label.setText("No Spectrum Data!")
            return
        
        if self.SpectrumManager.spectralpower_dat is None:
            self.result_label.setText("No Spectrum Data!")
            return
        
        self.getspectrum()
        
        x = self.spec_EeV
        y = self.spec_WpereV
        
        dx = x[1]-x[0]
        #Find peaks
        ax = self.ax_other
        num_peaks, peak_locs = find_peaks_in_data(x, y, height=0.001, distance=int(4000/dx))
        print(f"Number of peaks: {num_peaks}")
        print(f"Peak locations: {peak_locs}")
        
            
        
        ax.clear()

        ax.plot(x, y, '-')
        ax.scatter(peak_locs, [y[np.where(x == loc)[0][0]] for loc in peak_locs],
                color='red', label='Peaks', zorder=5)
        
        ax.set_title('Input Power Spectrum')
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("dP/dE (W/eV)")
        ax.grid(True)
        self.plot_canvas_other.draw()
        self.plottabs.setCurrentWidget(self.OtherPlots_Tab_GroupBox)

        self.result_label.setText("Plot updated successfully!")

    def select_file(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(
            self, 
            "Select CIF or H5 File", 
            "", 
            "CIF and H5 Files (*.cif *.h5);;CIF Files (*.cif);;H5 Files (*.h5);;All Files (*)"
        )
        
        if file_path:
            # Check if it's an H5 file
            if file_path.lower().endswith('.h5'):
                self.load_cif_from_h5(file_path)
            else:
                # Regular CIF file
                self.load_cif(file_path)
    
    def load_cif_from_h5(self, h5_file_path):
        """Extract CIF data from H5 file and load it"""
        try:
            with h5py.File(h5_file_path, 'r') as h5f:
                # Check if cif_data group exists
                if 'cif_data' not in h5f:
                    QMessageBox.warning(
                        self,
                        "No CIF Data",
                        "This H5 file does not contain CIF data.\n"
                        "Please select a CIF file or an H5 file with cif_data."
                    )
                    return
                
                # Extract CIF content
                cif_content = h5f['cif_data']['cif_file_content'][()]
                original_filename = h5f['cif_data'].attrs.get('cif_filename', 'extracted.cif')
                
                # Convert bytes to string if necessary
                if isinstance(cif_content, bytes):
                    cif_content = cif_content.decode('utf-8')
                
                # Create temporary CIF file
                import tempfile
                temp_dir = tempfile.gettempdir()
                temp_cif_path = os.path.join(temp_dir, f'from_h5_{original_filename}')
                
                # Write CIF content to temporary file
                with open(temp_cif_path, 'w', encoding='utf-8') as cif_file:
                    cif_file.write(cif_content)
                
                print(f"Extracted CIF from H5: {original_filename}")
                print(f"Temporary CIF file: {temp_cif_path}")
                
                # Load the extracted CIF file
                self.load_cif(temp_cif_path)
                
                # Update status
                self.result_label.setText(f"Loaded CIF from H5: {original_filename}")
                
        except Exception as e:
            QMessageBox.critical(
                self,
                "Error Loading H5 File",
                f"Failed to extract CIF from H5 file:\n{e}"
            )
            print(f"Error loading CIF from H5: {e}")

        
        
                
    def load_cif(self,file_path):
        if file_path:
            self.cif_file = file_path
            print(self.cif_file)
            # self.file_label.setText(f"CIF File: {file_path}")
            self.setWindowTitle(f"XRD Image Simulation: {file_path}")
            try:
                
                self.structure_data = read(file_path)
                
            except Exception as e:
                self.result_label.setText(f"Error reading CIF file: {e}")
                return
            elements = self.structure_data.get_chemical_symbols()
            chemform = chemformula(''.join(elements))                
            
            self.chemform_out.setText(chemform)
            self.density_out.setText(f'{calculate_density(self.structure_data):.3}')
            self.makeDebyeWallerGroupBox()
            DWfound = 0
            elements = self.structure_data.get_chemical_symbols()
            unique_elements = set(elements)
            if len(unique_elements)==1:
                DW = interpolate_debye_waller(list(unique_elements)[0], 300)
                if DW is not None: 
                    DWfound = 1
                    for DWatom in self.atomDWstructure:
                        DWatom['editbox'].setText(f'{DW:0.3f}')
            return DWfound
                    
                    
            

            
                
        
    def makeDebyeWallerGroupBox(self):
        if self.DWBox is not None:
            self.DWBox.setParent(None)
            delattr(self,'DWBox')
        self.DWBox = QGroupBox('Debye Parameter ( \'B\' of exp(-2B * s^2)) ')
        self.DWlayout.addWidget(self.DWBox)
        
        overallHbox = QHBoxLayout()
        self.DWBox.setLayout(overallHbox)
        
        uniqueel = set(self.structure_data.get_chemical_symbols())
        
        atomDW = []
        for atom in uniqueel:
            colatom = QVBoxLayout()
            overallHbox.addLayout(colatom)
            colatom.addWidget(QLabel(atom))
            DWedit = QLineEdit()
            DWedit.setPlaceholderText('DWFac')
            DWedit.setToolTip('Debye-Waller Factor')
            DWedit.setText(f'{0.5}')
            colatom.addWidget(DWedit)
            atomDWel = {'atom':atom}
            atomDWel['editbox'] = DWedit
            atomDW.append(atomDWel)
        self.atomDWstructure = atomDW
            
            
        
    def calculate_single(self):
        try:
            energy = float(self.energy_input.text())  # Energy in eV
        except ValueError:
            self.result_label.setText("Please enter a valid energy value.")
            return
        elements = set(atom.symbol for atom in self.structure_data)  # Extract unique elements
        f0_values = precompute_f0(xraylib, elements, x_min=0.0, x_max=8.0, num_points=100)
        h = 2
        k = 2
        l = 0
        structurefac=calculate_structure_factor(self.structure_data,h, k, l, energy, f0_values)
        
        print(f'F[{h},{k},{l}] = {structurefac}')
        
    def calctotalfluorescence(self):
        """
        Calculate the total number of fluorescence photons and their energies.
    
        Returns:
            fluorescence_dict: dict
                A dictionary with element symbols as keys and fluorescence photon counts and energies as values.
        """
        if self.SpectrumManager is None:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified')
            return            
        self.getspectrum()
        if self.spec_EeV is None:
            self.result_label.setText('No Spectrum Specified')
            return
        self.ax.clear()
        #find spectral peaks
        x = self.spec_EeV
        y = self.spec_WpereV
        dx = x[1]-x[0]
        self.num_peaks, self.peak_locs = find_peaks_in_data(x, y, height=0.001, distance=int(4000/dx))
        
        if self.num_peaks>1:
            EnergyBetweenPeaks = abs(self.peak_locs[1]-self.peak_locs[0])
        else: 
            EnergyBetweenPeaks = (np.max(x)-np.min(x))
            
        self.result_label.setText('Calculating Fluorescence...')    
        # Extract atomic information
        elements = self.structure_data.get_chemical_symbols()
        totaldensity = calculate_density(self.structure_data)
        unique_elements = set(elements)
        compstr = xraylib.CompoundParser(''.join(elements))
        self.attenfactorinmm_list = get_attenlengthinmm(self.structure_data.get_chemical_formula(),totaldensity, self.peak_locs)
        allemissionlines = np.flip(np.arange(-190,0))
        try:
            thicknessinmm = float(self.thickness_input.text())
            Alpharadians =  float(self.sampleangle_input.text())*np.pi/180
        except ValueError: 
            print('Bad Values for thickness/Alpha')
            self.result_label.setText('Bad Values for thickness/Alpha')

        pathinmaterial=thicknessinmm/np.cos(Alpharadians)
        emissionenergies = [] #emission energies in eV
        Nemitted = []
        
        for energy,attenlength0 in zip(self.peak_locs,self.attenfactorinmm_list):
            energy_keV = energy/1000
            
            
            minE = energy - EnergyBetweenPeaks/2.1
            maxE = energy + EnergyBetweenPeaks/2.1
            harmonic_E = x[(x>minE)&(x<maxE)] 
            harmonic_dPdE = y[(x>minE)&(x<maxE)] 
            harmonic_dNdE = harmonic_dPdE/harmonic_E
            
            N_E = np.trapezoid(harmonic_dNdE,harmonic_E)/self.bunchfreq * JouleIneV

            
            # Loop through unique elements and calculate fluorescence
            for Z,massfrac in zip(compstr['Elements'],compstr['massFractions']):
                
                density_thisatom=totaldensity*massfrac

                for line in allemissionlines:
                    good=0
                    try: 
                        fluor_production_CS_Kissel_Cascade=xraylib.CS_FluorLine_Kissel_Cascade(Z,line,energy_keV)
                        lineenergy = xraylib.LineEnergy(Z, line)
                        good = 1
                    except ValueError: 
                        pass
                    if good:
                        
                        producedperlength_inmm = (fluor_production_CS_Kissel_Cascade*density_thisatom)/10 #div by 10 due to 1/cm in units
                        totalNproduced = N_E*producedperlength_inmm*pathinmaterial*np.exp(-pathinmaterial/attenlength0)
                        emissionenergies.append(lineenergy*1000)#emission energies in eV
                        Nemitted.append(totalNproduced)
                        # print(f'Line#{-1*line}, N_emitted: {totalNproduced:.2e}')
                    
        # edat,Ndat = sum_y_for_same_x(emissionenergies,Nemitted)
        dE = 50
        edat,Ndat =  sum_y_within_maxsep(emissionenergies,Nemitted,dE, dE)
        edat, Ndat = (list(t) for t in zip(*sorted(zip(edat, Ndat))))
        totalinput_N = np.trapezoid(y/x,x)/self.bunchfreq * JouleIneV
        totalinput_E = np.trapezoid(y,x)/self.bunchfreq*1e6 # in µJ
        totalfluor_N = np.sum(np.array(Ndat))
        totalfluor_E = np.sum(np.array(edat)*np.array(Ndat))/JouleIneV*1e6 # in µJ
        fluor_energyfraction = totalfluor_E/totalinput_E
        print('Fluorescence:')
        print(f'#Phot_in: {totalinput_N:.3e}, E_in: {totalinput_E:.3e} µJ, #Fluor_emit: {totalfluor_N:.3e}, fluor_E: {totalfluor_E:.3e} µJ, E_fluorfraction: {fluor_energyfraction:.2f}')
        
        def fluorplot(edata,Ndata):
            x = np.linspace(0.8*np.min(edata),1.1*np.max(edata),200)
            y = np.zeros_like(x)
            # for e,N in zip(edata,Ndata):
            #     x.append(e)
            #     y.append(N)
            
            x=np.append(x,edata)
            y=np.append(y,Ndata)
            x, y = (list(t) for t in zip(*sorted(zip(x, y))))
            return x,y
        
        x,y = fluorplot(edat,Ndat)
        self.ax_other.clear()
        self.ax_other.plot(x,y, '-')
        # self.ax.set_yscale('log')
        self.plot_canvas_other.draw()
        self.plottabs.setCurrentWidget(self.OtherPlots_Tab_GroupBox)
        self.result_label.setText('Done') 
            

        '''
        Line Macros
        Transition	IUPAC macro
        K ← L1	KL1_LINE
        K ← L2	KL2_LINE
        K ← L3	KL3_LINE
        K ← M1	KM1_LINE
        K ← M2	KM2_LINE
        K ← M3	KM3_LINE
        K ← M4	KM4_LINE
        K ← M5	KM5_LINE
        K ← N1	KN1_LINE
        K ← N2	KN2_LINE
        K ← N3	KN3_LINE
        K ← N4	KN4_LINE
        K ← N5	KN5_LINE
        K ← N6	KN6_LINE
        K ← N7	KN7_LINE
        K ← O	KO_LINE
        K ← O1	KO1_LINE
        K ← O2	KO2_LINE
        K ← O3	KO3_LINE
        K ← O4	KO4_LINE
        K ← O5	KO5_LINE
        K ← O6	KO6_LINE
        K ← O7	KO7_LINE
        K ← P	KP_LINE
        K ← P1	KP1_LINE
        K ← P2	KP2_LINE
        K ← P3	KP3_LINE
        K ← P4	KP4_LINE
        K ← P5	KP5_LINE
        L1 ← L2	L1L2_LINE
        L1 ← L3	L1L3_LINE
        L1 ← M1	L1M1_LINE
        L1 ← M2	L1M2_LINE
        L1 ← M3	L1M3_LINE
        L1 ← M4	L1M4_LINE
        L1 ← M5	L1M5_LINE
        L1 ← N1	L1N1_LINE
        L1 ← N2	L1N2_LINE
        L1 ← N3	L1N3_LINE
        L1 ← N4	L1N4_LINE
        L1 ← N5	L1N5_LINE
        L1 ← N6	L1N6_LINE
        L1 ← N67	L1N67_LINE
        L1 ← N7	L1N7_LINE
        L1 ← O1	L1O1_LINE
        L1 ← O2	L1O2_LINE
        L1 ← O3	L1O3_LINE
        L1 ← O4	L1O4_LINE
        L1 ← O45	L1O45_LINE
        L1 ← O5	L1O5_LINE
        L1 ← O6	L1O6_LINE
        L1 ← O7	L1O7_LINE
        L1 ← P1	L1P1_LINE
        L1 ← P2	L1P2_LINE
        L1 ← P23	L1P23_LINE
        L1 ← P3	L1P3_LINE
        L1 ← P4	L1P4_LINE
        L1 ← P5	L1P5_LINE
        L2 ← L3	L2L3_LINE
        L2 ← M1	L2M1_LINE
        L2 ← M2	L2M2_LINE
        L2 ← M3	L2M3_LINE
        L2 ← M4	L2M4_LINE
        L2 ← M5	L2M5_LINE
        L2 ← N1	L2N1_LINE
        L2 ← N2	L2N2_LINE
        L2 ← N3	L2N3_LINE
        L2 ← N4	L2N4_LINE
        L2 ← N5	L2N5_LINE
        L2 ← N6	L2N6_LINE
        L2 ← N7	L2N7_LINE
        L2 ← O1	L2O1_LINE
        L2 ← O2	L2O2_LINE
        L2 ← O3	L2O3_LINE
        L2 ← O4	L2O4_LINE
        L2 ← O5	L2O5_LINE
        L2 ← O6	L2O6_LINE
        L2 ← O7	L2O7_LINE
        L2 ← P1	L2P1_LINE
        L2 ← P2	L2P2_LINE
        L2 ← P23	L2P23_LINE
        L2 ← P3	L2P3_LINE
        L2 ← P4	L2P4_LINE
        L2 ← P5	L2P5_LINE
        L2 ← Q1	L2Q1_LINE
        L3 ← M1	L3M1_LINE
        L3 ← M2	L3M2_LINE
        L3 ← M3	L3M3_LINE
        L3 ← M4	L3M4_LINE
        L3 ← M5	L3M5_LINE
        L3 ← N1	L3N1_LINE
        L3 ← N2	L3N2_LINE
        L3 ← N3	L3N3_LINE
        L3 ← N4	L3N4_LINE
        L3 ← N5	L3N5_LINE
        L3 ← N6	L3N6_LINE
        L3 ← N7	L3N7_LINE
        L3 ← O1	L3O1_LINE
        L3 ← O2	L3O2_LINE
        L3 ← O3	L3O3_LINE
        L3 ← O4	L3O4_LINE
        L3 ← O45	L3O45_LINE
        L3 ← O5	L3O5_LINE
        L3 ← O6	L3O6_LINE
        L3 ← O7	L3O7_LINE
        L3 ← P1	L3P1_LINE
        L3 ← P2	L3P2_LINE
        L3 ← P23	L3P23_LINE
        L3 ← P3	L3P3_LINE
        L3 ← P4	L3P4_LINE
        L3 ← P45	L3P45_LINE
        L3 ← P5	L3P5_LINE
        L3 ← Q1	L3Q1_LINE
        M1 ← M2	M1M2_LINE
        M1 ← M3	M1M3_LINE
        M1 ← M4	M1M4_LINE
        M1 ← M5	M1M5_LINE
        M1 ← N1	M1N1_LINE
        M1 ← N2	M1N2_LINE
        M1 ← N3	M1N3_LINE
        M1 ← N4	M1N4_LINE
        M1 ← N5	M1N5_LINE
        M1 ← N6	M1N6_LINE
        M1 ← N7	M1N7_LINE
        M1 ← O1	M1O1_LINE
        M1 ← O2	M1O2_LINE
        M1 ← O3	M1O3_LINE
        M1 ← O4	M1O4_LINE
        M1 ← O5	M1O5_LINE
        M1 ← O6	M1O6_LINE
        M1 ← O7	M1O7_LINE
        M1 ← P1	M1P1_LINE
        M1 ← P2	M1P2_LINE
        M1 ← P3	M1P3_LINE
        M1 ← P4	M1P4_LINE
        M1 ← P5	M1P5_LINE
        M2 ← M3	M2M3_LINE
        M2 ← M4	M2M4_LINE
        M2 ← M5	M2M5_LINE
        M2 ← N1	M2N1_LINE
        M2 ← N2	M2N2_LINE
        M2 ← N3	M2N3_LINE
        M2 ← N4	M2N4_LINE
        M2 ← N5	M2N5_LINE
        M2 ← N6	M2N6_LINE
        M2 ← N7	M2N7_LINE
        M2 ← O1	M2O1_LINE
        M2 ← O2	M2O2_LINE
        M2 ← O3	M2O3_LINE
        M2 ← O4	M2O4_LINE
        M2 ← O5	M2O5_LINE
        M2 ← O6	M2O6_LINE
        M2 ← O7	M2O7_LINE
        M2 ← P1	M2P1_LINE
        M2 ← P2	M2P2_LINE
        M2 ← P3	M2P3_LINE
        M2 ← P4	M2P4_LINE
        M2 ← P5	M2P5_LINE
        M3 ← M4	M3M4_LINE
        M3 ← M5	M3M5_LINE
        M3 ← N1	M3N1_LINE
        M3 ← N2	M3N2_LINE
        M3 ← N3	M3N3_LINE
        M3 ← N4	M3N4_LINE
        M3 ← N5	M3N5_LINE
        M3 ← N6	M3N6_LINE
        M3 ← N7	M3N7_LINE
        M3 ← O1	M3O1_LINE
        M3 ← O2	M3O2_LINE
        M3 ← O3	M3O3_LINE
        M3 ← O4	M3O4_LINE
        M3 ← O5	M3O5_LINE
        M3 ← O6	M3O6_LINE
        M3 ← O7	M3O7_LINE
        M3 ← P1	M3P1_LINE
        M3 ← P2	M3P2_LINE
        M3 ← P3	M3P3_LINE
        M3 ← P4	M3P4_LINE
        M3 ← P5	M3P5_LINE
        M3 ← Q1	M3Q1_LINE
        M4 ← M5	M4M5_LINE
        M4 ← N1	M4N1_LINE
        M4 ← N2	M4N2_LINE
        M4 ← N3	M4N3_LINE
        M4 ← N4	M4N4_LINE
        M4 ← N5	M4N5_LINE
        M4 ← N6	M4N6_LINE
        M4 ← N7	M4N7_LINE
        M4 ← O1	M4O1_LINE
        M4 ← O2	M4O2_LINE
'''

    def calcfluorescenceimage(self):
        """
        Calculate the total number of fluorescence photons and their energies.
    
        Returns:
            fluorescence_dict: dict
                A dictionary with element symbols as keys and fluorescence photon counts and energies as values.
        """
             
        self.getspectrum()
        self.ax.clear()
        #find spectral peaks
        x = self.spec_EeV
        y = self.spec_WpereV
                
        if self.num_peaks>1:
            EnergyBetweenPeaks = abs(self.peak_locs[1]-self.peak_locs[0])
        else: 
            EnergyBetweenPeaks = (np.max(x)-np.min(x))
            
        self.result_label.setText('Calculating Fluorescence...')     
        # Extract atomic information
        elements = self.structure_data.get_chemical_symbols()
        totaldensity = calculate_density(self.structure_data)
        compstr = xraylib.CompoundParser(''.join(elements))
        
        #instead of using all the xraylib macros for the lines, just make a list from -190 to -1
        allemissionlines = np.flip(np.arange(-190,0))
        
        try:
            thicknessinmm = float(self.thickness_input.text())
            Alpharadians =  float(self.sampleangle_input.text())*np.pi/180
        except ValueError: 
            print('Bad Values for thickness/Alpha')

        pathinmaterial=thicknessinmm/np.cos(Alpharadians)
        emissionenergies_tot = [] #emission energies in eV
        Nemitted_total = []
        
        #Calc Geometric Factors
        solidanglemat = self.dPhi * self.dtwotheta * np.sin(self.TwoTheta_rad_9x1[4])
        gamma,beta,boolmask_infattenuation = AttenuationAngleFactors(Alpharadians, self.TwoTheta_rad_9x1[4], self.Phi_Rad_9x1[4])
        # self.showimageinnewwindow(solidanglemat)
        
        self.FluorescentImage = []
        self.FluorLineEnergybins = []
        
        
        for energy,attenlength0 in zip(self.peak_locs,self.attenfactorinmm_list):
            NemittedperAngstrompersolidangle_thisexcitationenergy = []
            emissionenergies = []
            energy_keV = energy/1000
            
            minE = energy - EnergyBetweenPeaks/2.1
            maxE = energy + EnergyBetweenPeaks/2.1
            harmonic_E = x[(x>minE)&(x<maxE)] 
            harmonic_dPdE = y[(x>minE)&(x<maxE)] 
            harmonic_dNdE = harmonic_dPdE/harmonic_E
            
            N_E = np.trapezoid(harmonic_dNdE,harmonic_E)/self.bunchfreq * JouleIneV

            
            # Loop through unique elements and calculate fluorescence for each line
            for Z,massfrac in zip(compstr['Elements'],compstr['massFractions']):
                density_thisatom=totaldensity*massfrac

                for line in allemissionlines:
                    good=0
                    try: 
                        fluor_production_CS_Kissel_Cascade=xraylib.CS_FluorLine_Kissel_Cascade(Z,line,energy_keV)
                        lineenergy = xraylib.LineEnergy(Z, line)
                        good = 1
                    except ValueError: 
                        pass
                    if good:
                        
                        producedperAngstrom_persolidangle =(1/4/np.pi)*N_E*(fluor_production_CS_Kissel_Cascade*density_thisatom)/1e8 #div by 1e8 to go from 1/cm to 1/Angstrom
                        totalNproduced = 4*np.pi*producedperAngstrom_persolidangle * pathinmaterial
                        Nemitted_total.append(totalNproduced)
                        emissionenergies_tot.append(lineenergy*1000)

                        emissionenergies.append(lineenergy*1000)#emission energies in eV                        
                        NemittedperAngstrompersolidangle_thisexcitationenergy.append(producedperAngstrom_persolidangle)
                        # print(f'Line#{-1*line}, N_emitted: {totalNproduced:.2e}')
            
            #bin the fluorescence from neighboring lines:
            energybins ,NemitperAngstrompersolidangle =  sum_y_within_maxsep(emissionenergies,NemittedperAngstrompersolidangle_thisexcitationenergy,1000, 1000)
            
            #Calc Attenuation lengths for each fluorescence bin
            attenfactorinmm_list_fluor = get_attenlengthinmm(self.structure_data.get_chemical_formula(),totaldensity, energybins)
            
            
            for meanenergy_bin, NperAngstromperdphidth, attenlength_thisfluorline in zip(energybins,NemitperAngstrompersolidangle,attenfactorinmm_list_fluor):
                #Calc effective length of material (in Angstrom) for each energy and each fluorescence bin
                
                effectivelength_inAngstrom = AttenuationFactorAngledBeam2_Fluorescence(attenlength0*1000,attenlength_thisfluorline*1000,thicknessinmm*1000,Alpharadians, self.TwoTheta_rad_9x1[4],self.Phi_Rad_9x1[4],gamma,beta,boolmask_infattenuation)
                
                NEmittedFromSample_thisline_perpix = effectivelength_inAngstrom * NperAngstromperdphidth *solidanglemat
                
                self.FluorescentImage.append(NEmittedFromSample_thisline_perpix)
                self.FluorLineEnergybins.append(meanenergy_bin)
        self.result_label.setText('Done') 
                
            
        # print(np.shape(self.FluorescentImage))
        # print(np.shape(self.FluorLineEnergybins))
        # FluorLineEnergybins, FluorescentImage
                
                    
        # edat,Ndat = sum_y_for_same_x(emissionenergies,Nemitted)
        
                
        # energybinstot,Ndattot =  sum_y_within_maxsep(emissionenergies_tot,Nemitted_total,1000, 1000)
        

        # totalinput_N = np.trapezoid(y/x,x)/self.bunchfreq * JouleIneV
        # totalinput_E = np.trapezoid(y,x)/self.bunchfreq*1e6 # in µJ
        # totalfluor_N = np.sum(np.array(Ndattot))
        # totalfluor_E = np.sum(np.array(energybinstot)*np.array(Ndattot))/JouleIneV*1e6 # in µJ
        # fluor_energyfraction = totalfluor_E/totalinput_E

        
    def displayatompositions(self):
        atom_data = []
        for atom in self.structure_data:
            atom_type = atom.symbol
            fractional_coords = atom.position / self.structure_data.get_cell().lengths()  # Convert to fractional coordinates
            atom_data.append((atom_type, fractional_coords.tolist()))
        for atom in atom_data:
                print(f"Atom: {atom[0]}, Fractional Coordinates: {atom[1]}")
        self.show_atom_3d_viewer()
        
    def show_atom_3d_viewer(self):
        """
        Show the 3D viewer for atoms in the unit cell.

        """
    
        # Extract atom data
        positions = self.structure_data.get_positions()  # Cartesian positions
        symbols = self.structure_data.get_chemical_symbols()  # Atom types
        cell = self.structure_data.cell  # Unit cell vectors (3x3 matrix)
    
        # Combine positions and symbols
        atoms = [(x, y, z, symbol) for (x, y, z), symbol in zip(positions, symbols)]
    
        # Start the viewer
        # app = QApplication(sys.argv)
        self.viewer = Atom3DViewer(atoms, cell)
        self.viewer.show()
        # sys.exit(app.exec_())    

    

    def calculate_and_plot_1D(self):
        if self.SpectrumManager is None:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified') 
            return
        if not self.cif_file or not self.structure_data:
            self.result_label.setText("Please select a valid CIF file.")
            return
    
        ax = self.ax_other
    
    
        # try:
        #     energy = float(self.energy_input.text())  # Energy in eV
        # except ValueError:
        #     self.result_label.setText("Please enter a valid energy value.")
        #     return
        self.getspectrum()
        
        if self.spec_EeV is None:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified') 
            return
        ax.clear()
        #find spectral peaks
        x = self.spec_EeV
        y = self.spec_WpereV
        dx = x[1]-x[0]
        self.num_peaks, self.peak_locs = find_peaks_in_data(x, y, height=0.00001, distance=int(4000/dx))
        highestEpeak = np.max(self.peak_locs)
        
        # Callback to update progress bar
        def update_progress(current, total):
            progress = int((current / total) * 100)
            # print(progress)
            self.progress_bar.setValue(progress)
            QApplication.processEvents() 
        
        for energy in self.peak_locs:
            theta_max = 45.0  # Maximum two-theta in degrees
            wavelength = hc / energy  # Energy (eV) to wavelength (Å)
            self.progress_bar.setValue(0)
    
            # Compute structure factors
            hklvals, hkl_angles, structure_factors,multiplicities,elements,f0_values = compute_structure_factors(self.structure_data, theta_max, energy, update_progress)
    
    
            for hkl, th,F,m in zip(hklvals, hkl_angles,structure_factors,multiplicities):
            
                h,k,l = hkl
                dspacing  = compute_d_spacing(self.structure_data, h, k, l)
                print(f'hkl: {hkl}, d: {dspacing:.3f}, 2th:{th:.3f},m: {m}, F:{F:.3e}, F**2:{np.square(F):.3e}')

            # sort them and plot 
            combined = list(zip(hkl_angles, structure_factors, hklvals,multiplicities))
            combined_sorted = sorted(combined, key=lambda x: x[0])
            two_theta_sorted, structure_factors_sorted, hkl_sorted,multiplicities_sorted= zip(*combined_sorted)
            two_theta_sorted = np.squeeze(np.array(two_theta_sorted))
            multiplicities_sorted = np.squeeze(np.array(multiplicities_sorted))
            structure_factors_sorted = np.squeeze(np.array(structure_factors_sorted))
            hkl_sorted = np.squeeze(np.array(hkl_sorted))
            
            XRDsignalstrengthPerLength=_calcXRDsignalstrength(hkl_sorted, two_theta_sorted,structure_factors_sorted,multiplicities_sorted,1)

            instrfwhmindegrees = 180/np.pi * np.arctan(float(self.InstrFWHM_input.text())/float(self.SampleToDetDistance_input.text()))
    
            angles,xrdplot=self.makeXRDplot(two_theta_sorted,XRDsignalstrengthPerLength,instrfwhmindegrees,energy,minangle=1)
            # angles,xrdplot_tot=self.makeXRDplot(two_theta_sorted,XRDsignalstrength_Total,instrfwhmindegrees,energy)
            xrdplot = xrdplot/np.max(xrdplot)
            # xrdplot_tot = xrdplot_tot/np.max(xrdplot_tot)
            
            
            # self.ax.plot(two_theta_sorted, structure_factors_sorted*300000, 'o')
            ax.plot(angles, xrdplot, '-')
            # self.ax.plot(angles, xrdplot_tot, '-')
            
        ax.set_title(f'XRD Intensity vs 2Θ: {self.structure_data.get_chemical_formula()}')
        ax.set_xlabel("Two-Theta (degrees)")
        ax.set_ylabel("XRD Intensity")
        ax.grid(True)
        self.plot_canvas_other.draw()
        
        self.plottabs.setCurrentWidget(self.OtherPlots_Tab_GroupBox)

        # self.result_label.setText("Plot updated successfully!")
        # except Exception as e:
        #     self.result_label.setText(f"Error during calculation: {e}")
        
    def calculate_and_plot_2D(self):
        if self.SpectrumManager is None:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified') 
            return
        self.result_label.setText('')
        self.spec_EeV=None
        self.spec_WpereV=None
        volume = self.structure_data.get_volume()  # unit cell volume in Å³
        if self.spec_EeV is None:
            self.getspectrum()
            
        if self.spec_EeV is None:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified') 
            return

        t0 = time()

        #find spectral peaks
        x = self.spec_EeV
        y = self.spec_WpereV
        dx = x[1]-x[0]
        self.num_peaks, self.peak_locs = find_peaks_in_data(x, y, height=0.001, distance=int(4000/dx))

        # Define EnergyBetweenPeaks for harmonic width calculation
        if self.num_peaks > 1:
            EnergyBetweenPeaks = abs(self.peak_locs[1]-self.peak_locs[0])
        else: 
            EnergyBetweenPeaks = (np.max(x)-np.min(x))
        highestEpeak = np.max(self.peak_locs)
        density = calculate_density(self.structure_data)
        self.attenfactorinmm_list = get_attenlengthinmm(self.structure_data.get_chemical_formula(),density, self.peak_locs)
        # print(self.attenfactorinmm_list)
        
        
        if not self.cif_file or not self.structure_data:
            self.result_label.setText("Please select a valid CIF file.")
            return
    
        try:
            Mx= 2*int(self.MHorizontal_input.text())+1
            Ny = 2*int(self.NVertical_input.text())+1
            pixelsizeinmm = float(self.pixelsize_input.text())*0.001
            ycenterpix = int(self.ycenterpix_input.text())
            xcenterpix = int(self.xcenterpix_input.text())
            detdistinmm = float(self.SampleToDetDistance_input.text())
            # energy = float(self.energy_input.text())
            thicknessinmm = float(self.thickness_input.text())
            Alpharadians =  float(self.sampleangle_input.text())*np.pi/180
            if self.onlyfluorescence.checkState()==2:
                onlyfluor = 1
            else:
                onlyfluor = 0
            
        except ValueError:
            self.result_label.setText("Please enter a valid energy value.")
            return
    
        
        # Set # Maximum two-theta in degrees
        fourcorners =np.array([(0,0),(Ny-1,0),(0,Mx-1),(Ny-1,Mx-1)])
        center = np.array([(ycenterpix-1,xcenterpix-1)])
        cornerdists=[]
        for corner in fourcorners:
            distvec = corner-center
            cornerdists.append(np.sqrt(np.sum(distvec**2)))
        maxdisttocorner = np.max(cornerdists)* pixelsizeinmm # mm
        theta_max = np.arctan(maxdisttocorner/detdistinmm)*180/np.pi
        
        # wavelength = hc / energy  # Energy (eV) to wavelength (Å)
        self.progress_bar.setValue(0)

        # Callback to update progress bar
        def update_progress(current, total):
            progress = int((current / total) * 100)
            # print(progress)
            self.progress_bar.setValue(progress)
            QApplication.processEvents() 

            

        # Compute initial list of structure factors, simply to collect non-zero reflections and sort
        hklvals, hkl_angles, structure_factors,multiplicities,elements,f0_values = self.compute_structure_factors(self.structure_data, theta_max, highestEpeak, update_progress)
        
        if len(hklvals)==0:
            self.reportmessage('No reflections within the field of view')
            return
        # print(hklvals)

                    
        # sort them and plot 
        combined = list(zip(hkl_angles, structure_factors, hklvals,multiplicities))
        combined_sorted = sorted(combined, key=lambda x: x[0])
        two_theta_sorted, structure_factors_sorted, hkl_sorted,multiplicities_sorted= zip(*combined_sorted)
        # print(np.shape(np.array(two_theta_sorted)))
        two_theta_sorted = np.squeeze(np.array(two_theta_sorted))
        multiplicities_sorted = np.squeeze(np.array(multiplicities_sorted))
        structure_factors_sorted = np.squeeze(np.array(structure_factors_sorted))
        hkl_sorted = np.squeeze(np.array(hkl_sorted))
        # print(np.shape(hkl_sorted))
       
        # Filter reflections
        self.validreflections = []        
        
        if len(np.shape(hkl_sorted))==1: # if it got turned into a (3,) array instead of a (1,3) array
            two_theta_sorted = [two_theta_sorted]
            multiplicities_sorted = [multiplicities_sorted]
            structure_factors_sorted = [structure_factors_sorted]
            hkl_sorted = [hkl_sorted]
        
        for hkl, th,F,m in zip(hkl_sorted, two_theta_sorted,structure_factors_sorted,multiplicities_sorted):
            if F>1:
                h,k,l = hkl
                
                dspacing  = compute_d_spacing(self.structure_data, h, k, l)
                # print(f'hkl: {hkl}, d: {dspacing:.3f}, 2th:{th:.3f},m: {m: 3d}, F**2:{np.square(F):6.0f}')
                self.validreflections.append({'hkl':hkl,'twotheta':th,'F':F,'multiplicity':m,'dspacing':dspacing})
                



        self.TwoTheta_rad, self.Phi_Rad, self.rmat = self.calcthetaphi_eachPixel_1x1(Mx,Ny,pixelsizeinmm, xcenterpix, ycenterpix, detdistinmm)
        self.dPhi = pixelsizeinmm/self.rmat
        self.dPhi[self.dPhi>2*np.pi]=2*np.pi
        self.dtwotheta = pixelsizeinmm/detdistinmm/(1+np.square(self.rmat/detdistinmm))  
        self.maxtwotheta = np.max(self.TwoTheta_rad)
        
        self.Ntotimage = np.zeros((self.num_peaks,Ny,Mx))
        
        if self.polarizedcheck.checkState()==2:
            polarized = 1 
        else:
            polarized = 0
            
        LPFactor = calcLPFactor(self.TwoTheta_rad,PolarizedFlag=polarized,phi_angle_rad=self.Phi_Rad) 
        gamma,beta,boolmask_infattenuation = AttenuationAngleFactors(Alpharadians, self.TwoTheta_rad, self.Phi_Rad)
        attenuation_params = (gamma,beta,boolmask_infattenuation)

        self.progress_bar.setValue(0)
        totalrefldone = 0
        totvalidrefl = len(self.validreflections)
        
        if not onlyfluor:
            
            # Prepare Bdict and enableDW for threading
            if self.DebyeWallerCheckbox.checkState()==2:
                enableDW = 1
            else:
                enableDW = 0
                
            el_list = []
            B_list = []
            for atomtype in elements: 
                el_list.append(atomtype)
                for DWatom in self.atomDWstructure: 
                    if atomtype==DWatom['atom']:
                        B = float(DWatom['editbox'].text())
                        B_list.append(B)
            Bdict = {'el_list':el_list,'B_list':B_list}

            # Prepare common params for energy calculation
            energy_params = {
                'peak_locs': self.peak_locs,
                'attenfactorinmm_list': self.attenfactorinmm_list,
                'x': x,
                'y': y,
                'EnergyBetweenPeaks': EnergyBetweenPeaks
            }

            # Strategy: Parallelize over Energies (Harmonics).
            # This allows computing the Attenuation Map (heavy) only ONCE per energy thread.
            t0_para = perf_counter()
            with concurrent.futures.ThreadPoolExecutor() as executor:
                total_energies_done = 0

                # Better loop to track indices:
                future_to_idx = {
                    executor.submit(
                        self.calc_single_energy_contribution, 
                        energy_idx, 
                        energy_params, 
                        self.validreflections, 
                        elements, f0_values, pixelsizeinmm, detdistinmm, thicknessinmm, Alpharadians, 
                        attenuation_params, LPFactor, volume, Bdict, enableDW
                    ): energy_idx for energy_idx in range(self.num_peaks)
                }

                for future in concurrent.futures.as_completed(future_to_idx):
                    idx = future_to_idx[future]
                    try:
                        energy_image = future.result()
                        self.Ntotimage[idx] += energy_image
                        
                        total_energies_done += 1
                        # Update progress based on energies instead of reflections?
                        # Or just pulse it?
                        update_progress(total_energies_done, self.num_peaks)
                        
                    except Exception as e:
                        print(f"Error in energy thread {idx}: {e}")
            t1_para = perf_counter()
            print(f"Parallel Section Time: {t1_para-t0_para:.4f} s")
            
        if self.enableFluorescence.checkState()==2 or self.onlyfluorescence.checkState()==2:
            self.calcfluorescenceimage()
        # imgsum = np.sum(self.Ntotimage,0)
        
        
        self.noiseimage = self.generateXRDImageWithNoise()
        self.plot_image(self.noiseimage)
        self.plottabs.setCurrentWidget(self.XRDImageTab_GroupBox)
        # print(np.mean(imgsum))
        # 
        # for ij in range(4):
        #     self.showimageinnewwindow(self.Ntotimage[ij])
        t1 = time()
        print(f"Time to generate XRD image: {t1-t0} seconds")
        
    def generateXRDImageWithNoise(self):

        
        try:
            thickness = float(self.ScintThickness_input.text())
        except ValueError:
            print('Bad Thickness Value')
            return
        
        #make list of QE and DQE vs energy        
        QE_list = []
        DQE_list = []
        for energy in self.peak_locs:
            QE = self.XRDDetector.QEfunc(
                    mat_thickness_inmm=thickness,
                    mat_density=self.XRDDetector.QE_materialdensity,
                    mat_formula=self.XRDDetector.QE_mat,energyineV=energy)
            DQE = self.XRDDetector.DQEfunc(
                    mat_thickness_inmm=thickness,
                    mat_density=self.XRDDetector.QE_materialdensity,
                    mat_formula=self.XRDDetector.QE_mat,energyineV=energy)          
            QE_list.append(QE)
            DQE_list.append(DQE)

            

        if self.enableFluorescence.checkState()==2:
            enablefluor = 1
        else:
            enablefluor = 0
             
        try:

            readnoiseoffset = float(self.readnoiseoffset_input.text())
            readnoiseSTD = float(self.readnoiseSTD_input.text())
            gainADUperkeV = float(self.gainADUperkeV.text())
            InstrFWHM = float(self.InstrFWHM_input.text())
            pixelsize =float(self.pixelsize_input.text())/1000
            if self.onlyfluorescence.checkState()==2:
                onlyfluor = 1
            else:
                onlyfluor = 0
        except ValueError:
            print('Bad Values for Detector Params')
            self.result_label.setText('Bad Values for Detector Params')
            return np.zeros_like(self.Ntotimage[0])
        
        if self.SimNoiseCheckbox.checkState()>0:
        
        
            rng1 = np.random.default_rng(seed)
            # rng.poisson(lam=lam, size=size)
            # AttenuationFactorAngledBeam2_Filter(self,attenlength,thickness, gamma,beta,boolmask_infattenuation)
            
            ADUpereV = gainADUperkeV/1000
            fullimg = np.zeros_like(self.Ntotimage[0])
            rng = np.random.default_rng(seed)
            if not onlyfluor:
                for img,energy,QE,DQE in zip(self.Ntotimage,self.peak_locs,QE_list,DQE_list):
                    # for filt in self.PostSampleFilters_GroupBox.Active
                    if len(self.PostSampleFilters_GroupBox.activeAngles)>0:                        
                        self.img = self.AttenuateNTotImage(img,energy)
                    fullimg += ADUpereV*energy*QE/DQE*rng.poisson(DQE*img) #N photons in this pixel with this particular energy
            if enablefluor or onlyfluor:
                for img,energy,QE,DQE in zip(self.FluorescentImage,self.FluorLineEnergybins,QE_list,DQE_list):
                    if len(self.PostSampleFilters_GroupBox.activeAngles)>0:                        
                        self.img = self.AttenuateNTotImage(img,energy)
                    fullimg += ADUpereV*energy*QE/DQE*rng.poisson(DQE*img) #N photons in this pixel with this particular energy
                
            fullimg = XRDDetector.blurbygaussian(fullimg,InstrFWHM,pixelsize)
            fullimg += readnoiseoffset + rng.poisson(np.ones_like(self.Ntotimage[0])*readnoiseSTD)
        else:
            fullimg = np.zeros_like(self.Ntotimage[0])
            ADUpereV = gainADUperkeV/1000
            if not onlyfluor:
                for img,energy,QE in zip(self.Ntotimage,self.peak_locs,QE_list):
                    if len(self.PostSampleFilters_GroupBox.activeAngles)>0:                        
                        self.img = self.AttenuateNTotImage(img,energy)
                    fullimg += ADUpereV*energy*QE*img #N photons in this pixel with this particular energy
            if enablefluor or onlyfluor:
                for img,energy,QE in zip(self.FluorescentImage,self.FluorLineEnergybins,QE_list):
                    if len(self.PostSampleFilters_GroupBox.activeAngles)>0:                        
                        self.img = self.AttenuateNTotImage(img,energy)
                    fullimg += ADUpereV*energy*QE*img #N photons in this pixel with this particular energy
            fullimg = XRDDetector.blurbygaussian(fullimg,InstrFWHM,pixelsize)
        
        if self.maskfun is not None:
            try:      
                Mx= 2*int(self.MHorizontal_input.text())+1
                Ny = 2*int(self.NVertical_input.text())+1
                # yc = int(self.ycenterpix_input.text())
                # xc = int(self.xcenterpix_input.text())
            except ValueError:
                print('Bad values')
                self.result_label.setText('Bad Values for Detector Params')

            fullimg = self.maskfun(fullimg,Mx,Ny)
        return fullimg  
    
            
    def AttenuateNTotImage(self,img,energy):
        
        
        for chemformula,thickness,Density,Angle_deg in zip(self.PostSampleFilters_GroupBox.activeCompounds,self.PostSampleFilters_GroupBox.activeThicknesses, self.PostSampleFilters_GroupBox.activeDensities, self.PostSampleFilters_GroupBox.activeAngles):
            gamma,beta,boolmask_infattenuation  = AttenuationAngleFactors(np.pi/180*Angle_deg, self.TwoTheta_rad_9x1[4], self.Phi_Rad_9x1[4])
            attenlength=get_attenlengthinmm_single(chemformula, Density, energy)
            img*=AttenuationFactorAngledBeam2_Filter(attenlength,thickness, gamma,beta,boolmask_infattenuation)
        return img
            
        
    def simspon2dintegrate_perpixel(self,TwoTheta_rad_9x1,Phi_Rad_9x1,dtwotheta, dphi,functiontointegrate):
        weights_simpsonOneThird_9x1 = [1,4,1,4,16,4,1,4,1] #weights for simpson's 1/3 rule
        vals_9x1 = functiontointegrate(TwoTheta_rad_9x1,Phi_Rad_9x1)
        pixelvals = np.zeros_like(vals_9x1[0])
        for weight, image_vals in zip(weights_simpsonOneThird_9x1,vals_9x1):
            pixelvals += weight*image_vals*dphi*dtwotheta/9
            
        return pixelvals

        
    def makeXRDplot(self, two_theta_angles,xrdstrength,instrfwhmindegrees,energyineV,minangle=None,maxangle=None):
        if minangle is None: 
            minangle = np.min(two_theta_angles)
        if maxangle is None: 
            maxangle = np.max(two_theta_angles)
        angles = np.linspace(minangle, maxangle,10000)
        xrdplot = np.zeros_like(angles)
        
        try: 
            # crystalsizeinmicrons = float(self.crystalsize_input.text())
            crystalsizeinmicrons=.01
        except ValueError:
            crystalsizeinmicrons = 10
            print('Invalid Crystal Size, Using 10 µm')
        self.result_label.setText(f'Using Characteristic Crystal size: {crystalsizeinmicrons:.3f}µm')
        
        for twoth, xrdI in zip(two_theta_angles,xrdstrength):
            FWHMofeachpeak = FWHMFromInstrumentAndCrystalSize(instrfwhmindegrees, twoth, energyineV, crystalsizeinmicrons)*180/np.pi
            # print(FWHMofeachpeak)
            sig = FWHMofeachpeak/2.355
            xrdplot = xrdplot + xrdI*gaussian(angles,twoth,sig)
            
        return angles,xrdplot

        
        
        

    def plot_image(self,imagedata):
        self.ax.clear()
        self.im = self.ax.imshow(imagedata, cmap="gray")
        self.imgshown = imagedata
        maxval = np.nanmax(self.imgshown)
        minval = np.nanmin(self.imgshown)
        self.vmin_slider.setValue(int(minval/maxval*100))
        # self.vmax_slider.setValue(minval/maxval*100)
        self.update_image_display()
        
    def update_image_display(self):
        if self.imgshown is not None:
            self.aspectratio = 'auto'
            self.ax.set_aspect(self.aspectratio)
            maxval = np.nanmax(self.imgshown)
            minval = np.nanmin(self.imgshown)
            
            self.im.set_clim(vmin=maxval/100*np.min([self.vmin_slider.value(),self.vmax_slider.value()]), vmax=maxval/100*np.max([self.vmin_slider.value(),self.vmax_slider.value()]))
            self.ax.set_aspect(self.aspectratio)
            self.plot_canvas.draw() 
            
    def update_image_display_cake(self):
        if self.cakeimgshown is not None:
            aspectratio = 'auto'
            maxval = np.nanmax(self.cakeimgshown)
            minval = np.nanmin(self.cakeimgshown)
            
            self.im_2d.set_clim(vmin=maxval/100*np.min([self.vmin_slider_cake.value(),self.vmax_slider_cake.value()]), vmax=maxval/100*np.max([self.vmin_slider_cake.value(),self.vmax_slider_cake.value()]))
            self.ax_2d.set_aspect(aspectratio)
            self.plot_canvas_2d.draw() 
    

    
    def getspectrum(self):
        usefile=0
        #If user has defined a spectrum with the DCSPSpectrumTool this session, use that. otherwise use file
        if self.SpectrumManager:

            if len(self.SpectrumManager.spectralpower_dat)>0:
                self.spec_WpereV = self.SpectrumManager.spectralpower_dat
                self.spec_EeV = self.SpectrumManager.energy_dat
                
                if self.SpectrumManager.APSUcheck.checkState()>0:
                    self.bunchfreq = 13e6 #13 MHz in 48-bunch mode
                else:
                    self.bunchfreq = 6.5e6 #6.5 MHz in 24-bunch mode
            else:
                usefile=1
        else: 
            usefile = 1

        if usefile:
            print('No Spectrum Specified')
            self.result_label.setText('No Spectrum Specified')
            return
            # self.readspec()
    def CalcSpectrum(self):
        if self.SpectrumManager is None:
            self.SpectrumManager=DCSSpectrum.DCSSpectrum()
            self.SpectrumManager.Emax_edit.setText(str(30000))
            self.subwindows.append(self.SpectrumManager) 
        self.SpectrumManager.show()

    def LoadSpectrum(self):
        """Load spectrum from CSV or H5 file"""
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(
            self, 
            "Load Spectrum", 
            "", 
            "Spectrum Files (*.csv *.h5);;CSV Files (*.csv);;H5 Files (*.h5);;All Files (*)", 
            options=options
        )
        
        if file_name:
            try:
                # Determine file type from extension
                if file_name.lower().endswith('.h5'):
                    # Read H5 file
                    with h5py.File(file_name, 'r') as h5f:
                        # Check if the expected datasets exist
                        energy_path = '/energy_spectrum/energy_eV'
                        power_path = '/energy_spectrum/spectral_power_W_per_eV'
                        
                        if energy_path not in h5f or power_path not in h5f:
                            QMessageBox.critical(
                                self, 
                                "Error", 
                                f"Invalid H5 format. Expected datasets:\n"
                                f"  {energy_path}\n"
                                f"  {power_path}"
                            )
                            return
                        
                        # Load the datasets
                        energy_dat = np.array(h5f[energy_path])
                        spectralpower_dat = np.array(h5f[power_path])
                        
                elif file_name.lower().endswith('.csv'):
                    # Read the CSV file (two columns: energy, spectral power)
                    data = np.loadtxt(file_name, delimiter=',')
                    
                    if data.shape[1] != 2:
                        QMessageBox.critical(
                            self, 
                            "Error", 
                            "Invalid CSV format. Expected 2 columns: Energy (eV), Spectral Power (W/eV)"
                        )
                        return
                    
                    energy_dat = data[:, 0]
                    spectralpower_dat = data[:, 1]
                else:
                    QMessageBox.critical(
                        self, 
                        "Error", 
                        "Unsupported file format. Please select a .csv or .h5 file."
                    )
                    return
                
                # Create SpectrumManager if it doesn't exist (without showing GUI)
                if self.SpectrumManager is None:
                    self.SpectrumManager = DCSSpectrum.DCSSpectrum()
                    self.SpectrumManager.Emax_edit.setText(str(30000))
                    self.subwindows.append(self.SpectrumManager)
                
                # Set the spectrum data
                self.SpectrumManager.energy_dat = energy_dat
                self.SpectrumManager.spectralpower_dat = spectralpower_dat
                
                # Update the result label
                self.result_label.setText(
                    f"Spectrum loaded: {len(energy_dat)} points, "
                    f"E range: {energy_dat.min():.1f} - {energy_dat.max():.1f} eV"
                )
                
                print(f"Loaded spectrum from: {file_name}")
                print(f"  Energy range: {energy_dat.min():.1f} - {energy_dat.max():.1f} eV")
                print(f"  Number of points: {len(energy_dat)}")
                
                QMessageBox.information(
                    self, 
                    "Success", 
                    f"Spectrum loaded successfully!\n\n"
                    f"Points: {len(energy_dat)}\n"
                    f"Energy range: {energy_dat.min():.1f} - {energy_dat.max():.1f} eV"
                )
                
            except Exception as e:
                QMessageBox.critical(
                    self, 
                    "Error", 
                    f"Failed to load spectrum: {e}"
                )
                print(f"Error loading spectrum: {e}")
        
        
        
    # def calchklXRDintensity(self,twotheta_deg, )
        
    def get_LPFactor(self, twotheta_deg,phi_deg,polarized):
        '''
        I've defined Ez = E0*cos(phi) and Ey = E0*sin(phi), with scattering
        direction in the X-Y plane. And phi is the angle between E0 and the 
        direction perpendicular to the scattering plane
        
        With undulator polarization horizontal, phi should be defined to be 
        zero when the scattering has a vertical deflection, but it's 
        previously defined as zero on horizontal. So, we shift the arguments 
        below to (pi/2 - phi) 
        

        '''
        
        if polarized:
            LP = (np.square(np.pi/2 - np.sin((phi_deg*np.pi/180))) + np.square(np.cos(np.pi/2 - (phi_deg*np.pi/180)) * np.cos(twotheta_deg*np.pi/180))) / 2*np.sin(twotheta_deg/2*np.pi/180)
        else:
            LP = (1 + np.square(np.cos(twotheta_deg*np.pi/180))) / 2*np.sin(twotheta_deg/2*np.pi/180)
        return LP
            
        
    def compute_structure_factors(self, structure, theta_max, energy, progress_callback=None):
        wavelength = hc/energy
        dmax = max(main_window.structure_data.get_cell().lengths())
        """Compute structure factors for all (h, k, l) up to given two-theta at a specific energy."""
        max_index = int(2 * dmax / wavelength * np.sin(theta_max*np.pi/180/2))
                        
        
        unique_hkl_list = generate_unique_hkl(structure, max_index)  # Get unique (h, k, l)
        
        
        if self.DebyeWallerCheckbox.checkState()==2:
            enableDW = 1
        else:
            enableDW = 0
        
        hklvals = []
        angles = []
        structure_factors = []
        multiplicities = []
        total_steps = len(unique_hkl_list)
        elements = set(atom.symbol for atom in structure)  # Extract unique elements
        f0_values = precompute_f0(xraylib, elements, x_min=0.0, x_max=8.0, num_points=100)
        
        #Make dict of debye waller factors
        el_list = []
        B_list = []
        for atomtype in elements: 
            el_list.append(atomtype)
            for DWatom in self.atomDWstructure: #find matchin element in GUI
                if atomtype==DWatom['atom']:
                    B = float(DWatom['editbox'].text())
                    B_list.append(B)
        Bdict = {'el_list':el_list,'B_list':B_list}
            
        refl = 0
        for i, elem in enumerate(unique_hkl_list):
            h,k,l = elem['hkl']
            multiplicity=elem['multiplicity']
            d_spacing = compute_d_spacing(structure, h, k, l)
            if d_spacing > 0 and d_spacing > wavelength / 2:  # Ensure d_spacing is valid
                # print(f'2,i,h,k,l, 2')
                two_theta = 2 * np.degrees(np.arcsin(wavelength / (2 * d_spacing)))
                if two_theta <= theta_max:
                    # print(f'1,{i},{h},{k},{l}, 3')
                    if enableDW:
                        F_hkl = calculate_structure_factor_withDebyeWaller(structure, h, k, l, energy, f0_values,Bdict)  
                    else:
                        F_hkl = calculate_structure_factor(structure, h, k, l, energy, f0_values)  
                    
                    # if h==-1 and k==-1 and l==-1:
                    #     print(h,k,l)
                    #     print(np.abs(F_hkl))
                    

                    # if F_hkl>20:
                        # refl+=1
                        # print(f'{refl}: d: {d_spacing:2f}, th: {two_theta:2f}, {h},{k},{l}: F^2={F_hkl**2}, m = {multiplicity}, m*F = {multiplicity*abs(F_hkl)}')

                    angles.append(two_theta)
                    structure_factors.append(np.abs(F_hkl))
                    multiplicities.append(multiplicity)
                    hklvals.append(elem['hkl'])

            # Update progress
            if progress_callback:
                progress_callback(i + 1, total_steps)
                


        return hklvals,angles, structure_factors,multiplicities,elements,f0_values            
        
        
        
    def loaddetector(self,XRDDetObj):
        self.pixelsize_input.setText(f'{XRDDetObj.pixelsize:.3f}')
        self.MHorizontal_input.setText(str(XRDDetObj.Mxpixels_half))
        self.NVertical_input.setText(str(XRDDetObj.Nypixels_half))
        self.InstrFWHM_input.setText(str(XRDDetObj.InstrFWHM))
        self.ScintThickness_input.setText(str(XRDDetObj.QE_matthickness))
        self.readnoiseoffset_input.setText(str(XRDDetObj.ReadNoiseOffset))
        self.readnoiseSTD_input.setText(str(XRDDetObj.ReadNoiseSTDDev))
        self.gainADUperkeV.setText(f'{XRDDetObj.GainADUperkeV:.3f}')
        
        if hasattr(XRDDetObj,'maskfunction'):
            self.maskfun = XRDDetObj.maskfunction
        else:
            self.maskfun = None
            
        self.DetectorGroupBox.setTitle(f'Detector Params: {XRDDetObj.name}')
        
        
            
        
    def set_cod_sql_preference(self):
        self.use_cod_sql = self.use_cod_sql_action.isChecked()
        self.settings.setValue("CIF_Sources/COD_UseSQL", self.use_cod_sql)
        
    def select_mp_api_key_file(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Select API Key File", "", "Text Files (*.txt);;All Files (*)", options=options)
        if file_name:
            self.settings.setValue("CIF_Sources/MP_API_Key_File_Path", file_name)
            self.api_key_path = file_name
            
            # Validate and update UI
            is_valid = self.validate_mp_api_key(file_name)
            self.MP_APIkeyValid = is_valid
            self.matproject_checkbox.setEnabled(is_valid)
            if is_valid:
                 self.matproject_checkbox.setToolTip("Fetch from Materials Project")
                 self.matproject_checkbox.setChecked(False) 
                 QMessageBox.information(self, "Success", "API Key is valid and has been saved.")
            else:
                 self.matproject_checkbox.setToolTip("Valid API Key not found. Please configure in 'CIF Sources' menu.")
                 self.matproject_checkbox.setChecked(False)
                 QMessageBox.warning(self, "Invalid Key", "The selected file does not contain a valid API key or the test query failed.")

    def get_mp_api_key(self):
        try:
             with open(self.api_key_path, 'r') as file:
                content = file.readline().strip()
                if not content.isalnum():
                     return None
                return content
        except Exception:
            return None

    def validate_mp_api_key(self, api_key_path):
        try:
            with open(api_key_path, 'r') as file:
                content = file.readline().strip()
                if not content.isalnum():
                     return False
            
            # Test query
            with MPRester(content) as mpr:
                # Query for a simple material, e.g., Silicon (mp-149)
                mpr.materials.summary.search(material_ids=["mp-149"])
            return True
        except Exception as e:
            print(f"API Key Validation Failed. Add a valid API key using the menu bar. 'CIF Sources' -> 'Set Materials Project API Key'.")
            return False

    def check_mp_api_key(self):
        if self.validate_mp_api_key(self.api_key_path):
             QMessageBox.information(self, "Valid", f"The API key in {self.api_key_path} works correctly.")
        else:
             QMessageBox.critical(self, "Invalid", f"The API key in {self.api_key_path} is invalid or the connection failed.")

    def reportmessage(self,string):
        print(string)
        self.result_label.setText(string)
        
    def closeEvent(self,event):
        for wins in self.subwindows:
            if wins:
                wins.close()    
                
    
    def plotcake(self,radial,cakeintensity):

        self.cakeimgshown = cakeintensity
        self.ax_2d.clear()
        self.im_2d = self.ax_2d.imshow(cakeintensity,extent=[0, 3.1416, radial.max(), radial.min()],cmap='gray')
        self.ax_2d.set_aspect('auto')
        self.ax_2d.set_ylabel(f'2{chr(952)} (deg.)')
        self.plot_canvas_2d.draw()
        
    def write_csv(self,filename, col1, col2):

        if len(col1) != len(col2):
            raise ValueError("Columns must have the same length.")
    
        # Open the file in write mode
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
    
            # Write header
            # writer.writerow(['Energy (eV)', 'Spectral Power (W/eV)'])

            for value1, value2 in zip(col1, col2):
                writer.writerow([value1, value2])
        
    def SaveAzInt1DData(self):
        
        # self.plot1DIntegration(self.twotheta,self.IntvsTh)
        # self.plotcake(self.radial,self.cakeintensity)
        if self.twotheta is not None:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
    
            # If the user selects a file, write the CSV
            if file_name:
                try:
                    self.write_csv(file_name,self.twotheta,self.IntvsTh)
                    QMessageBox.information(self, "Success", f"CSV file saved as: {file_name}")
    
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to save CSV file: {e}")
        else: 
            QMessageBox.critical(self, "Error", "No Data Yet")
            
    def plot1DIntegration(self,twothetavals,Ivals):

        # cmap = plt.get_cmap('rainbow', 4)
        # c = cmap(1)
        self.ax_1d.clear()
        self.ax_1d.plot(twothetavals,Ivals)
        self.ax_1d.set_xlabel(f'2{chr(952)} (deg.)')
        self.plot_canvas_1d.draw()
        
    # def moveGroupBoxToWindow(self,GroupBox):
        
        
    def isvalidangle(self,widget):    
        good=0
        try:
            value = float(widget.text())
            good = 1 if float(value)>=-180 and float(value)<=180 else 0
        except:
            good= 0  
            
        if good==0:    
            blink_red(widget)
            self.reportmessage('Enter Valid Angle -180 to 180 degrees')
            widget.setText('0.0')
        
    def isvalidvalue(self,widget):
        good = 0
        try:
            value = float(widget.text())
            good = 1 if float(value)>=0 else 0
        except:
            good= 0
        if good==0:
            blink_red(widget)
            self.reportmessage('Enter Valid Value > 0')
            widget.setText('???')

    def calc_single_energy_contribution(self, energy_idx, energy_params, reflections_chunk, elements, f0_values, pixelsizeinmm, detdistinmm, thicknessinmm, Alpharadians, attenuation_params, LPFactor, volume, Bdict, enableDW):
        
        # t_start = perf_counter()
        
        peak_energy = energy_params['peak_locs'][energy_idx]
        attenlengthinmm = energy_params['attenfactorinmm_list'][energy_idx]
        x = energy_params['x']
        y = energy_params['y']
        EnergyBetweenPeaks = energy_params['EnergyBetweenPeaks']
        
        gamma, beta, boolmask_infattenuation = attenuation_params
        
        # attenuation map for this energy
        effectivelength_Angstrom_Map = AttenuationFactorAngledBeam2(attenlengthinmm*1000,thicknessinmm*1000,Alpharadians, self.TwoTheta_rad,self.Phi_Rad,gamma,beta,boolmask_infattenuation)
        
        
        minE = peak_energy - EnergyBetweenPeaks/2.1
        maxE = peak_energy + EnergyBetweenPeaks/2.1
        wavelength = hc/peak_energy 
        
        harmonic_E = x[(x>minE)&(x<maxE)] 
        harmonic_dPdE = y[(x>minE)&(x<maxE)] 
        
        scatted_image_sum = np.zeros_like(self.TwoTheta_rad)
        
        # reflections for this energy
        for validreflection in reflections_chunk:
            h,k,l = validreflection['hkl']
            multiplicity=validreflection['multiplicity']
            d_spacing = validreflection['dspacing']
           
            bragg_limit_factor = hc / (2 * d_spacing)
            
            # check if diffraction is possible at all for maxE
            val_min = bragg_limit_factor / maxE
            
            if val_min > 1.0:
                continue # Wavelength too long for this d-spacing
                
            val_max = bragg_limit_factor / minE
            if val_max > 1.0:
                val_max = 1.0
                
            min_2th_rad = 2 * np.arcsin(val_min)
            max_2th_rad = 2 * np.arcsin(val_max)
            
            # only process if some part of the ring is on the detector (within max theta)
            if min_2th_rad < self.maxtwotheta:
                
                # mask within the valid 2-theta range for this harmonic
                # add a small buffer 
                mask = (self.TwoTheta_rad >= min_2th_rad) & (self.TwoTheta_rad <= max_2th_rad)
                
                if not np.any(mask):
                    continue
                    
                if enableDW:
                    F_hkl = calculate_structure_factor_withDebyeWaller(self.structure_data, h, k, l, peak_energy, f0_values,Bdict)  
                else:
                    F_hkl = calculate_structure_factor(self.structure_data, h, k, l, peak_energy, f0_values)      
                
                G = (r0**2)*(wavelength**3)*multiplicity*(F_hkl**2)/4/(volume**2) 
                
                Wper2th_rad_interp, NperSecPer2th_rad_interp = convertspectrumRangeto2th(harmonic_E,harmonic_dPdE,d_spacing)
                
                # Apply calculations only on the masked pixels
                TwoTheta_masked = self.TwoTheta_rad[mask]
                
                NperSecper2th_rad_masked = NperSecPer2th_rad_interp(TwoTheta_masked)
                
                Nper2thetaperphi_rad_masked = NperSecper2th_rad_masked / 2/np.pi / self.bunchfreq * JouleIneV
                
                # element-wise operations on masked arrays
                ScatteredImageSingle_masked = (
                    Nper2thetaperphi_rad_masked * 
                    G * 
                    LPFactor[mask] * 
                    effectivelength_Angstrom_Map[mask] * 
                    self.dtwotheta[mask] * 
                    self.dPhi[mask]
                )
                
                ScatteredImageSingle_masked[ScatteredImageSingle_masked<0] = 0
                
                # accumulate
                scatted_image_sum[mask] += ScatteredImageSingle_masked
        
        return scatted_image_sum

                
    def calcXRDimage(self, validreflection,elements,f0_values,pixelsizeinmm,detdistinmm,thicknessinmm,Alpharadians, attenuation_params, LPFactor, volume, Bdict, enableDW):
        
        '''
        Strategy:
        0) separate spectrum into separate harmonics
        1) calc FHKL for the peak of a harmonic, then treat the entire harmonic range as having constant FHKL, attenuation, etc
            this is valid away from abs. edges
        '''
        
        # , self.Phi_Rad_9x1
        h,k,l = validreflection['hkl']
        multiplicity=validreflection['multiplicity']
        d_spacing =validreflection['dspacing']
        
        gamma,beta,boolmask_infattenuation = attenuation_params

        #find spectral peaks
        x = self.spec_EeV
        y = self.spec_WpereV

        if self.num_peaks>1:
            EnergyBetweenPeaks = abs(self.peak_locs[1]-self.peak_locs[0])
        else: 
            EnergyBetweenPeaks = (np.max(x)-np.min(x))
        
        ScatteredImage = []
        
        # gamma,beta,boolmask_infattenuation  = AttenuationAngleFactors(Alpharadians, self.TwoTheta_rad_9x1[4], self.Phi_Rad_9x1[4])
        
        for peak_energy,attenlengthinmm in zip(self.peak_locs,self.attenfactorinmm_list):
            
            minE = peak_energy - EnergyBetweenPeaks/2.1
            maxE = peak_energy + EnergyBetweenPeaks/2.1
            wavelength = hc/peak_energy #approximation, can substitute a matrix with lam = 2*d*sin(2th/2) later
            
            min_2th_rad = 2*np.arcsin(hc/2/d_spacing/maxE)
            
            if min_2th_rad<self.maxtwotheta:
                harmonic_E = x[(x>minE)&(x<maxE)] 
                harmonic_dPdE = y[(x>minE)&(x<maxE)] 
                effectivelength_Angstrom = AttenuationFactorAngledBeam2(attenlengthinmm*1000,thicknessinmm*1000,Alpharadians, self.TwoTheta_rad,self.Phi_Rad,gamma,beta,boolmask_infattenuation)
               
                if enableDW:
                    F_hkl = calculate_structure_factor_withDebyeWaller(self.structure_data, h, k, l, peak_energy, f0_values,Bdict)  
                else:
                    F_hkl = calculate_structure_factor(self.structure_data, h, k, l, peak_energy, f0_values)      
                
                G = (r0**2)*(wavelength**3)*multiplicity*(F_hkl**2)/4/(volume**2) # [A^2]*[A^3]/[A^6] = 1/Angstrom
                expectedatten = np.exp(-thicknessinmm/attenlengthinmm)

                Wper2th_rad_interp, NperSecPer2th_rad_interp= convertspectrumRangeto2th(harmonic_E,harmonic_dPdE,d_spacing)
                
                NperSecper2th_rad = NperSecPer2th_rad_interp(self.TwoTheta_rad)
                
                # Use local variable instead of self.Nper2thetaperphi_rad
                Nper2thetaperphi_rad = NperSecper2th_rad / 2/np.pi / self.bunchfreq * JouleIneV
                
                ScatteredImageSingle = Nper2thetaperphi_rad * G * LPFactor * effectivelength_Angstrom *self.dtwotheta * self.dPhi
                ScatteredImageSingle[ScatteredImageSingle<0] = 0
                
                ScatteredImage.append(ScatteredImageSingle)
                print(f'E: {peak_energy:.3f}, [{h:> 2d},{k:> 2d},{l:> 2d}], F: {F_hkl:.3f}, Exp. Trans: {expectedatten:.3f}, ringsum: {np.sum(ScatteredImageSingle):.3e}')
                
            else:
                print(f'{peak_energy:.3f}: entire harmonic range off image')
                ScatteredImage.append(np.zeros_like(self.rmat))
        
        return ScatteredImage  
     
def AziIntegrate2D(img = None,d = None,x0 = None,y0 = None,energy = None,pixsize = None):
    if img is None:
        return
    tiltval = 0
    tiltplanerotation = 0
    lam = hc/energy * 1e-10
    pyFAI_orientation=3
    
    detcustom = pyFAI.detectors.Detector(pixel1=pixsize/1e3, pixel2=pixsize/1e3,orientation=pyFAI_orientation)
    ai = AzimuthalIntegrator(detector=detcustom,wavelength=lam)
    ai.setFit2D(d, x0, y0,tilt=tiltval,tiltPlanRotation=tiltplanerotation)
    res2d = ai.integrate2d(img, 1000,360,unit="2th_deg")
    radial = res2d.radial;
    cakeintensity = res2d.intensity.T
    
    return radial, cakeintensity

def AziIntegrate1D_single(img = None,d = None,x0 = None,y0 = None,energy = None,pixsize = None):
    if img is None:
        return
    tiltval = 0
    tiltplanerotation = 0
    lam = hc/energy * 1e-10
    pyFAI_orientation=3
    detcustom = pyFAI.detectors.Detector(pixel1=pixsize/1e3, pixel2=pixsize/1e3,orientation=pyFAI_orientation)
    ai = AzimuthalIntegrator(detector=detcustom,wavelength=lam)
    ai.setFit2D(d, x0, y0,tilt=tiltval,tiltPlanRotation=tiltplanerotation)
    twotheta,Ival = ai.integrate1d(img, 1000,unit="2th_deg")
    twothetavals = np.array(twotheta)
    Ivals = np.array(Ival)
    return twothetavals,Ivals
    
def get_nested_attribute(obj, attr_path):
    """
    Retrieves the value of a nested attribute from an object using a dot-separated attribute path.

    """
    try:
        # split the attribute path by dots and navigate the attributes
        for attr in attr_path.split('.'):
            obj = getattr(obj, attr)
        return obj
    except AttributeError as e:
        raise AttributeError(f"Error accessing attribute path '{attr_path}': {e}")



    def convertspectrumto2th(self,dspacing):
        if self.spec_EeV is None:
            self.getspectrum()
            
            print(dspacing, self.spec_EeV[0])

        hc = 12398.4198 # electron volts * angstrom    
        def convertEto2theta(Energies):

            return 2 * np.arcsin(hc/2/Energies/dspacing) #in rad
        #hc is in electron volts * angstrom
        # dspacing is in Angstrom 
        
        def convertdPdEtodPdtheta(twothetas,PowerPereV):
            #Scale Power Per eV by width of a 2theta degree 
            dEd2th = hc/4/dspacing/np.sin(twothetas/2)/np.tan(twothetas/2)
            return PowerPereV*dEd2th

        twothetas_spectrum = convertEto2theta(self.spec_EeV)
        Wper2thetadegree = convertdPdEtodPdtheta(twothetas_spectrum,self.spec_WpereV)
        twothetas_spectrum = np.flip(twothetas_spectrum)
        Wper2thetadegree = np.flip(Wper2thetadegree)
        
        return twothetas_spectrum, Wper2thetadegree
    
    
        

class MaterialTableWidget_MatProject(QWidget):
    def __init__(self, parentobj, matlist, attributes, attributelabels,selectmethod_str='select_material_forDebyeTemp'):
        """
        This widget is used to display a list of materials from the Materials Project.

        Parameters:
            parentobj: The parent object whose selectedmat will be set.
            matlist: List of objects to display in the table.
            attributes: List of attribute names (strings) to display for each object.
        """
        super().__init__()
        self.parentobj = parentobj
        self.matlist = matlist
        self.attributes = attributes
        self.setGeometry(25, 50, 600, 600)
        self.setWindowTitle('Materials From the Materials Project')
        # Main layout
        layout = QVBoxLayout(self)

        # Create the table widget
        self.table = QTableWidget(self)
        self.table.setColumnCount(len(attributes) + 1)  # +1 for the 'Select' button
        self.table.setHorizontalHeaderLabels(['Select'] + attributelabels)
        self.table.setRowCount(len(matlist))


        func = getattr(self, selectmethod_str)
        # Populate the table
        for row, material in enumerate(matlist):
            # Add a button in the first column
            select_button = QPushButton("Select")
            select_button.clicked.connect(lambda _, r=row: func(r))
            self.table.setCellWidget(row, 0, select_button)

            # Add attributes to the table
            for col, attr_name in enumerate(attributes, start=1):
                attr_value = get_nested_attribute(material, attr_name)
                if isinstance(attr_value, float):
                    valuestr = f'{attr_value:0.3f}'
                else:
                    valuestr = f'{attr_value}'
                self.table.setItem(row, col, QTableWidgetItem(valuestr))

        # Add the table to the layout
        layout.addWidget(self.table)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    def getDebyeTemp(self,material_id):
        with MPRester(self.parentobj.API_KEY) as mpr:
            try:
                elasticitydata = mpr.materials.elasticity.search(material_ids=material_id)  
            except (ConnectionError,OSError):
                print('Connection Error. Perhaps Mat Project is down?')
                return
        # elastic_tensor = elasticitydata.elastic_tensor
        
        if elasticitydata==[]:
            self.parentobj.reportmessage('No Elastic Tensor Data, D-W Factors not written')
            return
        DbT = elasticitydata[0].debye_temperature
            
        if DbT is not None:
            self.parentobj.debye_temperature = DbT
            if isinstance(DbT, float):
                valuestr = f'{DbT:0.3f}'
            else:
                valuestr = f'{DbT}'
            print(f'Debye Temp: {valuestr}')
            self.parentobj.result_label.setText(f'Debye Temp: {valuestr}')
            T = 300 #K
            # B = 12*h2*T/m
            
            elements = self.parentobj.structure_data.get_chemical_symbols()
            compstr = xraylib.CompoundParser(''.join(elements))
            
            
            for Z in compstr['Elements']:
                # self.parentobj.Z = Z
                elname = xraylib.AtomicNumberToSymbol(Z)
                weightpermol = xraylib.AtomicWeight(Z)
                B = calcDWfromDebyeTemp(T, DbT,weightpermol)
                print(f'{elname}: weight: {weightpermol}, B: {B}')
                
                for DWatom in self.parentobj.atomDWstructure:

                    if elname==DWatom['atom']:
                        DWatom['editbox'].setText(f'{B:0.3f}')
                        
            self.parentobj.reportmessage('D-W Factors estimated based on Elasticity Tensor (assuming cubic crystal)')
        else:
            
            self.parentobj.reportmessage('No Debye Temp.')
            
    

    def select_material_forDebyeTemp(self, row):
        """Handles selecting a material."""
        selected_material = self.matlist[row]
        self.parentobj.selectedmat = selected_material.material_id
        self.getDebyeTemp(selected_material.material_id)
            
    def select_material_forCIF(self, row):
        """Handles selecting a material."""
        selected_material = self.matlist[row]
        self.parentobj.selectedmat = selected_material.material_id
        with MPRester(self.parentobj.API_KEY) as mpr:
            try:
                self.parentobj.data = mpr.materials.summary.get_data_by_id(selected_material.material_id)  
            except (ConnectionError,OSError):
                print('Connection Error. Perhaps Mat Project is down?')
                return
            filen = os.path.join(dir_path,'ciffiles',f'{selected_material.material_id}.cif')
            self.parentobj.data.structure.to(fmt="cif",filename=filen)
            if os.path.exists(filen):
                print(f'File: {filen} written')
                DWfound = self.parentobj.load_cif(filen)
                
                if DWfound==0:
                    self.parentobj.reportmessage('Chemical Not Found in Table of Elemental Material D-W Factors')
                    self.getDebyeTemp(selected_material.material_id)
                else:
                    self.parentobj.reportmessage('D-W Factors from Table of D-W Factors for Elemental Materials (Peng et al.)')
                


class MaterialTableWidget_COD(QWidget):
    def __init__(self, parentobj, dictlist,selectmethod_str='select_material_forCIF'):
        """
        Initialize the MaterialTableWidget_COD.

        Parameters:
            parentobj: The parent object whose selectedmat will be set.
            matlist: List of objects to display in the table.
            attributes: List of attribute names (strings) to display for each object.
        """
        super().__init__()
        self.parentobj = parentobj
        self.dictlist = dictlist
        self.setGeometry(700, 50, 900, 600)
        self.setWindowTitle('Materials From Crystallography Open Database')
        # Main layout
        layout = QVBoxLayout(self)
        

        # Create the table widget
        self.table = QTableWidget(self)
        self.table.setColumnCount(len(dictlist[0].keys()) + 1)  # +1 for the 'Select' button
        self.table.setHorizontalHeaderLabels(['Select'] + list(dictlist[0].keys()))
        self.table.setRowCount(len(dictlist))

        # def makemultilinewidget(list_of_text_strings):
        #     paramwidget = QWidget()
        #     vlayout = QVBoxLayout()
        #     paramwidget.setLayout(vlayout)
        #     for string in list_of_text_strings:
        #         vlayout.addWidget(QLabel(string))            
        #     return paramwidget
            
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        
        bibcol =None
        func = getattr(self, selectmethod_str)
        # Populate the table
        for row, materialdict in enumerate(dictlist):
            # Add a button in the first column
            select_button = QPushButton("Select")
            select_button.clicked.connect(lambda _, r=row: func(r))
            self.table.setCellWidget(row, 0, select_button)
            attributes = list(materialdict.keys())
            # Add attributes to the table
            for col, attr_name in enumerate(attributes, start=1):
                attr_value = materialdict[attr_name]
                if isinstance(attr_value, float):
                    valuestr = f'{attr_value:0.3f}'
                else:
                    valuestr = f'{attr_value}'
                    
                if 'biblio' in attr_name.casefold():
                    # valuestr = textwrap.fill(attr_value,80).strip() #wrap the long bibliography
                    # self.table.setCellWidget(row, col, QLabel(valuestr))
                    bibcol = col
                    
                if 'parameters' in attr_name.casefold():

                    # stringlist = valuestr.replace(';','\n')

                    # paramwidget = makemultilinewidget(stringlist) 
                    # self.table.setCellWidget(row, col, paramwidget)
                    self.table.setCellWidget(row, col, QLabel(valuestr.replace(';','\n').replace(' ','')))
                else:
                    self.table.setItem(row, col, QTableWidgetItem(valuestr))
                
            if bibcol is not None:
                self.table.horizontalHeader().setSectionResizeMode(bibcol, QHeaderView.ResizeMode.ResizeToContents)
        # Add the table to the layout
        layout.addWidget(self.table)
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)

    def select_material_forCIF(self, row):
        """Handles selecting a material."""
        selected_material = self.dictlist[row]
        # print(f'mat: {selected_material}')
        formula = selected_material['Chemical Formula']
        # print(f'selected chem formula: {formula}')
        IDnumber= selected_material['COD ID']
        print(IDnumber)
        self.parentobj.selectedmat = IDnumber
        
        ciffilename = getCODCif(IDnumber,formula)
        if ciffilename is not None:
            print(f'File: {ciffilename} written')
            DWfound = self.parentobj.load_cif(ciffilename)
            
            if DWfound==0:
                self.parentobj.reportmessage('Chemical Not Found in Table of Elemental Material D-W Factors. Enter it manually')

            else:
                self.parentobj.reportmessage('D-W Factors from Table of D-W Factors for Elemental Materials (Peng et al.)')
            


class ExperimentGeometry():
    def __init__(self,filterlist=[],sampleangle_withrespecttobeamaxis_deg=0.0,detectordistance=100):
        self.filterlist = filterlist
        self.sampleangle_deg = sampleangle_withrespecttobeamaxis_deg
        self.detdistance = detectordistance
        
    @staticmethod
    def makeImpactExperiment_Estation():
        exptgeo = ExperimentGeometry()
        exptgeo.filterlist=[
            {'mat':'LiF','density':2.64 ,'thickness':2.0,'angle':62,'zpos':0},
            {'mat':'C16H14O3','density':1.26 ,'thickness':8.0,'angle':-28,'zpos':80},
            ]
        exptgeo.sampleangle_deg = 62
        exptgeo.detdistance = 140
        return exptgeo
    
    @staticmethod
    def makeImpactExperiment_Dstation():
        exptgeo = ExperimentGeometry()
        exptgeo.filterlist=[
            {'mat':'LiF','density':2.64 ,'thickness':2.0,'angle':-62,'zpos':0},
            {'mat':'C16H14O3','density':1.26 ,'thickness':8.0,'angle':28,'zpos':80},
            ]
        exptgeo.sampleangle_deg = -62
        exptgeo.detdistance = 140
        return exptgeo
    
    @staticmethod
    def makeLaserShockExperiment():
        exptgeo = ExperimentGeometry()
        exptgeo.filterlist=[
            {'mat':'LiF','density':2.64 ,'thickness':0.5,'angle':-52,'zpos':0},
            {'mat':'C16H14O3','density':1.26 ,'thickness':8.0,'angle':0,'zpos':80},
            ]
        exptgeo.sampleangle_deg = -52
        exptgeo.detdistance = 100
        return exptgeo
    

    
class NewWindow(QMainWindow):
    def __init__(self, groupbox):
        super().__init__()
        self.setWindowTitle("New Window")
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        layout.addWidget(groupbox)        

class XRDDetector():
    def __init__(self,pixelsize=[],Nypixels_half=[],Mxpixels_half=[],InstrFWHM=[],QE=[],DQE=[],ReadNoiseOffset=[],ReadNoiseSTDDev=[],GainADUperkeV=[]):
            self.pixelsize = pixelsize
            self.Nypixels_half=Nypixels_half
            self.Mxpixels_half=Mxpixels_half
            self.InstrFWHM=InstrFWHM
            self.QE = QE
            self.DQE = DQE
            self.ReadNoiseOffset=ReadNoiseOffset
            self.ReadNoiseSTDDev=ReadNoiseSTDDev
            self.GainADUperkeV = GainADUperkeV
            
    @staticmethod
    def makeGenericSquareDetector():
        detobj = XRDDetector()
        detobj.pixelsize = 100
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 1
        detobj.DQE = 1
        detobj.InstrFWHM=0.01
        detobj.ReadNoiseOffset=10
        detobj.ReadNoiseSTDDev=0.01
        detobj.GainADUperkeV = 1 
        detobj.maskfunction = None
        detobj.name = 'Ideal Square Detector'
        
        def QEfunc(mat_thickness_inmm=None, mat_density=None,
        mat_formula=None,energyineV=None):
            return 1.0
        
        detobj.QE_mat = 'Unobtainium'
        detobj.QE_materialdensity = 999 #measured by Rayonix for settled screen
        detobj.QE_matthickness = 0.01 #mm
        
        detobj.QEfunc = QEfunc
        detobj.DQEfunc = QEfunc
        
        return detobj            
            
    @staticmethod
    def makeGOSRayonix():
        detobj = XRDDetector()
        detobj.pixelsize = 79
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 0.5
        detobj.DQE = detobj.QE
        detobj.InstrFWHM=0.2
        detobj.ReadNoiseOffset=11
        detobj.ReadNoiseSTDDev=2
        detobj.GainADUperkeV = 0.121 #(8 pe-/12 keV)/(5.5 pe/ADU)
        detobj.maskfunction = maskfun_circleinscribedinsquare
        detobj.name = 'GadOx Rayonix SX165'
        detobj.QE_mat = 'Gd2O2S'
        detobj.QE_materialdensity = 4.3 #measured by Rayonix for settled screen
        detobj.QE_matthickness = 0.08 #mm
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = XRDDetector.QEfunc_generic
        return detobj
    
    @staticmethod
    def makeIdealRayonix():
        detobj = XRDDetector()
        detobj.pixelsize = 79
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 1
        detobj.DQE = detobj.QE
        detobj.InstrFWHM=0.01
        detobj.ReadNoiseOffset=11
        detobj.ReadNoiseSTDDev=2
        detobj.GainADUperkeV = 0.121 #(8 pe-/12 keV)/(5.5 pe/ADU)
        detobj.maskfunction = maskfun_circleinscribedinsquare
        detobj.name = 'Ideal Rayonix SX165'
        detobj.QE_mat = 'Gd2O2S'
        detobj.QE_materialdensity = 4.3 #measured by Rayonix for settled screen
        detobj.QE_matthickness = 0.08 #mm
        def QEfunc(mat_thickness_inmm=None, mat_density=None,
        mat_formula=None,energyineV=None):
            return 1.0
        
        detobj.QEfunc = QEfunc
        detobj.DQEfunc = QEfunc
        return detobj
    
    @staticmethod        
    def makeCsIRayonix():
        detobj = XRDDetector()
        detobj.pixelsize = 79
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 0.5
        detobj.DQE = detobj.QE
        detobj.InstrFWHM=0.2
        detobj.ReadNoiseOffset=11
        detobj.ReadNoiseSTDDev=2
        detobj.GainADUperkeV = 0.076 #(5 pe-/12 keV)/(5.5 pe/ADU)      
        detobj.name = 'CsI Rayonix SX165'
        detobj.QE_mat = 'CsI'
        detobj.QE_materialdensity = 0.75*4.51 # estimated to make my 150 µm thick QE calc match that of Rayonix
        detobj.QE_matthickness = 0.300 #mm
        ''' Custom Screen Specifications: 
            Thickness: 300µ+/-30
            Amorphous Carbon Substrate: 169.5mm diameter, Thickness: 0.5mm +/- 0.1mm
            Additional Film coating: 0.025 mm thick
    '''
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = XRDDetector.QEfunc_generic
        
        detobj.maskfunction = maskfun_circleinscribedinsquare
        
        
        return detobj
    
    @staticmethod        
    def makeSiKeckPAD():
        detobj = XRDDetector()
        detobj.pixelsize = 150
        detobj.Nypixels_half=256
        detobj.Mxpixels_half=256
        detobj.QE = 0.2
        detobj.DQE = detobj.QE
        detobj.InstrFWHM=0.01
        detobj.ReadNoiseOffset=11
        detobj.ReadNoiseSTDDev=0.7 * 1.77 #((0.7 8 keV X-ray (high gain)*(1.77 ADU/8keV ph))
        detobj.GainADUperkeV = 0.22152 #(1.77 ADU/8 keV)
        
        
        detobj.name = 'Si Keck PAD'
        detobj.QE_mat = 'Si'
        detobj.QE_materialdensity = 2.33 # estimated to make my 150 µm thick QE calc match that of Rayonix
        detobj.QE_matthickness = 0.500 #mm
        '''  
    '''
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = XRDDetector.QEfunc_generic
        
        detobj.maskfunction = None
               
        return detobj
    
    @staticmethod        
    def makeCdTeKeckPAD():
        detobj = XRDDetector()
        detobj.pixelsize = 150
        detobj.Nypixels_half=256
        detobj.Mxpixels_half=256
        detobj.QE = 0.2
        detobj.DQE = detobj.QE
        detobj.InstrFWHM=0.01
        detobj.ReadNoiseOffset=11
        detobj.ReadNoiseSTDDev=0.7 * 1.77 #((0.7 8 keV X-ray (high gain)*(1.77 ADU/8keV ph))
        detobj.GainADUperkeV = 0.22152 #(1.77 ADU/8 keV)
        
        
        detobj.name = '750µm-CdTe Keck PAD'
        detobj.QE_mat = 'CdTe'
        detobj.QE_materialdensity = 5.85 # estimated to make my 150 µm thick QE calc match that of Rayonix
        detobj.QE_matthickness = 0.750 #mm
        '''  
    '''
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = XRDDetector.QEfunc_generic
        
        detobj.maskfunction = None
               
        return detobj
    
    @staticmethod
    def makeDCSFourFrameDetector_LSO_75mm():
        detobj = XRDDetector()
        detobj.pixelsize = 1000/23
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 0.8
        detobj.DQE = 0.3
        detobj.InstrFWHM=0.22
        detobj.ReadNoiseOffset=600
        detobj.ReadNoiseSTDDev=60
        detobj.GainADUperkeV = 25      
        
        detobj.QE_mat = 'Lu2SiO5'
        detobj.QE_materialdensity = 0.55 * 7.4 #Estimated packing density
        detobj.QE_matthickness = 0.100 #mm
        '''(7.4 (g / cc)) * 0.55 * (0.100 mm) = 40.7 mg / (cm^2)'''

        def DQEfun(**kwargs):
            return 0.3/0.8*XRDDetector.QEfunc_generic(**kwargs)
        
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = DQEfun
        
        def maskfun(img,Mx,Ny):
            maxradius_inpixels = 75*1000/detobj.pixelsize/2
            img = maskfun_generic(img, Mx, Ny, maxradius_inpixels)
            return img
        
        detobj.maskfunction = maskfun
        detobj.name = 'DCS 4-Frame Detector (75mm)'
        return detobj
    @staticmethod        
    def makeDCSFourFrameDetector_LSO_120mm():
        detobj = XRDDetector()
        detobj.pixelsize = 2000/29
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 0.8
        detobj.DQE = 0.14
        detobj.InstrFWHM=0.35
        detobj.ReadNoiseOffset=600
        detobj.ReadNoiseSTDDev=10
        detobj.GainADUperkeV = 25/2.7 
        
        detobj.QE_mat = 'Lu2SiO5'
        detobj.QE_materialdensity = 0.55 * 7.4 #Estimated packing density
        detobj.QE_matthickness = 0.100 #mm
        
        def DQEfun(**kwargs):
            return 0.14/0.8*XRDDetector.QEfunc_generic(**kwargs)
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = DQEfun
        
        def maskfun(img,Mx,Ny):
            maxradius_inpixels = 120*1000/detobj.pixelsize/2
            img = maskfun_generic(img, Mx, Ny, maxradius_inpixels)
            return img
        
        detobj.maskfunction = maskfun
        detobj.name = 'DCS 4-Frame Detector (120mm)'
        
        return detobj
    @staticmethod
    def makeDCSFourFrameDetector_LSO_150mm():
        detobj = XRDDetector()
        detobj.pixelsize = 2000/23
        detobj.Nypixels_half=1024
        detobj.Mxpixels_half=1024
        detobj.QE = 0.8
        detobj.DQE = 0.05
        detobj.InstrFWHM=0.35
        detobj.ReadNoiseOffset=600
        detobj.ReadNoiseSTDDev=10
        detobj.GainADUperkeV = 25/8       
        
        detobj.QE_mat = 'Lu2SiO5'
        detobj.QE_materialdensity = 0.55 * 7.4 #Estimated packing density
        detobj.QE_matthickness = 0.100 #mm
        
        def DQEfun(**kwargs):
            return 0.05/0.8*XRDDetector.QEfunc_generic(**kwargs)
        
        detobj.QEfunc = XRDDetector.QEfunc_generic
        detobj.DQEfunc = DQEfun
        
        def maskfun(img,Mx,Ny):
            maxradius_inpixels = 150*1000/detobj.pixelsize/2
            img = maskfun_generic(img, Mx, Ny, maxradius_inpixels)
            return img
        
        detobj.maskfunction = maskfun
        detobj.name = 'DCS 4-Frame Detector (150mm)'
        
        
        return detobj
    
    def blurbygaussian_self(self,img):
        
        # [N,M] = np.shape(img)
        # x = (np.arange(N)-int(N/2))  # 1-D x-array in pixels
        # y = (np.arange(M)-int(M/2)) # 1-D y-array in pixels
        # xx, yy = np.meshgrid(x, y)  # make 2-D array of x's and y's
        # stdsize = self.InstrFWHM/self.pixelsize # std dev in pixels
        # cw = np.sqrt(-((stdsize / 2.) ** 2) / np.log(0.5))  # for gaussian blur with std erfsize
        # filt_blur = np.exp(-(xx ** 2 + yy ** 2) / cw ** 2)  # camera blur in real space
        # filt_blur = np.abs(fft2(filt_blur))
        # filt_blur /= filt_blur.max()  # normalize so max = 1
        # return np.abs(ifft2(fft2(img) * filt_blur))
        return XRDDetector.blurbygaussian(img, self.InstrFWHM, self.pixelsize)
    
    @staticmethod
    def QEfunc_generic(mat_thickness_inmm=None,mat_density=None,mat_formula=None,energyineV=None):
        if any([arg is None for arg in [mat_thickness_inmm,mat_density,mat_formula,energyineV]]):
            return None
        
        attenlength_inmm=get_attenlengthinmm_single(mat_formula,mat_density, energyineV)
        QE = 1 - np.exp(-mat_thickness_inmm/attenlength_inmm)
        return QE
    
    @staticmethod
    def blurbygaussian(img,fwhminmm,pixelsizeinmm):
        [N,M] = np.shape(img)
        x = (np.arange(N)-int(N/2))  # 1-D x-array in pixels
        y = (np.arange(M)-int(M/2)) # 1-D y-array in pixels
        xx, yy = np.meshgrid(x, y)  # make 2-D array of x's and y's
        stdsize = fwhminmm/pixelsizeinmm # std dev in pixels
        cw = np.sqrt(-((stdsize / 2.) ** 2) / np.log(0.5))  # for gaussian blur with std erfsize
        filt_blur = np.exp(-(xx ** 2 + yy ** 2) / cw ** 2)  # camera blur in real space
        filt_blur = np.abs(fft2(filt_blur))
        filt_blur /= filt_blur.max()  # normalize so max = 1
        # plt.imshow(abs(ifft2(filt_blur)))
        # plt.show()
        return np.abs(ifft2(fft2(img) * filt_blur))
    
    
def maskfun_generic(img,Mx,Ny,maxradius_inpixels):
    # maxradius_inpixels = np.min([Mx,Ny])/2
    xc = (Mx-1)/2 +1
    yc = (Ny-1)/2 +1
    imagecolnumbers_x = np.tile(np.arange(Mx).reshape(1, -1), (Ny, 1))
    imagerownumbers_y = np.tile(np.arange(Ny).reshape(-1, 1), (1, Mx))
    rmat =np.sqrt(np.square(imagecolnumbers_x - xc) + np.square(imagerownumbers_y - yc))
    img[rmat>maxradius_inpixels]=0.0
    return img  

def maskfun_circleinscribedinsquare(img,Mx,Ny):
    maxradius_inpixels = np.min([Mx,Ny])/2
    img = maskfun_generic(img,Mx,Ny,maxradius_inpixels)
    return img        

def convertspectrumRangeto2th(spec_EeV,spec_WpereV,dspacing):
    #Note all values are in radians, contrasting with the degree convention of most of this script
    # print(f'd: {dspacing}')
    
    hc = 12398.4198 # electron volts * angstrom    
    def convertEto2theta(Energies):
        argval = hc/2/Energies/dspacing
        argval[argval>1]=np.nan
        twoth = 2 * np.arcsin(argval) #in rad
        return twoth
        #hc is in electron volts * angstrom
        # dspacing is in Angstrom 
    
    def convertdPdEtodPdtheta(twothetas,PowerPereV):
        #Scale Power Per eV by width of a 2theta degree 
        dEd2th = hc/4/dspacing/np.sin(twothetas/2)/np.tan(twothetas/2)
        return PowerPereV*dEd2th
    
    # print(f' E sum: {np.trapezoid(spec_WpereV,spec_EeV)}')


    minE = hc/2/dspacing
    
    spec_WpereV = spec_WpereV[spec_EeV>minE]
    spec_EeV = spec_EeV[spec_EeV>minE]
    twothetas_spectrum = convertEto2theta(spec_EeV)
    Wper2theta_rad = convertdPdEtodPdtheta(twothetas_spectrum,spec_WpereV)
    NperSecper2theta_rad = Wper2theta_rad/spec_EeV
    twothetas_spectrum = np.flip(twothetas_spectrum)
    Wper2theta_rad = np.flip(Wper2theta_rad)
    NperSecper2theta_rad = np.flip(NperSecper2theta_rad)
    
    # print(f'sum: {np.trapezoid(Wper2theta_rad,twothetas_spectrum)}')
    
    
    
    EquallySpacedtwothvals_rad = np.linspace(np.nanmax([0, np.min(twothetas_spectrum)]), np.nanmin([np.pi, np.max(twothetas_spectrum)]),5000)
    Wper2th_rad_interp = interp1d(twothetas_spectrum, Wper2theta_rad, kind='cubic', fill_value=0.0, bounds_error=False)
    EquallySpaceWper2th_rad = Wper2th_rad_interp(EquallySpacedtwothvals_rad)
    NperSecPer2th_rad_interp = interp1d(twothetas_spectrum, NperSecper2theta_rad, kind='cubic', fill_value=0.0, bounds_error=False)
    # print(f'sum: {np.trapezoid(EquallySpaceWper2th_rad,EquallySpacedtwothvals_rad)}')
    
    
    
    return Wper2th_rad_interp, NperSecPer2th_rad_interp
        
    
class ExtraFigureWindow(QMainWindow):
    def __init__(self, parent_obj):
        super().__init__()

        self.parent_obj = parent_obj
        self.init_ui()
        

    def init_ui(self):
        
        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QVBoxLayout(self._main)

        self.static_canvas = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(NavigationToolbar_sub(self.static_canvas, self))
        layout.addWidget(self.static_canvas)


        self.ax = self.static_canvas.figure.subplots()
        # self.ax.set_aspect(self.parent_obj.aspectratio)
        self.static_canvas.draw()
        
def  AttenuationFactorAngledBeam2(attenlength,thickness,Alpharadians, TwothetaRadians,PhiRadians, gamma,beta,boolmask_infattenuation):
    ''''%This one is not the average attenuation (1/thickness)*Int(dz), it's the integral
    %of Int(dN/dz) so it doesn't have the 1/thickness 
    
    %output has units of Angstrom if attenlength and thickness are in microns
    
    %Alpha is angle of incident beam w.r.t material surface normal;
    %thickness measured in direction of normal (as standard), attenlength and
    %thickness in micrometers

    '''
    
    
    thickness = np.float64(thickness)
    # b = (np.exp(-gamma*thickness/attenlength) - np.exp(-beta*thickness/attenlength))
    beta[boolmask_infattenuation]=0 #to eliminate exp overflow error, fixed in later boolmask assignment
    ratioOuttoIn = np.zeros_like(beta)
    ratioOuttoIn = np.float64(10000.0)*attenlength*gamma/(beta - gamma)*(np.exp(-gamma*thickness/attenlength) - np.exp(-beta*thickness/attenlength))
    ratioOuttoIn[gamma==beta] = np.float64(10000.0) * gamma * thickness * np.exp(-gamma*thickness/attenlength)
    ratioOuttoIn[boolmask_infattenuation]=0
    # self.showimageinnewwindow(a)
    # self.showimageinnewwindow(b)
    # self.showimageinnewwindow(a*b)
    ratioOuttoIn[np.isnan(ratioOuttoIn)]=1e12
    # ratioOuttoIn[np.isinf(ratioOuttoIn)]=1e12
    # self.showimageinnewwindow(ratioOuttoIn)
    return ratioOuttoIn
        
#placeholder for eventually using multiprocessing
# def calcXRD_singleE(package):
#     '''
#             package={'gamma':gamma,
#                      'beta':beta,
#                      'boolmask':boolmask_infattenuation,
#                      'energy':peak_energy,
#                      'attenlengthinmm':attenlengthinmm,
#                      'd_spacing':d_spacing,
#                      'minE':minE,
#                      'maxE':maxE,
#                      'maxtwotheta':self.maxtwotheta,
#                      'x':x,
#                      'y':y,
#                      'twothetaimg':self.TwoTheta_rad_9x1[4],
#                      'phiimg':self.Phi_Rad_9x1[4],
#                      'thicknessinmm':thicknessinmm,
#                      'Alpharadians':Alpharadians,
#                      'structuredata':structuredata,
#                      'h':h,
#                      'k':k,
#                      'l':l,
#                      'f0_values':f0_values,
                     
#                      'multiplicity':multiplicity,
#                      'volume':volume,
#                      'bunchfreq':bunchfreq,
#                      'dtwotheta':dtwotheta,
#                      'dPhi':dPhi,
#                      'LPFactor':LPFactor
#                      }
            
      
#             '''
#     gamma = package['gamma']            
#     beta = package['beta']                    
#     boolmask_infattenuation = package['boolmask']                    
#     peak_energy=package['peak_energy']
#     attenlengthinmm = package['attenlengthinmm']                    
#     d_spacing = package['d_spacing']                    
#     minE = package['minE']                        
#     maxtwotheta = package['maxtwotheta']            
#     x = package['x']  
#     y = package['y']  
#     twothetaimg = package['twothetaimg']  
#     phiimg = package['phiimg']  
#     maxE = package['maxE']  
#     thicknessinmm = package['thicknessinmm']  
#     Alpharadians = package['Alpharadians']  
#     structuredata = package['structuredata']  
#     h = package['h']  
#     k = package['k']  
#     l = package['l']  
#     f0_values = package['f0_values']
    
#     multiplicity = package['multiplicity']
#     volume = package['volume']
#     bunchfreq = package['bunchfreq']
#     dtwotheta = package['dtwotheta']
#     dPhi = package['dPhi']
    
#     LPFactor = package['LPFactor']
#     dPhi = package['dPhi']
#     dPhi = package['dPhi']
    
    
    
    
#     wavelength = hc/peak_energy #approximation, can substitute a matrix with lam = 2*d*sin(2th/2) later
#     min_2th_rad = 2*np.arcsin(hc/2/d_spacing/maxE)
    
#     if min_2th_rad<maxtwotheta:
#         harmonic_E = x[(x>minE)&(x<maxE)] 
#         harmonic_dPdE = y[(x>minE)&(x<maxE)] 
#         effectivelength_Angstrom = AttenuationFactorAngledBeam2(attenlengthinmm*1000,thicknessinmm*1000,Alpharadians, twothetaimg,phiimg,gamma,beta,boolmask_infattenuation)
       
#         F_hkl = calculate_structure_factor(structuredata, h, k, l, peak_energy, f0_values)  
#         G = (r0**2)*(wavelength**3)*multiplicity*(F_hkl**2)/4/(volume**2) # [A^2]*[A^3]/[A^6] = 1/Angstrom
#         expectedatten = np.exp(-thicknessinmm/attenlengthinmm)

#         Wper2th_rad_interp, NperSecPer2th_rad_interp= convertspectrumRangeto2th(harmonic_E,harmonic_dPdE,d_spacing)
        
#         def calcN(th,phi):
#             return NperSecPer2th_rad_interp(th) / 2 / np.pi / bunchfreq * JouleIneV
        
        
        
#         NperSecper2th_rad = NperSecPer2th_rad_interp(twothetaimg)
#         Nper2thetaperphi_rad = NperSecper2th_rad / 2/np.pi / bunchfreq * JouleIneV
        
                
#         ScatteredImageSingle = Nper2thetaperphi_rad * G * LPFactor * effectivelength_Angstrom *dtwotheta * dPhi
        
#         print(f'E: {peak_energy:.3f}, [{h},{k},{l}], F: {F_hkl:.3f},Exp. Trans: {expectedatten:.3f}, ringsum: {np.sum(ScatteredImageSingle):.3e}')
#         # print(f'LP: {np.max(LPFactor):.3e}')
#         queue.put(ScatteredImageSingle)
#     else:
#         print(f'{peak_energy:.3f}: entire harmonic range off image')
#         # ScatteredImage.append(np.zeros_like(self.rmat))
#         queue.put(np.zeros_like(self.rmat))
        

def find_peaks_in_data(x, y, height=None, distance=None):
    peaks, _ = find_peaks(y, height=height, distance=distance)
    peak_locations = x[peaks]

    return len(peaks), peak_locations.tolist()
    

def calculate_multiplicity(structure, h, k, l):
    """Calculate the multiplicity of a given Miller index (h, k, l)."""
    # Extract symmetry information using spglib
    cell = (
        structure.cell,  # Lattice vectors
        structure.get_scaled_positions(),  # Atomic positions (fractional coordinates)
        structure.get_atomic_numbers(),  # Atomic numbers
    )
    symmetry_dataset = spglib.get_symmetry_dataset(cell)

    # Get symmetry operations
    rotations = symmetry_dataset.rotations
    # translations = symmetry_dataset.translations

    # Apply symmetry operations to the given Miller index
    miller_indices = []
    original_index = np.array([h, k, l])

    for rotation in rotations:
        # Transform the Miller index
        transformed_index = np.dot(rotation, original_index)
        miller_indices.append(tuple(transformed_index))

    # Deduplicate equivalent reflections
    unique_indices = set(miller_indices)

    # Return the number of unique reflections (multiplicity)
    return len(unique_indices)


# def generate_unique_hkl_withmult(structure, max_index):
#     """Generate unique (h, k, l) indices up to max_index and their multiplicities."""
#     # Extract symmetry operations from the structure
#     cell = (
#         structure.cell,  # Lattice vectors
#         structure.get_scaled_positions(),  # Atomic positions (fractional coordinates)
#         structure.get_atomic_numbers(),  # Atomic numbers
#     )
#     symmetry_dataset = spglib.get_symmetry_dataset(cell)
#     rotations = symmetry_dataset.rotations

#     # Generate all possible (h, k, l) combinations
#     hkl_list = [(h, k, l) for h in range(-max_index, max_index + 1)
#                           for k in range(-max_index, max_index + 1)
#                           for l in range(-max_index, max_index + 1) if h != 0 or k != 0 or l != 0]

#     # Dictionary to store unique (h, k, l) and their multiplicities
#     unique_hkl = {}

#     for hkl in hkl_list:
#         # Apply all symmetry operations to this hkl
#         transformed_hkls = {tuple(np.dot(rotation, hkl)) for rotation in rotations}

#         # Normalize to ensure consistent representation
#         normalized_hkl = tuple(sorted(transformed_hkls, key=lambda x: (x[0], x[1], x[2]))[0])

#         # Count the multiplicity
#         if normalized_hkl in unique_hkl:
#             unique_hkl[normalized_hkl] += len(transformed_hkls)
#         else:
#             unique_hkl[normalized_hkl] = len(transformed_hkls)

#     # Return unique (h, k, l) as a sorted list with multiplicities
#     unique_hkl_list = sorted(unique_hkl.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))
#     return unique_hkl_list




# def generate_unique_hkl_ffff(structure, max_index):
#     """Generate unique (h, k, l) indices up to max_index and their multiplicities."""
#     # Extract symmetry operations from the structure
#     cell = (
#         structure.cell,  # Lattice vectors
#         structure.get_scaled_positions(),  # Atomic positions (fractional coordinates)
#         structure.get_atomic_numbers(),  # Atomic numbers
#     )
#     symmetry_dataset = spglib.get_symmetry_dataset(cell)
#     rotations = symmetry_dataset.rotations
#     print(rotations)

#     # Generate all possible (h, k, l) combinations
#     hkl_list = [(h, k, l) for h in range(-max_index, max_index + 1)
#                           for k in range(-max_index, max_index + 1)
#                           for l in range(-max_index, max_index + 1) if h != 0 or k != 0 or l != 0]

#     # Dictionary to store unique (h, k, l) and their multiplicities
#     unique_hkl = {}

#     for hkl in hkl_list:
#         # Apply all symmetry operations to this hkl
#         transformed_hkls = {tuple(np.dot(rotation, hkl)) for rotation in rotations}

#         # Normalize and canonicalize
#         normalized_hkls = {canonicalize_hkl(h) for h in transformed_hkls}

#         # Add to unique_hkl
#         for h in normalized_hkls:
#             if h in unique_hkl:
#                 unique_hkl[h] += 1
#             else:
#                 unique_hkl[h] = 1

#     # Return unique (h, k, l) as a sorted list with multiplicities
#     unique_hkl_list = sorted(unique_hkl.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))
#     return unique_hkl_list

# def remove_equivalent(hkllist,)

def generate_unique_hkl(structure, max_index):
    """Generate unique (h, k, l) indices up to max_index and their multiplicities."""
    # Extract symmetry operations from the structure
    cell = (
        structure.cell,  # Lattice vectors
        structure.get_scaled_positions(),  # Atomic positions (fractional coordinates)
        structure.get_atomic_numbers(),  # Atomic numbers
    )
    symmetry_dataset = spglib.get_symmetry_dataset(cell)
    rotations = symmetry_dataset.rotations
    identityrotation = np.array([[1,0,0],[0,1,0],[0,0,1]],type(rotations[0][0][0]))
    rotations = [r for r in rotations if not np.array_equal(r,identityrotation)]
    rotations.append(identityrotation) #only one identity rotation
    # print(rotations)

    # Generate all possible (h, k, l) combinations
    hkl_list = [(h, k, l) for h in range(-max_index, max_index + 1)
                          for k in range(-max_index, max_index + 1)
                          for l in range(-max_index, max_index + 1) if h != 0 or k != 0 or l != 0]
    
    # Dictionary to store unique (h, k, l) and their multiplicities
    unique_hkl_list =[]

    

    for ij in range(len(hkl_list)):
        hkl = hkl_list[ij]
        if not hkl == (None,None,None):
            # Apply all symmetry operations to this hkl
            transformed_hkls = {tuple(np.dot(rotation, hkl)) for rotation in rotations}
            #replace all equivalent reflections in the list with None (to skip)
            matches = transformed_hkls.intersection(hkl_list)  # Find matches
            hkl_list =[(None,None,None) if t in transformed_hkls else t for t in hkl_list]
            hkl_list[ij] = hkl #previous condition also removes the initial refl.

            multiplicity = len(matches) # Count matches
            refl = {'hkl': hkl,'multiplicity': multiplicity}
            unique_hkl_list.append(refl)
        
        
            # Normalize and canonicalize
            # normalized_hkls = {canonicalize_hkl(h) for h in transformed_hkls}
    


    # Return unique (h, k, l) as a sorted list with multiplicities
    # unique_hkl_list = sorted(unique_hkl.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))
    return unique_hkl_list


# def canonicalize_hkl(hkl):
#     """Canonicalize (h, k, l) by sorting absolute values while preserving the sign of the largest index."""
#     abs_hkl = sorted(hkl, key=lambda x: (-abs(x), x))  # Sort by absolute value, preserve sign of largest
#     return tuple(abs_hkl)



def compute_structure_factors(structure, theta_max, energy, progress_callback=None):
    wavelength = hc/energy
    dmax = max(main_window.structure_data.get_cell().lengths())
    """Compute structure factors for all (h, k, l) up to given two-theta at a specific energy."""
    max_index = int(2 * dmax / wavelength * np.sin(theta_max*np.pi/180/2))
                    
    
    unique_hkl_list = generate_unique_hkl(structure, max_index)  # Get unique (h, k, l)
    
    
    hklvals = []
    angles = []
    structure_factors = []
    multiplicities = []
    total_steps = len(unique_hkl_list)
    elements = set(atom.symbol for atom in structure)  # Extract unique elements
    f0_values = precompute_f0(xraylib, elements, x_min=0.0, x_max=8.0, num_points=100)

    refl = 0
    for i, elem in enumerate(unique_hkl_list):
        h,k,l = elem['hkl']
        multiplicity=elem['multiplicity']
        d_spacing = compute_d_spacing(structure, h, k, l)
        if d_spacing > 0 and d_spacing > wavelength / 2:  # Ensure d_spacing is valid
            two_theta = 2 * np.degrees(np.arcsin(wavelength / (2 * d_spacing)))
            if two_theta <= theta_max:
                F_hkl = calculate_structure_factor(structure, h, k, l, energy, f0_values)  
                angles.append(two_theta)
                structure_factors.append(np.abs(F_hkl))
                multiplicities.append(multiplicity)
                hklvals.append(elem['hkl'])

        # Update progress
        if progress_callback:
            progress_callback(i + 1, total_steps)


    return hklvals,angles, structure_factors,multiplicities,elements,f0_values

def calcLPFactor(two_theta_angles_rad,PolarizedFlag=0,phi_angle_rad=0):
        '''polarization factor for power emitted from entire crystal (not per unit length on diffraction circle)
        
        See Warren X-ray Diffraction Equation 1.2 and discussion immediately prior
        E^2 = Ez^2  + Ey^2*cos(phi)^2. where phi is scattering angle  in that discussion
        
        '''
        if PolarizedFlag ==0:
            LPFactor =0.5*(1 + np.square(np.cos(two_theta_angles_rad)))/np.sin(two_theta_angles_rad/2)
        elif PolarizedFlag ==1:
            LPFactor =(np.square(np.sin(phi_angle_rad)) + np.square(np.cos(phi_angle_rad))*np.square(np.cos(two_theta_angles_rad)))/np.sin(two_theta_angles_rad/2)/np.sin(two_theta_angles_rad)
        
        return LPFactor


def _calcXRDsignalstrength(hklvals, two_theta_angles,F_HKL,multiplicities,LPFlag=0):
        '''LP Flag, if calculating total energy in a reflection, LPFlag = 0
        if calculating itnensity per unit length on a diffraction circle   LPFlag = 1

        
        # LPFactor =0.5*(1 + np.square(np.cos(two_theta_angles*np.pi/180)))/np.sin(two_theta_angles*np.pi/180/2)
        '''
        if LPFlag ==0:
            LPFactor =(1 + np.square(np.cos(two_theta_angles*np.pi/180)))/np.sin(two_theta_angles*np.pi/180/2)
        elif LPFlag ==1:
            LPFactor =(1 + np.square(np.cos(two_theta_angles*np.pi/180)))/np.sin(two_theta_angles*np.pi/180/2)/np.sin(two_theta_angles*np.pi/180)
        
            
        for hkl, th,F,LP,m in zip(hklvals, two_theta_angles,F_HKL,LPFactor,multiplicities):
            if F>1:
                print(f'hkl: {hkl}, 2th:{th:.3f},m: {m: 3d}, F**2:{F**2:6.0f},LP:{LP: 6.3f}')
        return np.square(F_HKL) * LPFactor * multiplicities
    
def gaussian(x, mu, sig):
    return (
        2.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
    )    

def FWHMFromInstrumentAndCrystalSize(instrfwhmindegrees,twotheta,energyineV,crystalthicknessinmicrons):
    #
    """ From Scherrer Formula and instr broadening. 
    For further development: “Scherrer after sixty years: A survey and some new
    results in the determination of crystallite size,” J. Appl. Cryst. 11
    (1978) 102-113. """
    K = 0.9 #Shape Factor
    wavelengthinmicrons = 1.23984198/energyineV
    B = K*wavelengthinmicrons/crystalthicknessinmicrons/np.cos(twotheta/2*np.pi/180)
    return np.sqrt(B**2 + (instrfwhmindegrees*np.pi/180)**2)

def search_cod_fortabledata(chemical_formula):
    # Create a session object
    with requests.Session() as session:
        # URL for the initial search POST
        search_url = "https://www.crystallography.net/cod/result.php"

        # Form data for your search
        form_data = {
            'formula': chemical_formula,  # example chemical formula; replace as needed
        }

        # Headers (as used previously)
        headers = {
            'dnt': '1',
            'referer': 'https://www.crystallography.net/cod/search.html',
            'sec-ch-ua': '"Google Chrome";v="131", "Chromium";v="131", "Not_A Brand";v="24"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"Windows"',
            'upgrade-insecure-requests': '1',
            'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
        }

        response = session.post(search_url, data=form_data, headers=headers)
        if response.status_code != 200:
            print(f"Failed to fetch search results. Status code: {response.status_code}")
            return []

        soup = BeautifulSoup(response.text, 'html.parser')

        #  Find the pagination link with the maximum count
        linkresults = []
        counts = []
        for link in soup.find_all('a', href=True):
            href = link['href']
            if 'count=' in href and '&page=0' in href:
                start = href.find('count=') + len('count=')
                end = href.find('&page=0')
                count_value = int(href[start:end])
                counts.append(count_value)
                linkresults.append(href)

        if not counts:
            print("No pagination links found.")
            soup_full = soup #assume we've got all listed on one page
        else:
            
            max_link_index = np.argmax(counts)
            relative_url = linkresults[max_link_index]
            strtorep = 'count='+str(counts[max_link_index])
            relative_url_new = relative_url.replace(strtorep,'count=1000')
        
            # Combine with the base URL
            base_url = "https://www.crystallography.net/cod/"
            full_url = urllib.parse.urljoin(base_url, relative_url_new)
            print("Following URL:", full_url)
        
            full_results_response = session.get(full_url, headers=headers)
            if full_results_response.status_code != 200:
                print(f"Failed to fetch full results. Status code: {full_results_response.status_code}")
                exit()
        
            soup_full = BeautifulSoup(full_results_response.text, 'html.parser')
        
        # Extract table data from the full results page
        table = soup_full.find('table', attrs={'class': 'information'})
        if table:
            data_head = []
            data = []
            for row in table.find_all('tr'):
                headers_row = [th.get_text(strip=True) for th in row.find_all('th')]
                if headers_row:
                    data_head.append(headers_row)
            for row in table.find_all('tr'):
                cells = [td.get_text(strip=True) for td in row.find_all('td')]
                if cells:
                    data.append(cells)
            # print("Table headers:", data_head)
            # print("Table data:", data)
            return data, data_head
        else:
            print("Table not found in the full results page.")
            return []
    
def readtabledata(table_row_data,table_headers):
    
    for row in table_headers:
        if not row==[]:
            if 'ormula' in ''.join(row): #if 'formula' somewhere in the header row
                headerdata = row #then this is the header data
    
    def findstrandputindict(header,data,stringtomatch,newdict,label):
        if stringtomatch in header.casefold():
            newdict[label]=data
        return newdict

    
    dictlist=[]
    for row in table_row_data:
        if not row == []:
            newdict={}
            for coldata,colheader in zip(row,headerdata):
                # print(coldata,colheader)
                newdict = findstrandputindict(colheader,coldata,'id',newdict,'COD ID')
                newdict = findstrandputindict(colheader,coldata,'formula',newdict,'Chemical Formula')
                newdict = findstrandputindict(colheader,coldata,'param',newdict,'Cell Parameters')
                newdict = findstrandputindict(colheader,coldata,'volume',newdict,'Cell Vol.')
                newdict = findstrandputindict(colheader,coldata,'group',newdict,'Space Group')
                newdict = findstrandputindict(colheader,coldata,'biblio',newdict,'Bibliography')
            dictlist.append(newdict)
    return dictlist

def calculate_density(atoms):
    """
    Calculate the density of a structure from an ASE Atoms object.

    Returns:
        float: Density in g/cm³.
    """

    volume = atoms.get_volume()  # unit cell volume in Å³
    total_mass = sum(atoms.get_masses())  # Mass in amu
    total_mass_grams = total_mass * 1.66054e-24
    volume_cm3 = volume * 1e-24 #(1 Å³ = 1e-24 cm³)
    density = total_mass_grams / volume_cm3 #density in g/cm³
    return density


def get_attenlengthinmm_single(chemformula,density, energyineV):
    cs = xraylib.CS_Total_CP(chemformula,energyineV/1000)
    return 10/(density*cs)

def get_attenlengthinmm(chemformula,density, energies):
    
    cs_cm2perg = np.zeros(len(energies))
    for j,energy in enumerate(energies):
        cs_cm2perg[j] = xraylib.CS_Total_CP(chemformula,energy/1000.0)    
    return 10/(density*cs_cm2perg)



def compute_d_spacing(structure, h, k, l):
    """Compute d-spacing for a given (h, k, l)."""
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 

    # Metric tensor for the lattice
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    # V = (a * b * c) * np.sqrt(
    #     1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 +
    #     2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)
    # )
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    hkl = np.array([h, k, l])

    try:
        d_spacing = 1 / np.sqrt(hkl @ np.linalg.inv(G) @ hkl.T)
    except np.linalg.LinAlgError:
        return 0  # Return 0 if the metric tensor is singular (rare edge case)

    return d_spacing




def calculate_structure_factor(structure, h, k, l, energy, f0_values):
    """Calculate the structure factor for the given (h, k, l) at a specific energy."""
    F_hkl = 0.0 + 0.0j  # Initialize as complex

    # Calculate sin(theta)/lambda (x) for this (h, k, l)
    d_spacing = compute_d_spacing(structure, h, k, l)
    if d_spacing == 0:
        return 0  # Skip invalid reflections
    x = 1 / (2 * d_spacing)  # sin(theta)/lambda
    # print(f'x = {x}')

    for atom in structure:
        # fractional_coords = atom.position   # Convert to fractional coordinates
        x_pos, y_pos, z_pos = atom.position / structure.get_cell().lengths()
        element = atom.symbol
        # print(f'E: {energy}, hkl = {(h,k,l)}, x,y,z: {x_pos},{y_pos},{z_pos}: {element}')
        # Interpolate f0 for this x value
        f0 = f0_values[element](x)  # Interpolated f0
        atomic_number = atom.number  # Atomic number Z
        
        
        
        f1,f2 = f1f2_calc(element, energy, theta=None, F=0, density=None, rough=None, verbose=False,
                      material_constants_library=xraylib)
        f1_corrected = f1 - atomic_number
        # print(f'{element}: f0 = {f0}, f1 = {f1_corrected}, , f2 = {f2}')


        # Combined scattering factor
        f = f0 + f1_corrected + 1j * f2

        

        # Phase contribution
        phase = np.exp(2j * np.pi * (h * x_pos + k * y_pos + l * z_pos))

        # print(f'(x,y,z): ({x_pos},{y_pos},{z_pos}):f with phase: f = {f*phase}')
        # Accumulate structure factor
        F_hkl += f * phase
        # print(f'{element}: F={F_hkl}')
    

    return np.float64(abs(np.squeeze(F_hkl)))

def calculate_structure_factor_withDebyeWaller(structure, h, k, l, energy, f0_values, Bdict):
    """Calculate the structure factor for the given (h, k, l) at a specific energy."""
    F_hkl = 0.0 + 0.0j  # Initialize as complex

    # Calculate sin(theta)/lambda (x) for this (h, k, l)
    d_spacing = compute_d_spacing(structure, h, k, l)
    if d_spacing == 0:
        return 0  # Skip invalid reflections
    x = 1 / (2 * d_spacing)  # sin(theta)/lambda
    # print(f'x = {x}')
    

    for atom in structure:
        element = atom.symbol
        for atomtype, Bval in zip(Bdict['el_list'],Bdict['B_list']):
            if element==atomtype:
                B = Bval
                
        # fractional_coords = atom.position   # Convert to fractional coordinates
        x_pos, y_pos, z_pos = atom.position / structure.get_cell().lengths()
        
        # print(f'E: {energy}, hkl = {(h,k,l)}, x,y,z: {x_pos},{y_pos},{z_pos}: {element}')
        # Interpolate f0 for this x value
        f0 = f0_values[element](x)  # Interpolated f0
        atomic_number = atom.number  # Atomic number Z
        
        
        
        f1,f2 = f1f2_calc(element, energy, theta=None, F=0, density=None, rough=None, verbose=False,
                      material_constants_library=xraylib)
        f1_corrected = f1 - atomic_number
        # print(f'{element}: f0 = {f0}, f1 = {f1_corrected}, , f2 = {f2}')


        # Combined scattering factor
        f = f0 + f1_corrected + 1j * f2

        

        # Phase contribution
        phase = np.exp(2j * np.pi * (h * x_pos + k * y_pos + l * z_pos))

        # print(f'(x,y,z): ({x_pos},{y_pos},{z_pos}):f with phase: f = {f*phase}')
        # Accumulate structure factor
        F_hkl += f * phase * np.exp(-B*(x**2))
        # print(f'{element}: F={F_hkl}')
    

    return np.float64(abs(np.squeeze(F_hkl)))

def sum_y_for_same_x(x, y):
    """
    Sum all y values that correspond to the same x value.

    Parameters:
        x (list or array): List of x values.
        y (list or array): List of y values corresponding to the x values.

    Returns:
        result_x: List of unique x values.
        result_y: List of summed y values for each unique x.
    """
    xy_sums = defaultdict(float)

    for xi, yi in zip(x, y):
        xy_sums[xi] += yi

    # Extract results as two lists
    result_x = list(xy_sums.keys())
    result_y = list(xy_sums.values())

    return result_x, result_y


def sum_y_within_maxsep(x, y, maxsep, maxgroupspan):
    
    '''
    group x,y data for fluorescence energy, Nphotons, so we don't have to
    calculated attenuation map for each line
    '''
    totenergy = np.array(x)*np.array(y)
    data = pd.DataFrame({'x': x, 'y': y,'totenergy': totenergy}).sort_values('x').reset_index(drop=True)

    # Initialize grouping
    groups = [0]  # First item is always in group 0
    current_group = 0
    groupstartx = x[0]
    
    # Assign groups based on maxsep
    for i in range(1, len(data)):
        if (data['x'][i] - data['x'][i - 1] > maxsep) or (data['x'][i] - groupstartx > maxgroupspan) :
            current_group += 1
            groupstartx = data['x'][i]
        groups.append(current_group)
    
    data['group'] = groups
    
    
    # Aggregate by group
    result_df = data.groupby('group', as_index=False).agg({'x': 'mean', 'y': 'sum','totenergy': 'sum'})
    weightedmeanenergy = result_df['totenergy']/result_df['y']
    # return result_df['x'],result_df['y']
    return np.array(weightedmeanenergy), np.array(result_df['y'])


def AttenuationAngleFactors(Alpharadians, TwothetaRadians,PhiRadians):
    
    gamma = 1/np.cos(Alpharadians)
    denom = 1 - (np.tan(Alpharadians)*np.tan(TwothetaRadians)*np.cos(PhiRadians));
    beta = (gamma/np.cos(TwothetaRadians))/denom
    # ycenterpix = int(self.ycenterpix_input.text())
    # xcenterpix = int(self.xcenterpix_input.text())
    beta[beta==gamma]=gamma+1e-7
    boolmask_infattenuation = denom<1e-12
    
    # beta[ycenterpix,xcenterpix]=gamma+1e-7
    # self.showimageinnewwindow(1-beta)
    return gamma,beta,boolmask_infattenuation

def  AttenuationFactorAngledBeam2_Filter(attenlength,thickness, gamma,beta,boolmask_infattenuation):
   
    beta[boolmask_infattenuation]=0 #to eliminate exp overflow error, fixed in later boolmask assignment
    ratioOuttoIn = np.exp(-beta*thickness/attenlength)
    ratioOuttoIn[boolmask_infattenuation]=0
    ratioOuttoIn[np.isnan(ratioOuttoIn)]=0
    # self.showimageinnewwindow(ratioOuttoIn)
    return ratioOuttoIn

def to_hill_notation(formula):
    # Regular expression to extract elements and counts
    element_pattern = r"([A-Z][a-z]*)(\d*)"
    
    # Parse the formula into elements and counts
    elements = re.findall(element_pattern, formula)
    element_counts = Counter({elem: int(count) if count else 1 for elem, count in elements})
    
    # Find the greatest common divisor (GCD) of all counts
    counts = list(element_counts.values())
    common_divisor = reduce(gcd, counts) if counts else 1
    
    # Reduce all counts by the GCD
    for element in element_counts:
        element_counts[element] //= common_divisor
    
    # Separate carbon, hydrogen, and others
    carbon_count = element_counts.pop('C', 0)
    hydrogen_count = element_counts.pop('H', 0)
    others = sorted(element_counts.items())
    
    # Format the output string
    result = []
    if carbon_count > 0:
        if carbon_count==1:
            carbon_count=''
        result.append(f"C{carbon_count}")
    if hydrogen_count > 0:
        if hydrogen_count==1:
            hydrogen_count=''
        result.append(f"H{hydrogen_count}")
    
    for elem,count in others:
        if count==1:
            count=''
        result.append(f"{elem}{count}")
    
    return ' '.join(result)


def  AttenuationFactorAngledBeam2_Fluorescence(attenlength0,attenlength_fluor,thickness,Alpharadians, TwothetaRadians,PhiRadians,gamma,beta,boolmask_infattenuation):
    ''' Same as AttenuationFactorAngledBeam2,except two different atten
    lengths, 1 for input beam and 1 for fluorescence

    '''
        
    thickness = np.float64(thickness)
    # b = (np.exp(-gamma*thickness/attenlength) - np.exp(-beta*thickness/attenlength))
    beta[boolmask_infattenuation]=0 #to eliminate exp overflow error, fixed in later boolmask assignment
    ratioOuttoIn = np.zeros_like(beta)
    ratioOuttoIn = np.float64(10000.0)*gamma/(beta/attenlength_fluor - gamma/attenlength0)*(np.exp(-gamma*thickness/attenlength0) - np.exp(-beta*thickness/attenlength_fluor))
    ratioOuttoIn[gamma==beta] = np.float64(10000.0) * gamma * thickness * np.exp(-gamma*thickness/attenlength0)
    ratioOuttoIn[boolmask_infattenuation]=0
    # self.showimageinnewwindow(a)
    # self.showimageinnewwindow(b)
    # self.showimageinnewwindow(a*b)
    ratioOuttoIn[np.isnan(ratioOuttoIn)]=1e12
    # ratioOuttoIn[np.isinf(ratioOuttoIn)]=1e12
    # self.showimageinnewwindow(ratioOuttoIn)
    return ratioOuttoIn

def read_API_key():

    try:
        with open(API_KEY_LOCATION, 'r') as file:
            # Read the first (and only) line of the file
            content = file.readline().strip()
            # Ensure the content is alphanumeric
            if not content.isalnum():
                raise ValueError("The file contains non-alphanumeric characters.")
            return content
    except FileNotFoundError:
        print(f"Error: File not found at {API_KEY_LOCATION}")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")



def precompute_f0(material_constants_library, elements, x_min, x_max, num_points):
    """Precompute f0 values for elements using xoppy_calc_f0."""
    f0_values = {}

    for element in elements:
        # Call xoppy_calc_f0 to compute f0 values on a grid
        f0_data = xoppy_calc_f0(
            descriptor=element,
            MAT_FLAG=0,
            GRIDSTART=x_min,
            GRIDEND=x_max,
            GRIDN=num_points,
            DUMP_TO_FILE=0,
            FILE_NAME="f0.dat",
            CHARGE=0.0,
            material_constants_library=xraylib,
        )


        # Extract x (grid) and f0 values
        x_grid = f0_data['data'][0,:]  # sin(theta)/lambda
        f0_grid = f0_data['data'][1,:]  # f0 values

        # Interpolate f0 for this element
        f0_values[element] = interp1d(x_grid, f0_grid, kind='cubic', fill_value="extrapolate")

    return f0_values

def debye_1(x):
    """
    Computes the first-order Debye function D1(x).

    :param x: Input value (dimensionless temperature).
    :return: First-order Debye function D1(x).
    """
    if x == 0:
        return 0.0  # Handle the edge case for x = 0

    # Define the integrand
    def integrand(t):
        return t / (np.exp(t) - 1)

    # Perform the integral from 0 to x
    result, _ = quad(integrand, 0, x)

    # Compute D1(x)
    return result / x

def interpolate_debye_waller(element: str, temperature: float):

    # Load the data
    data = pd.read_csv(Debye_Waller_elementalcrystalpath,header=0, sep=',')
    # return data

    # Ensure the element exists in the data
    if element not in data.columns:
        print('Not in elemental crystal Debye-Waller Data. Get DW from elsewhere.')
        return None

    # Extract temperature and element values, skipping hte first element which is text
    
    temp_values = np.array(data['T'].values[1:],dtype='float64')
    element_values = np.array(data[element].values[1:],dtype='float64')

    # Interpolate the value
    interpolated_value = np.interp(temperature, temp_values, element_values)

    return interpolated_value

def chemformula(atoms_list_string):
    compound = xraylib.CompoundParser(atoms_list_string)
    chemformlist=[]
    gcdenom = math.gcd(*np.int8(compound['nAtoms']))
    for atomicnum,n in zip(compound['Elements'],compound['nAtoms']):
        n = n/gcdenom
        if n>1:
            chemformlist.append(f'{xraylib.AtomicNumberToSymbol(atomicnum)}{int(n):d}')
        else:
            chemformlist.append(f'{xraylib.AtomicNumberToSymbol(atomicnum)}')
    return ''.join(chemformlist)

def calcDWfromDebyeTemp(tempinK, DtempinK,atomicmass_gpermol):
    x=tempinK/DtempinK
    B = 6*(6.62e-27)**2*(tempinK) / ((atomicmass_gpermol)/(6.02e23)) * 1e16 / (1.38e-16) / (DtempinK)**2 * (debye_1(x)+x/4) #in Angstrom^2
    return B

class Atom3DViewer(QMainWindow):
    def __init__(self, atoms, cell, parent=None):
        """
        Initialize a 3D viewer for atom positions in the unit cell.

        Args:
            atoms (list): List of tuples [(x, y, z, symbol), ...] for atom positions and types.
            cell (np.ndarray): 3x3 array representing the unit cell vectors.
        """
        super().__init__(parent)
        self.setWindowTitle("3D Atom Viewer")
        self.setGeometry(100, 100, 800, 600)

        # Matplotlib Figure and Canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)

        # Central widget
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        layout.addWidget(self.canvas)
        self.setCentralWidget(central_widget)

        # Plot atoms
        self.plot_atoms(atoms, cell)

    def plot_atoms(self, atoms, cell):
        """Plot atoms in 3D along with the unit cell."""
        ax = self.figure.add_subplot(111, projection='3d')
        ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

        # Unit cell edges
        self.plot_unit_cell(ax, cell)

        # Atom positions and types
        for x, y, z, symbol in atoms:
            ax.scatter(x, y, z, s=300, label=symbol, alpha=0.8)  # Atom as sphere

        # Adjust plot
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.legend(loc='upper left', fontsize='small')
        self.canvas.draw()

    @staticmethod
    def plot_unit_cell(ax, cell):
        """Plot the edges of the unit cell."""
        # Create corners of the unit cell
        corners = np.array([
            [0, 0, 0],
            cell[0],
            cell[1],
            cell[2],
            cell[0] + cell[1],
            cell[0] + cell[2],
            cell[1] + cell[2],
            cell[0] + cell[1] + cell[2]
        ])

        # Define edges of the unit cell
        edges = [
            (0, 1), (0, 2), (0, 3),  # Edges from origin
            (1, 4), (1, 5),  # Edges from cell[0]
            (2, 4), (2, 6),  # Edges from cell[1]
            (3, 5), (3, 6),  # Edges from cell[2]
            (4, 7), (5, 7), (6, 7)  # Edges to opposite corner
        ]

        # Plot edges
        for edge in edges:
            ax.plot3D(
                [corners[edge[0], 0], corners[edge[1], 0]],
                [corners[edge[0], 1], corners[edge[1], 1]],
                [corners[edge[0], 2], corners[edge[1], 2]],
                color="black", linewidth=1
            )

class NavigationToolbar_sub(NavigationToolbar2QT):
    # only display the buttons we need
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] not in ('Back', 'Forward')]
    
    
def setFusionpalette(app):
    app.setStyle("Fusion")
    # Now use a palette to switch to dark colors:
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(53, 53, 53))
    palette.setColor(QPalette.WindowText, Qt.white)
    palette.setColor(QPalette.Base, QColor(25, 25, 25))
    palette.setColor(QPalette.AlternateBase, QColor(53, 53, 53))
    palette.setColor(QPalette.ToolTipBase, Qt.black)
    palette.setColor(QPalette.ToolTipText, Qt.white)
    palette.setColor(QPalette.Text, Qt.white)
    palette.setColor(QPalette.Button, QColor(53, 53, 53))
    palette.setColor(QPalette.ButtonText, Qt.white)
    palette.setColor(QPalette.BrightText, Qt.red)
    palette.setColor(QPalette.Link, QColor(42, 130, 218))
    palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    palette.setColor(QPalette.HighlightedText, Qt.black)
    app.setPalette(palette)

def blink_red(widgettoblink):
    # Set the edit box to red
    widgettoblink.setStyleSheet("QLineEdit { background-color: red; }")

    # Create a timer to revert the color back
    QTimer.singleShot(500, lambda: revert_color(widgettoblink))  # Change color back after 500 ms

def revert_color(widgettoblink):
    # Revert to the original style
    widgettoblink.setStyleSheet("")
    
# def search_cod(chemical_formula):
    
#     url = "https://www.crystallography.net/cod/result.php"
    
#     form_data = {
#         'formula': chemical_formula,  # Input for the chemical formula
#     }
    
#     # Define the headers (to mimic the browser's request)
#     headers = {
#         'dnt': '1',
#         'referer': 'https://www.crystallography.net/cod/search.html',  # Referer URL as seen in the browser
#         'sec-ch-ua': '"Google Chrome";v="131", "Chromium";v="131", "Not_A Brand";v="24"',
#         'sec-ch-ua-mobile': '?0',
#         'sec-ch-ua-platform': '"Windows"',
#         'upgrade-insecure-requests': '1',
#         'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
#     }
    
#     # Send the POST request
#     response = requests.post(url, data=form_data, headers=headers)
    
#     # Check if the request was successful
#     if response.status_code == 200:
#         # Parse the HTML response using BeautifulSoup
#         soup = BeautifulSoup(response.text, 'html.parser')
        
#         # Extract CIF file links
#         results = []
#         for link in soup.find_all('a', href=True):
#             if link['href'].endswith('.cif'):  # Look for CIF file links
#                 results.append(link['href'])
        
#         return results
#     else:
#         print(f"Failed to fetch data. Status code: {response.status_code}")
#         return []
    
# def search_cod(chemical_formula):
    
    
    
    
def getCODCif(IDnumberString,chemformula):
    url = f'http://www.crystallography.net/cod/{IDnumberString}.cif'
    output_directory = os.path.join(dir_path,'ciffiles')
    outfile = os.path.join(output_directory,chemformula+'_'+IDnumberString+'.cif')
    filename = wget.download(url, out=outfile)
    
    if os.path.exists(outfile):
        # print(f'File Written: {outfile}')
        return outfile
    else:
        return None
    

if __name__ == "__main__":
    app = QApplication(sys.argv)
    setFusionpalette(app)
    
    ciffile = os.path.join(dir_path,'ciffiles','Ta.cif')
    if not os.path.exists(ciffile):
        ciffile=None
    main_window = XRayScatteringApp(ciffile)
    
    main_window.show()
    
    
    sys.exit(app.exec_())
