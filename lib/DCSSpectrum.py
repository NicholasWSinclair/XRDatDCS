import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,QComboBox,
    QCheckBox, QLineEdit, QLabel, QMainWindow, QFormLayout, QAction, QMenuBar, 
    QFileDialog, QMessageBox, QGroupBox, QProgressBar, QSlider,QSizePolicy, 
    QSystemTrayIcon, QMenu,qApp, QDialog, QDialogButtonBox
)
from PyQt5.QtGui import QIcon
import itertools
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QThreadPool, QRunnable, pyqtSignal, QObject


import os
import h5py
from xoppylib.power.xoppy_calc_power import xoppy_calc_power
from xoppylib.sources.srundplug import calc1d_srw
try:
    import oasys_srw.srwlib as srwlib
except:
    from srwpy import srwlib
import xraylib
import scipy.constants as codata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.figure import Figure
from collections import OrderedDict

import csv
import scipy.interpolate
import time

try: # if using numpy < 2.0
    from numpy import trapz as trapezoid
except ImportError: # if using numpy >= 2.0, but still monkey patching for xoppylib which uses trapz
    from numpy import trapezoid
    np.trapz = trapezoid
    
# Local Imports -- These are in the same folder as this file
import FilterGUI
import MultilayerConfig
import MirrorGUI
import ConfigImportDialog
    
class DCSSpectrum(QWidget):
    def __init__(self):
        super().__init__()
        # print(sys.path)
        self.setWindowTitle("DCS Spectrum Calc GUI")
        self.layout = QVBoxLayout()
        self.initUI()

         
        

    def initUI(self):

        #File menu
        menu_bar = QMenuBar(self)
        file_menu = menu_bar.addMenu('File')
        Util_menu = menu_bar.addMenu('Utilities')
        Parameters_menu = menu_bar.addMenu('Parameters')
        
        
        
        SlitParams_action = QAction('Set White Beam Slits', self)        
        SlitParams_action.triggered.connect(self.open_gui_to_adjust_WBS)
        Parameters_menu.addAction(SlitParams_action)
        

        
        save_action = QAction('Save Spectrum To CSV', self)        
        save_action.triggered.connect(self.save_spectrum_toCSV)
        file_menu.addAction(save_action)

        save_config_action = QAction('Save Configuration', self)
        save_config_action.triggered.connect(self.save_config)
        file_menu.addAction(save_config_action)

        load_config_action = QAction('Load Configuration', self)
        load_config_action.triggered.connect(self.load_config)
        file_menu.addAction(load_config_action)
        
        TwoDCalc_action = QAction('2D Undulator Calc', self)
        TwoDCalc_action.triggered.connect(self.create2Dcalc)
        Util_menu.addAction(TwoDCalc_action)        
        
        
         
        main_layout = QVBoxLayout()
        main_layout.setMenuBar(menu_bar)
        # hboxlabel = QHBoxLayout()        
        # hboxlabel.addWidget(QLabel('------ Spectrum ------'))
        
        # main_layout.addLayout(hboxlabel)
        hboxErange = QHBoxLayout()
        hboxErange.addWidget(QLabel('E_min:'))
        self.Emin_edit = QLineEdit()
        self.Emin_edit.setPlaceholderText('Minimum Energy (eV)')
        self.Emin_edit.setText(str(4000.0))
        hboxErange.addWidget(self.Emin_edit)
        
        hboxErange.addWidget(QLabel('E_max:'))
        self.Emax_edit = QLineEdit()
        self.Emax_edit.setPlaceholderText('Maximum Energy (eV)')
        self.Emax_edit.setText(str(100000.0))
        hboxErange.addWidget(self.Emax_edit)
        
        hboxErange.addWidget(QLabel('# Points:'))
        self.NPts_edit = QLineEdit()
        self.NPts_edit.setText(str(500))
        self.NPts_edit.setPlaceholderText("# of E Pts in Spectrum")
        hboxErange.addWidget(self.NPts_edit)
        
        main_layout.addLayout(hboxErange)
         
        hbox1 = QHBoxLayout()
        self.APSUcheck = QCheckBox('APS-U')
        self.APSUcheck.setCheckState(2)
        self.APSUcheck.stateChanged.connect(self.APSUcheckchanged)
        hbox1.addWidget(self.APSUcheck)
        
   
        #  Undulator combo box
        self.combo = QComboBox()
        self.combo.addItems(["U23", "U14"])
        defaultundulator = 'U14'
        self.combo.setCurrentText(defaultundulator)
        self.combo.currentIndexChanged.connect(self.undulatorchanged)
        hbox1.addWidget(self.combo)
    
        self.harmonic_edit = QLineEdit()
        self.harmonic_edit.setPlaceholderText("Harmonic Number")
        self.harmonic_edit.setText(str(1))
        hbox1.addWidget(self.harmonic_edit)
    
        self.energy_edit = QLineEdit()
        self.energy_edit.setPlaceholderText("Energy (eV)")
        self.energy_edit.setText(str(22160.0))
        hbox1.addWidget(self.energy_edit)
         
         
        self.Go_button = QPushButton('Calc')
        hbox1.addWidget(self.Go_button)
        self.Go_button.clicked.connect(self.CalcSpectrum)
   
        main_layout.addLayout(hbox1)
        
        # line for setting K value
        hbox3 = QHBoxLayout()
        main_layout.addLayout(hbox3) 
        self.useKcheckbox = QCheckBox()
        self.useKcheckbox.setCheckState(0)
        hbox3.addWidget(QLabel('Use K:'))
        hbox3.addWidget(self.useKcheckbox)
        hbox3.addWidget(QLabel('K Value:'))
        self.KVal_edit = QLineEdit()
        self.KVal_edit.setText('0.436')
        hbox3.addWidget(self.KVal_edit)
        self.KVal_edit.setPlaceholderText('Set Undulator K-Value With Optimize Button, U23 Closed Gap = 1.498')
        
        #SetK button
        self.OptimizeK_button = QPushButton('Optimize K')
        hbox3.addWidget(self.OptimizeK_button)
        self.OptimizeK_button.clicked.connect(self.optimizeKclicked)  
         
        hbox2 = QHBoxLayout()
        
        formlayout = QFormLayout()
        main_layout.addLayout(formlayout)
        self.ringcurrent_edit = QLineEdit()
        formlayout.addRow(QLabel('Ring Current (mA)'),self.ringcurrent_edit)
        self.ringcurrent_edit.setPlaceholderText('Ring Current (mA)')
        self.ringcurrent_edit.setText(str(120.0))
        
    
        self.Mirror_button = QPushButton("Mirrors")
        hbox2.addWidget(self.Mirror_button)
        self.Mirror_button.clicked.connect(self.makeMirrorGUI)
         
        self.multilayer_button = QPushButton("Multilayer")
        hbox2.addWidget(self.multilayer_button)
        self.multilayer_button.clicked.connect(self.openMultilayerConfig)
         
        self.filter_button = QPushButton("Filters")
        hbox2.addWidget(self.filter_button)
        self.filter_button.clicked.connect(self.makefilterGUI)
   
        main_layout.addLayout(hbox2)

        hbox4 = QHBoxLayout()
        main_layout.addLayout(hbox4)        
        self.statustext = QLabel('Ready...')
        hbox4.addWidget(self.statustext)
   
        
        
    
        self.setLayout(main_layout)
    
        self.setWindowTitle("DCS Spectrum Calculation")
        self.setup_tray_icon()
        self.setGeometry(100, 100, 400, 200)
        
        #### Other Parameters
        self.SpectrumManager = None
        self.mirrormanager = []
        self.filtermanager =[]
        self.multilayerconfig = []
        self.multilayerenabled = 0
        self.subwindows = []
        self.spectralpower_dat=[]
        self.energy_dat=[]
        self.ringenergy=6.0
        self.WBScenterH = 0.0
        self.WBScenterV = 0.0
        self.WBSsizeH = 1.5
        self.WBSsizeV = 1.0
        self.esliceviewer = None
        
        
        self.slitparams =  ('WBSsizeH','WBSsizeV','WBScenterH','WBScenterV')
        

    def setup_tray_icon(self):
        if not QSystemTrayIcon.isSystemTrayAvailable():
            return

        self.tray_icon = QSystemTrayIcon(self)
        
        # Load icon 
        icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ui','DCSSpectrumIcon4.ico')
        if os.path.exists(icon_path):
            icon = QIcon(icon_path)
            self.tray_icon.setIcon(icon)
            self.setWindowIcon(icon)
        
        # Context Menu
        tray_menu = QMenu()
        
        restore_action = QAction("Restore", self)
        restore_action.triggered.connect(self.show_normal_and_raise)
        tray_menu.addAction(restore_action)
        
        quit_action = QAction("Quit", self)
        quit_action.triggered.connect(qApp.quit)
        tray_menu.addAction(quit_action)
        
        self.tray_icon.setContextMenu(tray_menu)
        self.tray_icon.show()
        
        self.tray_icon.activated.connect(self.on_tray_icon_activated)

    def show_normal_and_raise(self):
        self.showNormal()
        self.raise_()
        self.activateWindow()

    def on_tray_icon_activated(self, reason):
        if reason == QtWidgets.QSystemTrayIcon.Trigger:
            self.show_normal_and_raise() 
    def save_config(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getSaveFileName(self,"Save Configuration","","HDF5 Files (*.h5)", options=options)
        if fileName:
            try:
                with h5py.File(fileName, 'w') as f:
                    self.save_to_h5_group(f)
                QMessageBox.information(self, "Success", "Configuration saved successfully.")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save configuration: {e}")

    def save_to_h5_group(self, f):
        # Main Window Parameters
        f.attrs['Emin'] = self.Emin_edit.text()
        f.attrs['Emax'] = self.Emax_edit.text()
        f.attrs['NPts'] = self.NPts_edit.text()
        f.attrs['APSU'] = self.APSUcheck.checkState()
        f.attrs['Undulator'] = self.combo.currentText()
        f.attrs['Harmonic'] = self.harmonic_edit.text()
        f.attrs['Energy'] = self.energy_edit.text()
        f.attrs['UseK'] = self.useKcheckbox.checkState()
        f.attrs['KValue'] = self.KVal_edit.text()
        f.attrs['RingCurrent'] = self.ringcurrent_edit.text()
        
        # WBS Parameters
        f.attrs['WBSsizeH'] = self.WBSsizeH
        f.attrs['WBSsizeV'] = self.WBSsizeV
        f.attrs['WBScenterH'] = self.WBScenterH
        f.attrs['WBScenterV'] = self.WBScenterV

        # Mirror Manager
        # Always create the group to define the state (even if empty)
        grp = f.create_group("MirrorManager")
        
        m_stripes = []
        m_angles = []
        m_densities = []
        
        if self.mirrormanager:
             self.mirrormanager.update_lists()
             m_stripes = self.mirrormanager.activeMirrorStripes
             m_angles = self.mirrormanager.activeAngles
             m_densities = self.mirrormanager.activeDensities
        
        grp.create_dataset("activeMirrorStripes", data=np.array(m_stripes, dtype='S'))
        grp.create_dataset("activeAngles", data=np.array(m_angles))
        
        densities_str = [str(d) for d in m_densities]
        grp.create_dataset("activeDensities", data=np.array(densities_str, dtype='S'))

        # Filter Manager
        grp = f.create_group("FilterManager")
        f_compounds = []
        f_thicknesses = []
        f_densities = []

        if self.filtermanager:
            self.filtermanager.update_lists()
            f_compounds = self.filtermanager.activeCompounds
            f_thicknesses = self.filtermanager.activeThicknesses
            f_densities = self.filtermanager.activeDensities
        
        grp.create_dataset("activeCompounds", data=np.array(f_compounds, dtype='S'))
        grp.create_dataset("activeThicknesses", data=np.array(f_thicknesses))
        densities_str = [str(d) for d in f_densities]
        grp.create_dataset("activeDensities", data=np.array(densities_str, dtype='S'))

        # Multilayer Config
        grp = f.create_group("MultilayerConfig")
        if self.multilayerconfig:
            grp.attrs['enabled'] = self.multilayerenabled
            grp.attrs['ML_type'] = self.multilayerconfig.combo.currentText()
            grp.attrs['ML_energy'] = self.multilayerconfig.MLenergyedit.text()
            grp.attrs['ML_angle'] = self.multilayerconfig.angleedit.text()
            grp.attrs['checkbox_state'] = self.multilayerconfig.enableMLcheck.checkState()
        else:
            grp.attrs['enabled'] = 0
            grp.attrs['ML_type'] = '[None]'
            grp.attrs['ML_energy'] = '0'
            grp.attrs['ML_angle'] = '0'
            grp.attrs['checkbox_state'] = 0

    def scan_config(self, f):
        found = []
        # Handle nesting
        if 'SpectrumManager' in f:
             f = f['SpectrumManager']
             
        # Check Main Params (heuristic: check for a few key attrs)
        if 'Emin' in f.attrs or 'Undulator' in f.attrs:
            found.append("Main Parameters")
            
        if "MirrorManager" in f:
            found.append("Mirrors")
            
        if "FilterManager" in f:
            found.append("Filters")
            
        if "MultilayerConfig" in f:
            found.append("Multilayer")
            
        return found

    def load_config(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"Load Configuration","","HDF5 Files (*.h5)", options=options)
        if fileName:
            try:
                with h5py.File(fileName, 'r') as f:
                    available = self.scan_config(f)
                    if not available:
                         QMessageBox.information(self, "Info", "No recognized configuration data found.")
                         return
                         
                    dlg = ConfigImportDialog.ConfigImportDialog(available, self)
                    if dlg.exec_():
                        selected = dlg.get_selected()
                        summary = self.load_from_h5_group(f, selected)
                        
                        msg = "Configuration loaded successfully."
                        if summary:
                            msg += f"\n\nDetails:\n{summary}"
                        QMessageBox.information(self, "Success", msg)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load configuration: {e}")

    def load_from_h5_group(self, f, selected_components=None):
        if selected_components is None:
            # Default to all if not specified (legacy or direct call)
            selected_components = ["Main Parameters", "Mirrors", "Filters", "Multilayer"]
            
        loaded_components = []
        # Check if we are loading from a file that has the config in a subgroup 'SpectrumManager'
        # or if 'f' itself is the group/file containing the config
        
        # If f is a file or group, check if it has 'SpectrumManager' subgroup
        if 'SpectrumManager' in f:
            f = f['SpectrumManager']
        
        # Main Window Parameters
        if "Main Parameters" in selected_components:
            if 'Emin' in f.attrs: 
                self.Emin_edit.setText(str(f.attrs['Emin']))
                loaded_components.append("Main Parameters")
            if 'Emax' in f.attrs: self.Emax_edit.setText(str(f.attrs['Emax']))
            if 'NPts' in f.attrs: self.NPts_edit.setText(str(f.attrs['NPts']))
            if 'APSU' in f.attrs: self.APSUcheck.setCheckState(int(f.attrs['APSU']))
            if 'Undulator' in f.attrs: self.combo.setCurrentText(f.attrs['Undulator'])
            if 'Harmonic' in f.attrs: self.harmonic_edit.setText(str(f.attrs['Harmonic']))
            if 'Energy' in f.attrs: self.energy_edit.setText(str(f.attrs['Energy']))
            if 'UseK' in f.attrs: self.useKcheckbox.setCheckState(int(f.attrs['UseK']))
            if 'KValue' in f.attrs: self.KVal_edit.setText(str(f.attrs['KValue']))
            if 'RingCurrent' in f.attrs: self.ringcurrent_edit.setText(str(f.attrs['RingCurrent']))
            
            if 'WBSsizeH' in f.attrs: self.WBSsizeH = f.attrs['WBSsizeH']
            if 'WBSsizeV' in f.attrs: self.WBSsizeV = f.attrs['WBSsizeV']
            if 'WBScenterH' in f.attrs: self.WBScenterH = f.attrs['WBScenterH']
            if 'WBScenterV' in f.attrs: self.WBScenterV = f.attrs['WBScenterV']


        # Mirror Manager - Clear first if exists
        if "Mirrors" in selected_components:
            if self.mirrormanager:
                while self.mirrormanager.Mirrors_container.count() > 0:
                    item = self.mirrormanager.Mirrors_container.itemAt(0)
                    widget = item.widget()
                    if widget:
                        self.mirrormanager.remove_row(widget)
                self.mirrormanager.update_lists() # Clear internal lists

            if "MirrorManager" in f:
                print('Found MirrorManager in config file')
                grp = f["MirrorManager"]
                
                stripes = []
                if "activeMirrorStripes" in grp:
                    raw_stripes = grp["activeMirrorStripes"][:]
                    for s in raw_stripes:
                        if hasattr(s, 'decode'):
                            s = s.decode('utf-8')
                        else:
                            s = str(s)
                        stripes.append(s.strip())
                
                angles = []
                if "activeAngles" in grp:
                    angles = grp["activeAngles"][:]
                
                densities = []
                if "activeDensities" in grp:
                    raw_densities = grp["activeDensities"][:]
                    for s in raw_densities:
                        if hasattr(s, 'decode'):
                            s = s.decode('utf-8')
                        else:
                            s = str(s)
                        densities.append(s.strip())
                
                print(f'Loading Mirrors: {len(stripes)} found.')
                
                if len(stripes) > 0:
                    self.makeMirrorGUI() # Ensure GUI exists and show it if we have mirrors
                    
                    # Ensure lengths match or notify
                    if not (len(stripes) == len(angles) == len(densities)):
                        print(f"Warning: properties length mismatch: S:{len(stripes)} A:{len(angles)} D:{len(densities)}")
                    
                    for s, a, d in zip(stripes, angles, densities):
                        if hasattr(a, 'item'):
                            a = a.item()
                        print(f"DEBUG ADD MIRROR: Stripe='{s}' Angle={a} Density='{d}'")
                        self.mirrormanager.add_Mirror(s, a, d, enabled=True)
                    self.mirrormanager.update_lists()
                    loaded_components.append(f"Mirrors: {len(stripes)} loaded")
                else: 
                     # Mirrors were selected, but none found in file.
                     # We still cleared the UI, so effectively we loaded "0 mirrors"
                     loaded_components.append("Mirrors: 0 loaded (Cleared)")


        # Filter Manager
        if "Filters" in selected_components:
            if "FilterManager" in f:
                self.makefilterGUI()
                while self.filtermanager.filters_container.count() > 0:
                    item = self.filtermanager.filters_container.itemAt(0)
                    widget = item.widget()
                    if widget:
                        self.filtermanager.remove_row(widget)
                
                grp = f["FilterManager"]
                compounds = [s.decode('utf-8') for s in grp["activeCompounds"][:]]
                thicknesses = grp["activeThicknesses"][:]
                densities = [s.decode('utf-8') for s in grp["activeDensities"][:]]
                
                for c, t, d in zip(compounds, thicknesses, densities):
                    self.filtermanager.add_filter(c, t, d, enable=True)
                self.filtermanager.update_lists()
                loaded_components.append(f"Filters: {len(compounds)} loaded")

        # Multilayer Config
        if "Multilayer" in selected_components:
            if "MultilayerConfig" in f:
                self.openMultilayerConfig()
                grp = f["MultilayerConfig"]
                if 'ML_type' in grp.attrs: self.multilayerconfig.combo.setCurrentText(grp.attrs['ML_type'])
                if 'ML_energy' in grp.attrs: self.multilayerconfig.MLenergyedit.setText(str(grp.attrs['ML_energy']))
                if 'ML_angle' in grp.attrs: self.multilayerconfig.angleedit.setText(str(grp.attrs['ML_angle']))
                if 'checkbox_state' in grp.attrs: self.multilayerconfig.enableMLcheck.setCheckState(int(grp.attrs['checkbox_state']))
                if 'enabled' in grp.attrs: self.multilayerenabled = grp.attrs['enabled']
                loaded_components.append("Multilayer Configuration")

        return "\n".join(loaded_components)

    def undulatorchanged(self):
        val = self.combo.currentText()
        if val == 'U14':
            self.harmonic_edit.setText(str(1))
            self.energy_edit.setText(str(22160.0))
            self.KVal_edit.setText('0.436')

        elif val == 'U23':
            self.harmonic_edit.setText(str(3))
            self.energy_edit.setText(str(20925))
            self.KVal_edit.setText('1.498')
                
    def open_gui_to_adjust_WBS(self):
        """Opens a GUI to adjust specific parameters."""
        self.subwindows.append(ParameterEditor(self, self.slitparams))
        self.subwindows[-1].show()   
        
    def APSUcheckchanged(self):
        if self.APSUcheck.checkState()<2:
            self.combo.clear()
            self.combo.addItems(["U27", "U17"])
            self.ringcurrent_edit.setText(str(100.0))
            self.ringenergy=7.0
        else:
            self.combo.clear()
            self.combo.addItems(["U23", "U14"])    
            self.ringcurrent_edit.setText(str(120.0))
            self.ringenergy=6.0

    def create2Dcalc(self):
        self.twodcalcGUI = PlotterWidget(self)
        self.twodcalcGUI.show()            
            
        
    def closeEvent(self, event):
        for wins in self.subwindows:
            if wins:
                wins.close()                 
        
    def openMultilayerConfig(self):
        if not self.multilayerconfig:
            if self.energy_edit.text():
                self.multilayerconfig = MultilayerConfig.MultilayerConfig(self,float(self.energy_edit.text()))
            else:
                self.multilayerconfig = MultilayerConfig.MultilayerConfig(self)
        self.multilayerconfig.resize(200, 200)
        self.multilayerconfig.show()
        
    def getcurrentscreenvalues(self):
       undulatorvalue = self.combo.currentText()
       
       if self.APSUcheck.checkState()>0: 
           if undulatorvalue=='U23':
               period = 0.023
               minE = 5000
               maxE = 14000
           elif undulatorvalue=='U14':
               period = 0.014
               minE = 19000
               maxE = 23000
       else:
           if undulatorvalue=='U27':
               period = 0.027
               minE = 5000
               maxE = 14000
           elif undulatorvalue=='U17':
               period = 0.017
               minE = 19000
               maxE = 24000 
       harval = self.harmonic_edit.text() 
       if harval:
           try:
               harmonicnumber = int(harval)
           except:
               harmonicnumber = []
               print('Invalid Inputs: Harmonic Number !')
               return 0,0,0,0,0,0,0
       else:
           print('Required Input: Harmonic Number.')
           print(type(harval))
           return 0,0,0,0,0,0,0            
       
       harmonicenergy = self.energy_edit.text()
    
       if harmonicenergy:
            try:
                energy = float(harmonicenergy)
            except ValueError:
                print('Invalid Inputs: harmonic energy')
                return 0,0,0,0,0,0,0
       else:
                print('Required Input: harmonic energy')
                return 0,0,0,0,0,0,0
       
       if self.useKcheckbox.checkState()<2:
           if energy<(minE*float(harval)):
               print('Harmonic Energy too low')
               self.statustext.setText('Harmonic Energy too low')
               return 0,0,0,0,0,0,0
           if energy>(maxE*float(harval)):
               print('Harmonic Energy too high')
               self.statustext.setText('Harmonic Energy too high')
               return 0,0,0,0,0,0,0
       
       npts = int(self.NPts_edit.text())
       emin = float(self.Emin_edit.text())
       emax = float(self.Emax_edit.text())
       return undulatorvalue,period,harmonicnumber,energy,npts,emin,emax
   
    def CalcSpectrum(self): 
        undulatorvalue,period,harmonicnumber,energy,npts,emin,emax = self.getcurrentscreenvalues()
        
        if undulatorvalue==0:
            return
        slitsize_h = self.WBSsizeH#0.001457*1000
        slitsize_v= self.WBSsizeV#0.000973*1000
        slitcenter_h = self.WBScenterH
        slitcenter_v = self.WBScenterV
        
        if self.APSUcheck.checkState()>0:
            APSU=1
            offset = 1.00375
        else:
            APSU=0
            offset = 1
                    
        ringcurrent = float(self.ringcurrent_edit.text())
        
        if self.useKcheckbox.checkState()==2:
            K = float(self.KVal_edit.text())
        else:
            offsetenergy = energy*offset
            K = self.KfromRingEnergyAndUndulatorPeriod(self.ringenergy,period,offsetenergy,harmonicnumber)
             
        UndulatorConfig = UndulatorCalcConfig(self,undulatorvalue,period,harmonicnumber,energy,npts,emin,emax,K,ringcurrent,APSU)
        
        edat,specpower0 = undulatorcalc(UndulatorConfig,slitcenter_h,slitcenter_v,slitsize_h,slitsize_v)
             
        self.subwindows.append(ExtraFigureWindow(self))
        
        pltlist=[]
        newfig = self.subwindows[-1]        
        newfig.ax.set_aspect('auto')
        newfig.ax.plot(edat,specpower0,label='Undulator Output')
           

      
        #Mirrors
        if self.mirrormanager:
            if len(self.mirrormanager.activeMirrorStripes)>0:
                edat,specpower1=self.XrayMirror(edat,specpower0)
                newfig.ax.plot(edat,specpower1,label='After Mirrors')
            else:
                specpower1 = specpower0
        else:
            specpower1 = specpower0
    
            
            
        #Filters
        if self.filtermanager:
            if len(self.filtermanager.activeCompounds)>0:
                edat,specpower2=self.XrayFilter(edat,specpower1)
                newfig.ax.plot(edat,specpower2,label='After Filters')
                # newfig.xdata.append(edat)
                # newfig.ydata.append(specpower2)
            else:
                specpower2 = specpower1
        else:
            specpower2 = specpower1            
       
        

        if self.multilayerenabled:
            
            try:
                MLangle = float(self.multilayerconfig.angleedit.text())
            except:
                print('Invalid Value For ML Angle')
                self.statustext.setText('Invalid Value For ML Angle')
                return 
            
            if not self.multilayerconfig.MLReflectivityInterpolator:
                self.multilayerconfig.readMLData()
                
            reflectivityvals = self.multilayerconfig.MLReflectivityInterpolator((MLangle,edat))
            specpower3 = specpower2*reflectivityvals 
            newfig.ax.plot(edat,specpower3,linestyle='solid',label='After Multilayer')

        else:
            specpower3 = specpower2
            
        lines = newfig.ax.lines
            
        if len(lines)>1:
            for l in lines[:-1]:
                l.set(linestyle='--')
            
        newfig.ax.set_xlabel('Energy (eV)', fontsize=8)
        newfig.ax.set_ylabel('Spectral Power (W/eV)', fontsize=8)
        newfig.ax.tick_params(axis='both', labelsize=8)
        newfig.ax.legend(fontsize=8)    
        newfig.show() 

        
        newfig.xdata.append(edat)
        newfig.ydata.append(specpower3)
        self.spectralpower_dat = specpower3
        self.energy_dat = edat
        
        
        photons_pereV_persec = specpower3/edat/(1.602e-19)
        photons_persec = trapezoid(photons_pereV_persec,edat)
        photons_perbunch24bunch = trapezoid(photons_pereV_persec,edat)/(271647*24)
        photons_perbunch48bunch = trapezoid(photons_pereV_persec,edat)/(271647*48)
        totalpower = trapezoid(self.spectralpower_dat,self.energy_dat)
        
        Eperbunch_24_uJ = totalpower/(271647*24)*1e6
        Eperbunch_48_uJ = totalpower/(271647*48)*1e6
        # (70)/(24000)/(6.52e6)/(1.602e-19)
        self.statustext.setText(f'Power: {totalpower:.3f} W, 24-bunch ph.: {photons_perbunch24bunch:.2e} ph. ({Eperbunch_24_uJ:.2f} µJ), 48-bunch ph.: {photons_perbunch48bunch:.2e} ph. ({Eperbunch_48_uJ:.2f} µJ)')
         
        
    def save_spectrum_toCSV(self):
        if self.energy_dat is not None:
            options = QFileDialog.Options()
            file_name, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
    
            # If the user selects a file, write the CSV
            if file_name:
                try:
                    self.write_csv(file_name,self.energy_dat,self.spectralpower_dat)
                    QMessageBox.information(self, "Success", f"CSV file saved as: {file_name}")
    
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to save CSV file: {e}")
        else: 
            QMessageBox.critical(self, "Error", "No Data Yet")

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
    
    def addstandardfilters(self):
        self.filtermanager.add_filter("C", 0.2,3.53,enable=True)
        
        self.filtermanager.add_filter("Be", 0.508+0.508+0.254,'?',enable=True)
        self.filtermanager.add_filter("Kapton Polyimide Film",0.25, '?',enable=True)
        self.filtermanager.add_filter("He", 300,'?',enable=True)
        self.filtermanager.add_filter("C0.000124N0.755268O0.231781Ar0.012827", 100,0.00120479,enable=True)
        

    def makefilterGUI(self):
        if not self.filtermanager:
            self.filtermanager = FilterGUI.FilterManager(self)
            self.filtermanager.resize(600, 200)
            self.filtermanager.show()
            self.addstandardfilters()
        else:
            self.filtermanager.show()

    def makeMirrorGUI(self):
        if not self.mirrormanager:
            self.mirrormanager = MirrorGUI.MirrorManager(self)
            self.mirrormanager.resize(600, 200)
            self.mirrormanager.show()
            self.mirrormanager.add_Mirror("Pt", 2.1,'?')
            self.mirrormanager.add_Mirror("Pt", 2.1,'?')
            self.mirrormanager.add_Mirror("Rh", 2.1,'?')
            self.mirrormanager.add_Mirror("Rh", 2.1,'?')
            self.mirrormanager.add_Mirror("Si", 2.1,'?')
            self.mirrormanager.add_Mirror("Si", 2.1,'?')

        else:
            self.mirrormanager.show()            
        
        
    def KfromRingEnergyAndUndulatorPeriod(self,ringEnergyGeV,periodlength_inm,photon_energy_eV,harmonic=1):
        m2ev = codata.c * codata.h / codata.e
        wavelength = harmonic * m2ev / photon_energy_eV
        gamma = 1e9*ringEnergyGeV/(510998.95) 
        val = 2 * ((wavelength * 2 * gamma**2 / periodlength_inm) - 1)
        if val>0:
            K =np.sqrt(val)
        else:
            K = []
        return K
    


    
    def initialguessMLAngle(self,energyeV,dspacing_nm):
        return np.asin(1239.84198/(energyeV)/2/(dspacing_nm)) 
    
    def XrayFilter(self,energy,spectral_power):
        
        out_dictionary = xoppy_calc_power(
                energy,
                spectral_power,
                substance = self.filtermanager.activeCompounds,
                thick     = self.filtermanager.activeThicknesses,
                angle     = np.zeros(len(self.filtermanager.activeThicknesses)), # in mrad (for mirrors)
                dens      = self.filtermanager.activeDensities,
                roughness = np.zeros(len(self.filtermanager.activeThicknesses)), # in A (for mirrors)
                flags     = np.zeros(len(self.filtermanager.activeThicknesses)), # 0=Filter, 1=Mirror
                nelements = len(self.filtermanager.activeCompounds),
                FILE_DUMP = 0,
                material_constants_library = xraylib,
                )

 
        newenergy = out_dictionary["data"][0,:]
        newspectral_power = out_dictionary["data"][-1,:]
        return newenergy,newspectral_power
    
    def XrayMirror(self,energy,spectral_power):
        
        out_dictionary = xoppy_calc_power(
                energy,
                spectral_power,
                substance = self.mirrormanager.activeMirrorStripes,
                thick     = np.zeros(len(self.mirrormanager.activeAngles)),
                angle     = self.mirrormanager.activeAngles, # in mrad (for mirrors)
                dens      = self.mirrormanager.activeDensities,
                roughness = np.zeros(len(self.mirrormanager.activeAngles)), # in A (for mirrors)
                flags     = np.ones(len(self.mirrormanager.activeAngles)), # 0=Filter, 1=Mirror
                nelements = len(self.mirrormanager.activeMirrorStripes),
                FILE_DUMP = 0,
                material_constants_library = xraylib,
                )
 
        newenergy = out_dictionary["data"][0,:]
        newspectral_power = out_dictionary["data"][-1,:]
        return newenergy,newspectral_power
    
        

        
    def plotinnewwindow(self,x,y):
        self.subwindows.append(ExtraFigureWindow(self))
        newfig = self.subwindows[-1]
        newfig.ax.plot(x,y)
        newfig.ax.set_aspect('auto')
        newfig.show() 
        
    def plot_filter_attenuation(self, energy, attenuation):
        """Plot transmission of the current filter stack over the requested energy range."""
        self.subwindows.append(ExtraFigureWindow(self))
        newfig = self.subwindows[-1]
        newfig.ax.set_aspect('auto')
        newfig.ax.plot(energy, attenuation, label='Filter Attenuation')
        newfig.ax.set_xlabel('Energy (eV)')
        newfig.ax.set_ylabel('Transmission')
        newfig.ax.legend()
        newfig.show()
        newfig.xdata.append(energy)
        newfig.ydata.append(attenuation)
        
    def optimizeKclicked(self):
        undulatorvalue,period_inm,harmonic,DesiredPeakPhotonEnergy,npts,emin,emax = self.getcurrentscreenvalues()
        if undulatorvalue==0:
            return
        if self.APSUcheck.checkState()>0:

            if period_inm==0.023:
                numperiods = 206
            elif period_inm==0.014:
                numperiods = 338
            else: 
                undulatorlength = 4.738 #meters
                numperiods = undulatorlength/period_inm        
            offsetenergy = DesiredPeakPhotonEnergy*1.00375
            Kstart = self.KfromRingEnergyAndUndulatorPeriod(self.ringenergy,period_inm,offsetenergy,harmonic)
            self.KVal_edit.setText('Working...')
            self.show()
            K = self.Undulator_findpeakK_APSU(period_inm,DesiredPeakPhotonEnergy,harmonic)       
        else:
            
            if period_inm==0.027:
                numperiods = 88
            elif period_inm==0.017:
                numperiods = 137
            else: 
                undulatorlength = 0.027*88 #meters
                numperiods = undulatorlength/period_inm    
            offsetenergy = DesiredPeakPhotonEnergy*1.00375
            Kstart = self.KfromRingEnergyAndUndulatorPeriod(self.ringenergy,period_inm,offsetenergy,harmonic)
            self.KVal_edit.setText('Working...')
            self.show()
            K = self.Undulator_findpeakK_APS(period_inm,DesiredPeakPhotonEnergy,harmonic)   
        print(K)
        self.KVal_edit.setText(f'{K:.5f}')
        if K:
            self.useKcheckbox.setCheckState(2)
  
    
    def xoppy_calc_undulator_spectrum_streamline(self,ELECTRONENERGY=6.04,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
                                  ELECTRONBEAMSIZEH=0.000395,ELECTRONBEAMSIZEV=9.9e-06,\
                                  ELECTRONBEAMDIVERGENCEH=1.05e-05,ELECTRONBEAMDIVERGENCEV=3.9e-06,\
                                  PERIODID=0.018,NPERIODS=222,KV=1.68,KH=0.0,KPHASE=0.0,DISTANCE=30.0,
                                  GAPH=0.001,GAPV=0.001,GAPH_CENTER=0.0,GAPV_CENTER=0.0,\
                                  PHOTONENERGYMIN=3000.0,PHOTONENERGYMAX=55000.0,PHOTONENERGYPOINTS=500,
                                  USEEMITTANCES=1):
        METHOD=2
        # print("Inside xoppy_calc_undulator_spectrum. ")

        bl = OrderedDict()
        bl['ElectronBeamDivergenceH'] = ELECTRONBEAMDIVERGENCEH
        bl['ElectronBeamDivergenceV'] = ELECTRONBEAMDIVERGENCEV
        bl['ElectronBeamSizeH'] = ELECTRONBEAMSIZEH
        bl['ElectronBeamSizeV'] = ELECTRONBEAMSIZEV
        bl['ElectronCurrent'] = ELECTRONCURRENT
        bl['ElectronEnergy'] = ELECTRONENERGY
        bl['ElectronEnergySpread'] = ELECTRONENERGYSPREAD
        bl['Kv'] = KV
        bl['Kh'] = KH
        bl['Kphase'] = KPHASE
        bl['NPeriods'] = NPERIODS
        bl['PeriodID'] = PERIODID
        bl['distance'] = DISTANCE
        bl['gapH'] = GAPH
        bl['gapV'] = GAPV
        bl['gapHcenter'] = GAPH_CENTER
        bl['gapVcenter'] = GAPV_CENTER

        if USEEMITTANCES:
            zero_emittance = False
        else:
            zero_emittance = True

        outFile = None
        print(bl["PeriodID"])

        codata_mee = codata.m_e * codata.c**2 / codata.e # electron mass in eV
        gamma = bl['ElectronEnergy'] * 1e9 / codata_mee

        m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)
        resonance_wavelength = (1 + (bl['Kv']**2 + bl['Kh']**2) / 2.0) / 2 / gamma**2 * bl["PeriodID"]
        resonance_energy = m2ev / resonance_wavelength
        # print ("Gamma: %f \n"%(gamma))
        # print ("Resonance wavelength [A]: %g \n"%(1e10*resonance_wavelength))
        # print ("Resonance energy [eV]: %g \n"%(resonance_energy))


        ptot = (NPERIODS/6) * codata.value('characteristic impedance of vacuum') * \
               ELECTRONCURRENT * codata.e * 2 * np.pi * codata.c * gamma**2 * (KV**2+KH**2) / PERIODID
        # print ("\nTotal power radiated by the undulator with fully opened slits [W]: %g \n"%(ptot))


        
        if METHOD == 2:
            # get the maximum harmonic number
            h_max = int(2.5*PHOTONENERGYMAX/resonance_energy)
            e, f = calc1d_srw(bl,photonEnergyMin=PHOTONENERGYMIN,photonEnergyMax=PHOTONENERGYMAX,
                  photonEnergyPoints=PHOTONENERGYPOINTS,fileName=outFile,fileAppend=False,zero_emittance=zero_emittance,
                  srw_max_harmonic_number=h_max)


        # if zero_emittance:
            # print("\nNo emittance calculation")

        if METHOD == 1 and len(e) == 0: raise Exception("Invalid Input Parameters")

        # power_in_spectrum = f.sum()*1e3*codata.e*(e[1]-e[0])
        # print("\nPower from integral of spectrum: %8.3f W"%(power_in_spectrum))
        # print("\nRatio Power from integral of spectrum over Total emitted power: %5.4f"%(power_in_spectrum / ptot))

        spectral_power = f * codata.e * 1e3
        # try:
        #     cumulated_power = spectral_power.cumsum() * np.abs(e[0] - e[1])
        # except:
        #     cumulated_power = 0.0

        return e, f, spectral_power
    
    def Undulator_findpeakK_APSU(self,period_inm,DesiredPeakPhotonEnergy,harmonic):
        print('--------------Starting Optimization----------')
    
        ringenergy = 6.0 #GeV
        if period_inm==0.023:
            numperiods = 206
        elif period_inm==0.014:
            numperiods = 338
        else: 
            undulatorlength = 4.738 #meters
            numperiods = undulatorlength/period_inm        
        Kstart = self.KfromRingEnergyAndUndulatorPeriod(ringenergy,period_inm,DesiredPeakPhotonEnergy,harmonic)
        print(Kstart)
        K = Kstart
        # PhotonEnergyMin = DesiredPeakPhotonEnergy*0.95
        # PhotonEnergyMax = DesiredPeakPhotonEnergy*1.05
        numpts = 50
        Nitermax = 50
        Niter = 0
        ringcurr =float(self.ringcurrent_edit.text())/1000
        slitcenter_h = self.WBScenterH
        slitcenter_v = self.WBScenterV
        
        def getpeakE(Kval,EnergyCenter,erange,npts=numpts):
            energy, flux, spectral_power =self.xoppy_calc_undulator_spectrum_streamline(ELECTRONENERGY=ringenergy,
            ELECTRONENERGYSPREAD=0.00098,
            ELECTRONCURRENT=ringcurr,
            ELECTRONBEAMSIZEH=1.29e-05,
            ELECTRONBEAMSIZEV=2.5e-06,
            ELECTRONBEAMDIVERGENCEH=8.7e-06,
            ELECTRONBEAMDIVERGENCEV=3.6e-06,
            PERIODID=period_inm,
            NPERIODS=numperiods,
            KV=Kval,
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=31.0,
            GAPH=self.WBSsizeH*.001,#0.001457,
            GAPV=self.WBSsizeV*.001,#0.000973,
            GAPH_CENTER=slitcenter_h*.001,
            GAPV_CENTER=slitcenter_v*.001,
            PHOTONENERGYMIN=EnergyCenter-erange,
            PHOTONENERGYMAX=EnergyCenter+erange,
            PHOTONENERGYPOINTS=numpts,
            USEEMITTANCES=1)
            
            spl =scipy.interpolate.CubicSpline(energy, spectral_power)
            Enew = np.arange(EnergyCenter-erange,EnergyCenter+erange,1)
            pow_int = spl(Enew)
            # pow_deriv = spl(Enew,1)
            
            
            return Enew[np.argmax(pow_int)]
            
        

        
        # dEdK = -4000 #This isn't real, but it's a bit more robust than measuring it the way I was


        tolE = 2
        errorscaling = 1

        while True:
            if Niter ==0:
                #Find gradient dE/dK  
                dK = 0.002*Kstart/harmonic
                peakE = getpeakE(K,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy,100)
                errordE = peakE-DesiredPeakPhotonEnergy
                lastE = peakE
                print(f'Peak E:{peakE:.1f}, DesiredE: {DesiredPeakPhotonEnergy:.1f}, error: {errordE:.1f}')
                # peakE1 = getpeakE(K+dK,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy)
                # peakE2 = getpeakE(K-dK,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy)
                # dEdK = (peakE1-peakE2)/2/dK
                # print(f'E1 {peakE1}')
                # print(f'E2 {peakE2}')
                # print(f'dEdK measured {dEdK}')
                

            K = K-dK
            Erange = max([np.abs(errordE)*3,3000.0])
            # Erange = np.abs(errordE)*3
            
            npt= np.min([200,int(Erange)])
            peakE = getpeakE(K,DesiredPeakPhotonEnergy,Erange,npt)    
            errordE = peakE-DesiredPeakPhotonEnergy
            print(f'{Niter}: Erange: {Erange},{npt} pts, K={K:0.4f}')
            print(f'Peak E:{peakE:.1f}, DesiredE: {DesiredPeakPhotonEnergy:.1f}, error: {errordE:.1f}')
            
            if Niter>=Nitermax:
                print('MAX ITERATIONS EXCEEDED')
                return K
            
            if np.abs(errordE)<tolE:
                return K 
            else:
                
                dEdK = -(peakE - lastE)/dK
                if np.abs(dEdK)>0:
                    dK = (errordE/dEdK)*errorscaling
                    lastE = peakE
                    print(f'dE/dK:{dEdK:.1f}, K: {K:.5f}, dK: {(errordE/dEdK)*errorscaling:.5f}')
            Niter +=1
            
    # def set50keVconfig(self):
    #     undulatorconfig=dict()
    #     undulatorconfig['checkboxes'] = dict()
    #     undulatorconfig['editboxes'] = dict()
    #     undulatorconfig['main_window_values'] = dict()
    #     undulatorconfig['comboboxes'] = dict()
        
    #     undulatorconfig['comboboxes']['combo']='U23'
        
        
    #     undulatorconfig['checkboxes']['useKcheckbox']=1
    #     undulatorconfig['checkboxes']['APSUcheck']=1
        
        
    #     undulatorconfig['editboxes']['KVal_edit'] = 1.61
        
    #     hbox3.addWidget(self.KVal_edit)
    #     self.KVal_edit.setPlaceholderText('Set Undulator K-Value With Optimize Button, U23 Closed Gap = 1.61')
        
    #     #SetK button
    #     self.OptimizeK_button = QPushButton('Optimize K')
    #     hbox3.addWidget(self.OptimizeK_button)
    #     self.OptimizeK_button.clicked.connect(self.optimizeKclicked)  
        
    #     # Second horizontal layout for buttons
    #     hbox2 = QHBoxLayout()
        
    #     formlayout = QFormLayout()
    #     main_layout.addLayout(formlayout)
    #     self.ringcurrent_edit = QLineEdit()
    #     formlayout.addRow(QLabel('Ring Current (mA)'),self.ringcurrent_edit)
    #     self.ringcurrent_edit.setPlaceholderText('Ring Current (mA)')
    #     self.ringcurrent_edit.setText(str(200.0))
        
   
    #     # Mirror button
    #     self.Mirror_button = QPushButton("Mirrors")
    #     hbox2.addWidget(self.Mirror_button)
    #     self.Mirror_button.clicked.connect(self.makeMirrorGUI)
        
    #     # Multilayer button
    #     self.multilayer_button = QPushButton("Multilayer")
    #     hbox2.addWidget(self.multilayer_button)
    #     self.multilayer_button.clicked.connect(self.openMultilayerConfig)
        
    #     # Filter Manager button
    #     self.filter_button = QPushButton("Filters")
    #     hbox2.addWidget(self.filter_button)
    #     self.filter_button.clicked.connect(self.makefilterGUI)
   
    #     main_layout.addLayout(hbox2)

    #     hbox4 = QHBoxLayout()
    #     main_layout.addLayout(hbox4)        
    #     self.statustext = QLabel('Ready...')
    #     hbox4.addWidget(self.statustext)
   
        
        
   
    #     # Set the main layout for the widget
    #     self.setLayout(main_layout)
   
    #     # Window settings
    #     self.setWindowTitle("DCS Spectrum Calculation")
    #     self.setGeometry(100, 100, 400, 200)
        
    #     #### Other Parameters
    #     self.SpectrumManager = None
    #     self.mirrormanager = []
    #     self.filtermanager =[]
    #     self.multilayerconfig = []
    #     self.multilayerenabled = 0
    #     self.subwindows = []
    #     self.spectralpower_dat=[]
    #     self.energy_dat=[]
    #     self.ringenergy=6.0
    #     self.WBSsizeH = 1.5
    #     self.WBSsizeV = 1.0
        
        
        
        
        
        # self.setSpectrumConfig(filterconfig=None,mlconfig=None,mirrorconfig=None,undulatorconfig=undulatorconfig)
        

    def setSpectrumConfig(self,filterconfig=None,mlconfig=None,mirrorconfig=None,undulatorconfig=None):
        self.useKcheckbox.setCheckState(undulatorconfig.useKcheckbox)
        self.KVal_edit.setText(str(undulatorconfig['KVal_edit']))


    def Undulator_findpeakK_APS(self,period_inm,DesiredPeakPhotonEnergy,harmonic):
        print('--------------Starting Optimization----------')
    
        ringenergy = 7.0 #GeV
        if period_inm==0.027:
            numperiods = 88
        elif period_inm==0.017:
            numperiods = 137
        else: 
            undulatorlength = 0.027*88 #meters
            numperiods = undulatorlength/period_inm     
        Kstart = self.KfromRingEnergyAndUndulatorPeriod(ringenergy,period_inm,DesiredPeakPhotonEnergy,harmonic)
        K = Kstart
        # PhotonEnergyMin = DesiredPeakPhotonEnergy*0.95
        # PhotonEnergyMax = DesiredPeakPhotonEnergy*1.05
        numpts = 50
        Nitermax = 50
        Niter = 0
        ringcurr =float(self.ringcurrent_edit.text())/1000
        
        def getpeakE(Kval,EnergyCenter,erange,npts=numpts):
            energy, flux, spectral_power =self.xoppy_calc_undulator_spectrum_streamline(ELECTRONENERGY=ringenergy,
            ELECTRONENERGYSPREAD=0.00098,
            ELECTRONCURRENT=ringcurr,
            ELECTRONBEAMSIZEH=0.0002805,
            ELECTRONBEAMSIZEV=1.02e-05,
            ELECTRONBEAMDIVERGENCEH=1.18e-05,
            ELECTRONBEAMDIVERGENCEV=3.4e-06,
            PERIODID=period_inm,
            NPERIODS=numperiods,
            KV=Kval,
            KH=0.0,
            KPHASE=0.0,
            DISTANCE=31.0,
            GAPH=self.WBSsizeH*.001,#0.001457,
            GAPV=self.WBSsizeV*.001,#0.000973,
            GAPH_CENTER=0.0,
            GAPV_CENTER=0.0,
            PHOTONENERGYMIN=EnergyCenter-erange,
            PHOTONENERGYMAX=EnergyCenter+erange,
            PHOTONENERGYPOINTS=npts,
            USEEMITTANCES=1)
            spl =scipy.interpolate.CubicSpline(energy, spectral_power)
            Enew = np.arange(EnergyCenter-erange,EnergyCenter+erange,1)
            pow_int = spl(Enew)
            return Enew[np.argmax(pow_int)]
            
        tolE = 2
        errorscaling = 1

        while True:
            if Niter ==0:
                #Find gradient dE/dK  
                dK = 0.002*Kstart/harmonic
                peakE = getpeakE(K,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy,100)
                errordE = peakE-DesiredPeakPhotonEnergy
                lastE = peakE
                print(f'Peak E:{peakE:.1f}, DesiredE: {DesiredPeakPhotonEnergy:.1f}, error: {errordE:.1f}')
                # peakE1 = getpeakE(K+dK,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy)
                # peakE2 = getpeakE(K-dK,DesiredPeakPhotonEnergy,0.1*DesiredPeakPhotonEnergy)
                # dEdK = (peakE1-peakE2)/2/dK
                # print(f'E1 {peakE1}')
                # print(f'E2 {peakE2}')
                # print(f'dEdK measured {dEdK}')
                

            K = K-dK
            Erange = max([np.abs(errordE)*3,3000.0])
            # Erange = np.abs(errordE)*3
            
            npt= np.min([200,int(Erange)])
            peakE = getpeakE(K,DesiredPeakPhotonEnergy,Erange,npt)    
            errordE = peakE-DesiredPeakPhotonEnergy
            print(f'{Niter}: Erange: {Erange},{npt} pts, K={K:0.4f}')
            print(f'Peak E:{peakE:.1f}, DesiredE: {DesiredPeakPhotonEnergy:.1f}, error: {errordE:.1f}')
            
            if Niter>=Nitermax:
                print('MAX ITERATIONS EXCEEDED')
                return K
            
            if np.abs(errordE)<tolE:
                return K 
            else:
                
                dEdK = -(peakE - lastE)/dK
                if np.abs(dEdK)>0:
                    dK = (errordE/dEdK)*errorscaling
                    lastE = peakE
                    print(f'dE/dK:{dEdK:.1f}, K: {K:.5f}, dK: {(errordE/dEdK)*errorscaling:.5f}')
            Niter +=1
            
class ParameterEditor(QWidget):
    def __init__(self, parent_obj, params_to_edit):
        super().__init__()
        
        self.parent_obj = parent_obj
        self.params_to_edit = params_to_edit
        
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Edit Params")
        layout = QVBoxLayout()

        self.edit_boxes = {}
        

        #   create labels and edit boxes for each parameter
        for param in self.params_to_edit:
            # print(param)
            if hasattr(self.parent_obj, param):
                layout.addWidget(QLabel(param))
                currval = getattr(self.parent_obj, param) 
                edit_box = QLineEdit()
                current_value = f'{currval:,.4f}'
                edit_box.setText(current_value)
                edit_box.editingFinished.connect(self.save_changes)
                edit_box.setAlignment(Qt.AlignRight)

                layout.addWidget(edit_box)
                self.edit_boxes[param] = edit_box

        self.setLayout(layout)
    def objchanged(self):
        selected_option = self.combo_box.currentText()
        self.parent_obj.objecttype = selected_option
        

    def save_changes(self): 
        for param, edit_box in self.edit_boxes.items():
            if hasattr(self.parent_obj, param):
                new_value = float(edit_box.text())
                setattr(self.parent_obj, param, new_value)

class ExtraFigureWindow(QMainWindow):
    def __init__(self, parent_obj):
        super().__init__()

        self.parent_obj = parent_obj
        self.init_ui()
        
    def check_panzoom_mode(self):
        if (self.toolbar.mode == 'zoom rect') or (self.toolbar.mode == 'pan/zoom'):
            return 1
        else:
            return 0
        

    def init_ui(self):
        
        self._main = QWidget()
        self.setWindowTitle("Plot")
        self.setCentralWidget(self._main)
        layout = QVBoxLayout(self._main)

        self.static_canvas = FigureCanvas(Figure(figsize=(5, 4)))
        self.toolbar = NavigationToolbar2QT(self.static_canvas, self)
        
        # Create a horizontal layout for toolbar and options
        toolbar_row = QHBoxLayout()
        toolbar_row.addWidget(self.toolbar)
        
        # Plot Options
        gb = QGroupBox("Options")
        gb_layout = QHBoxLayout()
        gb_layout.setContentsMargins(2, 2, 2, 2) # Compact margins
        self.logx_cb = QCheckBox("Log X")
        self.logx_cb.toggled.connect(self.update_scales)
        gb_layout.addWidget(self.logx_cb)
        
        self.logy_cb = QCheckBox("Log Y")
        self.logy_cb.toggled.connect(self.update_scales)
        gb_layout.addWidget(self.logy_cb)
        
        gb.setLayout(gb_layout)
        # Add groupbox to the horizontal row
        toolbar_row.addWidget(gb)
        
        # Add the row layout to the main layout
        layout.addLayout(toolbar_row)

        layout.addWidget(self.static_canvas)
        self.PowerLabel = QLabel(' ')
        layout.addWidget(self.PowerLabel)


        self.ax = self.static_canvas.figure.subplots()

        self.static_canvas.figure.subplots_adjust(right=0.95,top=0.95)
        self.static_canvas.figure.set_facecolor('white')
        self.static_canvas.figure.set_edgecolor('black')
        self.ax.grid(True, which='both',ls='--',lw=0.5)
        self.ax.set_facecolor('white')
        self.ax.spines['bottom'].set_color('black')
        self.ax.spines['top'].set_color('black')
        self.ax.spines['left'].set_color('black')
        self.ax.spines['right'].set_color('black')
        self.ax.tick_params(axis='both', colors='black')
        self.ax.xaxis.label.set_color('black')
        self.ax.yaxis.label.set_color('black')      
        # self.ax.set_aspect(self.parent_obj.aspectratio)
        self.static_canvas.draw()
        

        

        self.zoom_start_pixel = None
        self.zoom_initial_xlim = None
        self.zoom_initial_ylim = None
        self.zoom_pivot = None # Initialize zoom_pivot

        self.highlight = None  #   highlighting rectangle
        self.start_point = None  # start point for  x-range selection
        self.xdata = []
        self.ydata = []
        def on_click(event):
            """Called when the mouse button is pressed."""
            if not self.check_panzoom_mode():
                
                if event.inaxes != self.ax:
                    return

                # Left Click: Region Selection
                if event.button == 1:
                    self.start_point = event.xdata   
                
                    # remove existing highlight  
                    if self.highlight:
                        self.highlight.remove()
                    self.highlight = self.ax.axvspan(self.start_point, self.start_point, color='yellow', alpha=0.3)
                    self.static_canvas.draw()
                
                # Right Click: Zoom Init
                elif event.button == 3:
                    self.zoom_start_pixel = (event.x, event.y)
                    self.zoom_initial_xlim = self.ax.get_xlim()
                    self.zoom_initial_ylim = self.ax.get_ylim()
                    self.zoom_pivot = (event.xdata, event.ydata) # Store pivot in data coordinates

        def on_drag(event):
            if not self.check_panzoom_mode():
                """Called when the mouse is dragged."""
    
                if event.inaxes != self.ax:
                    return
                
                # Left Drag: Update Region Selection
                if event.button == 1 and self.start_point is not None:
                    self.highlight.remove()
                    self.highlight = self.ax.axvspan(self.start_point, event.xdata, color='yellow', alpha=0.3)
                    self.static_canvas.draw()
                
                # Right Drag: Zoom
                elif event.button == 3 and self.zoom_start_pixel is not None:
                    # Calculate drag distance in pixels
                    dx_pix = event.x - self.zoom_start_pixel[0]
                    dy_pix = event.y - self.zoom_start_pixel[1]
                    
                    # Define sensitivity
                    # Using exponential scaling for smoother zoom
                    sx = 0.99 ** dx_pix # Drag Right (Positive dx) -> Zoom In (Smaller Range)
                    sy = 0.99 ** dy_pix # Drag Up (Positive dy) -> Zoom In (Smaller Range)
                    
                    # Get Initial Layout from when the right-click started
                    x0, x1 = self.zoom_initial_xlim
                    y0, y1 = self.zoom_initial_ylim
                    
                    if self.zoom_pivot and self.zoom_pivot[0] is not None and self.zoom_pivot[1] is not None:
                        px, py = self.zoom_pivot
                        
                        # Calculate new limits maintaining the pivot point's relative position
                        new_x0 = px - (px - x0) * sx
                        new_x1 = px + (x1 - px) * sx
                        new_y0 = py - (py - y0) * sy
                        new_y1 = py + (y1 - py) * sy
                        
                        self.ax.set_xlim(new_x0, new_x1)
                        self.ax.set_ylim(new_y0, new_y1)
                        self.static_canvas.draw()


        def on_release(event):
            """Called when the mouse button is released."""
            if not self.check_panzoom_mode():
                
                # Left Click Release: Finalize Selection
                if event.button == 1:
                    if self.start_point is None or event.inaxes != self.ax:
                        return
                    
                    x_min = min(self.start_point, event.xdata)
                    x_max = max(self.start_point, event.xdata)
                    
                    if self.xdata:
                        print(f'From {x_min:0.3f} eV to {x_max:0.3f} eV')
    
                        listofpowers = ''
                        for xvec,yvec in zip(self.xdata,self.ydata):
                        
                            mask = (xvec >= x_min) & (xvec <= x_max)         
                            x = xvec[mask]
                            y = yvec[mask]
                            power = trapezoid(y, x)
                            powertot = trapezoid(yvec, xvec)
                            
                            photons_pereV_persec = y/x/(1.602e-19)
                            photons_perbunch24bunch = trapezoid(photons_pereV_persec,x)/(271647*24)
                            photons_perbunch48bunch = trapezoid(photons_pereV_persec,x)/(271647*48)
                            
                            Eperbunch_24_uJ = power/(271647*24)*1e6
                            Eperbunch_48_uJ = power/(271647*48)*1e6
                            # (70)/(24000)/(6.52e6)/(1.602e-19)
                            FractionOfTotal = power/powertot
                            
                            
                            listofpowers=listofpowers + f'{power:0.3f} W, ph/153ns: {photons_perbunch24bunch:0.2e}({Eperbunch_24_uJ:.2f} µJ), ph/77ns: {photons_perbunch48bunch:0.2e}({Eperbunch_48_uJ:.2f} µJ), Fraction of Total: {FractionOfTotal:0.2f}'
                        self.PowerLabel.setText(listofpowers)
                        print(listofpowers) 
                
                    # reset for next selection
                    self.start_point = None
                
                # Right Click Release
                elif event.button == 3:
                     self.zoom_start_pixel = None
                     self.zoom_initial_xlim = None
                     self.zoom_initial_ylim = None
                     self.zoom_pivot = None

        def on_scroll(event):
            if not self.check_panzoom_mode() and event.inaxes == self.ax:
                 # Zoom Factor
                base_scale = 1.2
                if event.button == 'up':
                    scale_factor = 1/base_scale
                else: # 'down'
                    scale_factor = base_scale
                
                # Get current range
                x0, x1 = self.ax.get_xlim()
                y0, y1 = self.ax.get_ylim()
                
                # Pivot is current mouse position
                px, py = event.xdata, event.ydata
                
                # Calculate new range centered on pivot
                new_x0 = px - (px - x0) * scale_factor
                new_x1 = px + (x1 - px) * scale_factor
                new_y0 = py - (py - y0) * scale_factor
                new_y1 = py + (y1 - py) * scale_factor
                
                self.ax.set_xlim(new_x0, new_x1)
                self.ax.set_ylim(new_y0, new_y1)
                self.static_canvas.draw()

        self.static_canvas.mpl_connect('button_press_event', on_click)
        self.static_canvas.mpl_connect('motion_notify_event', on_drag)
        self.static_canvas.mpl_connect('button_release_event', on_release)
        self.static_canvas.mpl_connect('scroll_event', on_scroll)
        
        plt.show()

    def update_scales(self):
        if self.logx_cb.isChecked():
            self.ax.set_xscale('log')
        else:
            self.ax.set_xscale('linear')
        
        if self.logy_cb.isChecked():
            self.ax.set_yscale('log')
        else:
            self.ax.set_yscale('linear')
        self.static_canvas.draw()

class UndulatorCalcConfig(QObject):
    def __init__(self,parent=None,Undulator=None,period=None,harmonicnumber=None,energy=None,npts=None,emin=None,emax=None,K=None,ringcurrent=None,APSU=1):
        super().__init__(parent)
        self.undulatorvalue = Undulator
        self.period = period
        self.harmonicnumber = harmonicnumber
        self.energy = energy
        self.npts = npts
        self.emin = emin
        self.emax = emax
        self.K  = K
        self.ringcurrent = ringcurrent
        self.APSU = APSU
        
        
        if self.APSU==1:
            self.sourcesizermsH=12.29e-6
            self.sourcesizermsV=8.26e-6
            self.sourcedivergenceH = 2.336e-6
            self.sourcedivergenceV = 3.474e-6
            self.ringenergy = 6.0 #GeV
            if period==0.023:
                self.numperiods = 206
            elif period==0.014:
                self.numperiods = 338
            else: 
                undulatorlength = 4.738 #meters
                self.numperiods = undulatorlength/period        
        else:
            self.sourcesizermsH=280.3e-6
            self.sourcesizermsV=13.9e-6
            self.sourcedivergenceH = 11.7e-6
            self.sourcedivergenceV = 3.6e-6
            self.ringenergy = 7.0 #GeV
            if period==0.027:
                self.numperiods = 88
            elif period==0.017:
                self.numperiods = 137
            else: 
                undulatorlength = 0.027*88 #meters
                self.numperiods = undulatorlength/period    
                
    def print(self):
        attrdict = self.__dict__
        print('-------Undulator Configuration-------')
        for k in attrdict.keys():
            print(f'\t{k}: {attrdict[k]}')
        print('-------------------------------------')
        
        

        
    
# %%  2D Undulator Calc


#### Start of 2D Undulator Calc
 
def process_chunk(chunk_args):
    """
    Process a chunk of pixels
    
    Parameters:
        chunk_args:  tuple containing:
          - x_range: the range of x-indices to process
          - y_range: the range of y-indices (usually all rows)
          - xvals: 1D array of x positions
          - yvals: 1D array of y positions
          - UndulatorObj: configuration object (e.g., UndulatorCalcConfig)
    
    Returns:  A list of tuples (x, y, specpower0) for each pixel.
    """
    x_range, y_range, xvals, yvals, UndulatorObj = chunk_args
    results = []
    dx = xvals[1] - xvals[0]
    dy = yvals[1] - yvals[0]
    for x in x_range:
        for y in y_range:
            if UndulatorObj.undulatorvalue == 0:
                continue  # Skip if undulator is off
            edat, specpower0 = undulatorcalc(UndulatorObj, xvals[x], yvals[y], dx, dy)
            results.append((x, y, specpower0))
    return results



class WorkerSignals(QObject):
    progress = pyqtSignal(int, int, np.ndarray)  # x, y coordinates and pixel value
    finished = pyqtSignal()

class PixelCalculationTask(QRunnable):
    def __init__(self, xvals,yvals, x_range, y_range, signals, UndulatorObject):
        super().__init__()
        self.x_range = x_range  # range to calculate
        self.y_range = y_range  
        self.signals = signals  
        self.xvals = xvals
        self.yvals = yvals
        self.UndulatorObj = UndulatorObject
    
    def run(self):
        
        dx = self.xvals[1]-self.xvals[0]
        dy = self.yvals[1]-self.yvals[0]
        for x in self.x_range:
            for y in self.y_range:
  

                
                if self.UndulatorObj.undulatorvalue==0:
                    return

                edat,specpower0 = undulatorcalc(self.UndulatorObj,self.xvals[x],self.yvals[y],dx,dy)

                self.signals.progress.emit(x, y, specpower0) #progress
        self.signals.finished.emit()
                  
    
                
            
        

def xoppy_calc_undulator_spectrum_streamline(ELECTRONENERGY=6.04,ELECTRONENERGYSPREAD=0.001,ELECTRONCURRENT=0.2,\
                              ELECTRONBEAMSIZEH=0.000395,ELECTRONBEAMSIZEV=9.9e-06,\
                              ELECTRONBEAMDIVERGENCEH=1.05e-05,ELECTRONBEAMDIVERGENCEV=3.9e-06,\
                              PERIODID=0.018,NPERIODS=222,KV=1.68,KH=0.0,KPHASE=0.0,DISTANCE=30.0,
                              GAPH=0.001,GAPV=0.001,GAPH_CENTER=0.0,GAPV_CENTER=0.0,\
                              PHOTONENERGYMIN=3000.0,PHOTONENERGYMAX=55000.0,PHOTONENERGYPOINTS=500,
                              USEEMITTANCES=1):
    METHOD=2
    

    bl = OrderedDict()
    bl['ElectronBeamDivergenceH'] = ELECTRONBEAMDIVERGENCEH
    bl['ElectronBeamDivergenceV'] = ELECTRONBEAMDIVERGENCEV
    bl['ElectronBeamSizeH'] = ELECTRONBEAMSIZEH
    bl['ElectronBeamSizeV'] = ELECTRONBEAMSIZEV
    bl['ElectronCurrent'] = ELECTRONCURRENT
    bl['ElectronEnergy'] = ELECTRONENERGY
    bl['ElectronEnergySpread'] = ELECTRONENERGYSPREAD
    bl['Kv'] = KV
    bl['Kh'] = KH
    bl['Kphase'] = KPHASE
    bl['NPeriods'] = NPERIODS
    bl['PeriodID'] = PERIODID
    bl['distance'] = DISTANCE
    bl['gapH'] = GAPH
    bl['gapV'] = GAPV
    bl['gapHcenter'] = GAPH_CENTER
    bl['gapVcenter'] = GAPV_CENTER

    if USEEMITTANCES:
        zero_emittance = False
    else:
        zero_emittance = True

    outFile = None

    codata_mee = codata.m_e * codata.c**2 / codata.e # electron mass in eV
    gamma = bl['ElectronEnergy'] * 1e9 / codata_mee
    
    m2ev = codata.c * codata.h / codata.e      
    resonance_wavelength = (1 + (bl['Kv']**2 + bl['Kh']**2) / 2.0) / 2 / gamma**2 * bl["PeriodID"]
    resonance_energy = m2ev / resonance_wavelength



    ptot = (NPERIODS/6) * codata.value('characteristic impedance of vacuum') * \
           ELECTRONCURRENT * codata.e * 2 * np.pi * codata.c * gamma**2 * (KV**2+KH**2) / PERIODID
    
    if METHOD == 2:
        # get the maximum harmonic number
        h_max = int(2.5*PHOTONENERGYMAX/resonance_energy)
        e, f = calc1d_srw_removeoutput(bl,photonEnergyMin=PHOTONENERGYMIN,photonEnergyMax=PHOTONENERGYMAX,
              photonEnergyPoints=PHOTONENERGYPOINTS,fileName=outFile,fileAppend=False,zero_emittance=zero_emittance,
              srw_max_harmonic_number=h_max)


    # if zero_emittance:
        # print("\nNo emittance calculation")

    if METHOD == 1 and len(e) == 0: raise Exception("Invalid Input Parameters")

    # power_in_spectrum = f.sum()*1e3*codata.e*(e[1]-e[0])
    # print("\nPower from integral of spectrum: %8.3f W"%(power_in_spectrum))
    # print("\nRatio Power from integral of spectrum over Total emitted power: %5.4f"%(power_in_spectrum / ptot))

    spectral_power = f * codata.e * 1e3
    # try:
    #     cumulated_power = spectral_power.cumsum() * np.abs(e[0] - e[1])
    # except:
    #     cumulated_power = 0.0

    return e, f, spectral_power
    
def undulatorcalc(UndulatorObject,x,y,dx,dy):
    UndulatorObject.print()
    # print(UndulatorObject)
    energy, flux, spectral_power = xoppy_calc_undulator_spectrum_streamline(
    ELECTRONENERGY=UndulatorObject.ringenergy,
    ELECTRONENERGYSPREAD=0.00098,
    ELECTRONCURRENT=UndulatorObject.ringcurrent/1000,
    ELECTRONBEAMSIZEH=UndulatorObject.sourcesizermsH,# ELECTRONBEAMSIZEH=1.29e-05,
    ELECTRONBEAMSIZEV=UndulatorObject.sourcesizermsV,#ELECTRONBEAMSIZEV=2.5e-06,
    ELECTRONBEAMDIVERGENCEH=UndulatorObject.sourcedivergenceH, #        ELECTRONBEAMDIVERGENCEH=8.7e-06,
    ELECTRONBEAMDIVERGENCEV=UndulatorObject.sourcedivergenceV,   #     ELECTRONBEAMDIVERGENCEV=3.6e-06,
    PERIODID=UndulatorObject.period,
    NPERIODS=UndulatorObject.numperiods,
    KV=UndulatorObject.K,
    KH=0.0,
    KPHASE=0.0,
    DISTANCE=31.0,
    GAPH=dx*1e-3,
    GAPV=dy*1e-3,
    GAPH_CENTER=x*1e-3,
    GAPV_CENTER=y*1e-3,
    PHOTONENERGYMIN=UndulatorObject.emin,
    PHOTONENERGYMAX=UndulatorObject.emax,
    PHOTONENERGYPOINTS=UndulatorObject.npts,
    USEEMITTANCES=1);
    return energy,spectral_power


    
def calc1d_srw_removeoutput(bl,photonEnergyMin=3000.0,photonEnergyMax=55000.0,photonEnergyPoints=500,zero_emittance=False,
              srw_max_harmonic_number=None,fileName=None,fileAppend=False):
    codata_mee = np.array(codata.physical_constants["electron mass energy equivalent in MeV"][0])
    m2ev = codata.c * codata.h / codata.e      # lambda(m)  = m2eV / energy(eV)

    r"""
        run SRW for calculating flux

        input: a dictionary with beamline
        output: file name with results
    """

    global scanCounter
 

    cte = codata.e/(2*np.pi*codata.electron_mass*codata.c)
    B0 = bl['Kv']/bl['PeriodID']/cte

    try:
        B0x = bl['Kh']/bl['PeriodID']/cte
    except:
        B0x = 0.0

    try:
        Kphase = bl['Kphase']
    except:
        Kphase = 0.0

    if srw_max_harmonic_number == None:
        gamma = bl['ElectronEnergy'] / (codata_mee * 1e-3)

        try:
            Kh = bl['Kh']
        except:
            Kh = 0.0

        resonance_wavelength = (1 + (bl['Kv']**2 + Kh**2) / 2.0) / 2 / gamma**2 * bl["PeriodID"]
        resonance_energy = m2ev / resonance_wavelength

        srw_max_harmonic_number = int(photonEnergyMax / resonance_energy * 2.5)
        print ("Max harmonic considered:%d ; Resonance energy: %g eV\n"%(srw_max_harmonic_number,resonance_energy))


    Nmax = srw_max_harmonic_number # 21,61
 
    if B0x == 0:    #*********** Conventional Undulator
        harmB = srwlib.SRWLMagFldH() #magnetic field harmonic
        harmB.n = 1 #harmonic number ??? Mostly asymmetry
        harmB.h_or_v = 'v' #magnetic field plane: horzontal ('h') or vertical ('v')
        harmB.B = B0 #magnetic field amplitude [T]
        und = srwlib.SRWLMagFldU([harmB])
        und.per = bl['PeriodID'] #period length [m]
        print( und.per)
        
        und.nPer = bl['NPeriods'] #number of periods (will be rounded to integer)
        print( und.nPer)
        #Container of all magnetic field elements
        magFldCnt = srwlib.SRWLMagFldC([und], srwlib.array('d', [0]), srwlib.array('d', [0]), srwlib.array('d', [0]))
    else:  #***********Undulator (elliptical)
        magnetic_fields = []
        magnetic_fields.append(srwlib.SRWLMagFldH(1, 'v',
                                           _B=B0,
                                           _ph=0.0,
                                           _s=1, # 1=symmetrical, -1=antisymmetrical
                                           _a=1.0))
        magnetic_fields.append(srwlib.SRWLMagFldH(1, 'h',
                                           _B=B0x,
                                           _ph=Kphase,
                                           _s=1,
                                           _a=1.0))
        und = srwlib.SRWLMagFldU(_arHarm=magnetic_fields, _per=bl['PeriodID'], _nPer=bl['NPeriods'])

        magFldCnt = srwlib.SRWLMagFldC(_arMagFld=[und],
                                        _arXc=srwlib.array('d', [0.0]),
                                        _arYc=srwlib.array('d', [0.0]),
                                        _arZc=srwlib.array('d', [0.0]))


    #***********Electron Beam
    eBeam = srwlib.SRWLPartBeam()
    eBeam.Iavg = bl['ElectronCurrent'] #average current [A]
    eBeam.partStatMom1.x = 0. #initial transverse positions [m]
    eBeam.partStatMom1.y = 0.
    # eBeam.partStatMom1.z = 0 #initial longitudinal positions (set in the middle of undulator)
    eBeam.partStatMom1.z = - bl['PeriodID']*(bl['NPeriods']+4)/2 # initial longitudinal positions
    eBeam.partStatMom1.xp = 0 #initial relative transverse velocities
    eBeam.partStatMom1.yp = 0
    eBeam.partStatMom1.gamma = bl['ElectronEnergy']*1e3/codata_mee #relative energy

    if zero_emittance:
        sigX     = 1e-25
        sigXp    = 1e-25
        sigY     = 1e-25
        sigYp    = 1e-25
        sigEperE = 1e-25
    else:
        sigX =  bl['ElectronBeamSizeH'] #horizontal RMS size of e-beam [m]
        sigXp = bl['ElectronBeamDivergenceH'] #horizontal RMS angular divergence [rad]
        sigY =  bl['ElectronBeamSizeV'] #vertical RMS size of e-beam [m]
        sigYp = bl['ElectronBeamDivergenceV'] #vertical RMS angular divergence [rad]
        sigEperE = bl['ElectronEnergySpread']

    # print("calc1d_srw: starting calculation using ElectronEnergySpead=%e \n"%((sigEperE)))

    #2nd order stat. moments:
    eBeam.arStatMom2[0] = sigX*sigX #<(x-<x>)^2>
    eBeam.arStatMom2[1] = 0 #<(x-<x>)(x'-<x'>)>
    eBeam.arStatMom2[2] = sigXp*sigXp #<(x'-<x'>)^2>
    eBeam.arStatMom2[3] = sigY*sigY #<(y-<y>)^2>
    eBeam.arStatMom2[4] = 0 #<(y-<y>)(y'-<y'>)>
    eBeam.arStatMom2[5] = sigYp*sigYp #<(y'-<y'>)^2>
    eBeam.arStatMom2[10] = sigEperE*sigEperE #<(E-<E>)^2>/<E>^2

    #***********Precision Parameters
    arPrecF = [0]*5 #for spectral flux vs photon energy
    arPrecF[0] = 1 #initial UR harmonic to take into account
    arPrecF[1] = Nmax #final UR harmonic to take into account
    arPrecF[2] = 1.5 #longitudinal integration precision parameter
    arPrecF[3] = 1.5 #azimuthal integration precision parameter
    arPrecF[4] = 1 #calculate flux (1) or flux per unit surface (2)

    #***********UR Stokes Parameters (mesh) for Spectral Flux
    stkF = srwlib.SRWLStokes() #for spectral flux vs photon energy
    #srio stkF.allocate(10000, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
    stkF.allocate(photonEnergyPoints, 1, 1) #numbers of points vs photon energy, horizontal and vertical positions
    stkF.mesh.zStart = bl['distance'] #longitudinal position [m] at which UR has to be calculated
    stkF.mesh.eStart = photonEnergyMin #initial photon energy [eV]
    stkF.mesh.eFin =   photonEnergyMax #final photon energy [eV]
    stkF.mesh.xStart = bl['gapHcenter'] - bl['gapH']/2 #initial horizontal position [m]
    stkF.mesh.xFin =   bl['gapHcenter'] + bl['gapH']/2 #final horizontal position [m]
    stkF.mesh.yStart = bl['gapVcenter'] - bl['gapV']/2 #initial vertical position [m]
    stkF.mesh.yFin =   bl['gapVcenter'] + bl['gapV']/2 #final vertical position [m]

    #**********************Calculation (SRWLIB function calls)
    # print('Performing Spectral Flux (Stokes parameters) calculation ... ') # , end='')
     
    srwlib.srwl.CalcStokesUR(stkF, eBeam, und, arPrecF)

    # print('Done calc1d_srw calculation in %10.3f s'%(time.time()-t0))
    #**********************Saving results

    # if fileName is not None:
    #     if fileAppend:
    #         f = open(fileName,"a")
    #     else:
    #         scanCounter = 0
    #         f = open(fileName,"w")
    #         f.write("#F "+fileName+"\n")

    #     f.write("\n")
    #     scanCounter +=1
    #     f.write("#S %d Undulator spectrum calculation using SRW\n"%(scanCounter))

    #     for i,j in bl.items(): # write bl values
    #         f.write ("#UD %s = %s\n" % (i,j) )
    #     f.write("#UD photonEnergyMin =  %f\n"%(photonEnergyMin))
    #     f.write("#UD photonEnergyMax =  %f\n"%(photonEnergyMax))
    #     f.write("#UD photonEnergyPoints =  %d\n"%(photonEnergyPoints))
    #     f.write("#UD B0 =  %f\n"%(B0))

    #     #
    #     # write flux to file
    #     #
    #     header="#N 4 \n#L PhotonEnergy[eV]  PhotonWavelength[A]  Flux[phot/sec/0.1%bw]  Spectral Power[W/eV]\n"
    #     f.write(header)

    eArray = np.zeros(photonEnergyPoints)
    intensArray = np.zeros(photonEnergyPoints)
    for i in range(stkF.mesh.ne):
        ener = stkF.mesh.eStart+i*(stkF.mesh.eFin-stkF.mesh.eStart)/np.array((stkF.mesh.ne-1)).clip(min=1)
        # if fileName is not None: f.write(' ' + repr(ener) + '   ' + repr(m2ev/ener*1e10) + '    ' +
        #         repr(stkF.arS[i]) + '    ' +
        #         repr(stkF.arS[i]*codata.e*1e3) + '\n')
        eArray[i] = ener
        intensArray[i] = stkF.arS[i]

    # if fileName is not None:
    #     f.close()

    #     if fileAppend:
    #         print("Data appended to file: %s"%(os.path.join(os.getcwd(),fileName)))
    #     else:
    #         print("File written to disk: %s"%(os.path.join(os.getcwd(),fileName)))


    return (eArray,intensArray)



class ArrayViewer(QWidget):
    def __init__(self,xyspectrum=None, edat=None,xvals=None,yvals=None):
        if xyspectrum is None:
            print('Insufficient data passed')
            self.data_loaded = False
            return None
        else:
            self.xyspectrum = xyspectrum
            self.edat = edat
            self.xvals = xvals
            self.yvals = yvals
            self.data_loaded = True
            self.index = 0
            
        self.cmin,self.cmax = np.percentile(xyspectrum,[1,99])

        super().__init__()
        self.setWindowTitle("Undulator, 2D Energy Slice Viewer")
        



        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)

        # self.load_button = QPushButton("Load .npz File")
        # self.load_button.clicked.connect(self.load_file)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self.on_slider_changed)

        self.energy_input = QLineEdit()
        self.energy_input.setPlaceholderText("Enter energy value")
        self.energy_input.editingFinished.connect(self.on_energy_entered)

        layout = QVBoxLayout()
        # layout.addWidget(self.load_button)
        layout.addWidget(self.canvas)
        layout.addWidget(self.slider)
        layout.addWidget(self.energy_input)

        self.setLayout(layout)
        
        if self.data_loaded:
            self.slider.setMinimum(0)
            self.slider.setMaximum(self.xyspectrum.shape[2] - 1)
            self.slider.setValue(self.index)
            self.slider.setEnabled(True)
        self.plot_slice()

    # def load_file(self):
    #     file_path, _ = QFileDialog.getOpenFileName(self, "Open .npz file", "", "NumPy files (*.npz)")
    #     if not file_path:
    #         return

    #     try:
    #         data = np.load(file_path)
    #         self.xyspectrum = data["array_3d"]
    #         self.edat = data["array_2d_1"].flatten()
    #         self.xvals = data["array_2d_2"].flatten()
    #         self.yvals = data["array_2d_3"].flatten()

    #         self.index = 0
    #         self.slider.setMinimum(0)
    #         self.slider.setMaximum(self.xyspectrum.shape[2] - 1)
    #         self.slider.setValue(self.index)
    #         self.slider.setEnabled(True)

    #         self.data_loaded = True
    #         self.plot_slice()

        # except Exception as e:
        #     QMessageBox.critical(self, "Error", f"Failed to load file:\n{str(e)}")

    def plot_slice(self):
        if not self.data_loaded:
            return

        if self.index < 0 or self.index >= self.xyspectrum.shape[2]:
            return

        self.fig.clear()
        ax = self.fig.add_subplot(111)
        im = ax.imshow(self.xyspectrum[:, :, self.index],
                       extent=[self.xvals[0], self.xvals[-1], self.yvals[-1], self.yvals[0]],
                       aspect='auto')
        self.fig.colorbar(im, ax=ax)
        ax.set_title(f"Energy: {self.edat[self.index]:.3f} eV")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        self.canvas.draw()

    def on_slider_changed(self, value):
        if self.data_loaded:
            self.index = value
            self.plot_slice()

    def on_energy_entered(self):
        if not self.data_loaded:
            return

        try:
            energy_val = float(self.energy_input.text())
            closest_index = np.abs(self.edat - energy_val).argmin()
            self.slider.setValue(closest_index)  # This will trigger on_slider_changed
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid number for energy.")
            
class PlotterWidget(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        
        self.centralWidg = QWidget()
        self.setWindowTitle("Plot")
        self.setCentralWidget(self.centralWidg)
         
        main_layout = QVBoxLayout(self)
        self.centralWidg.setLayout(main_layout)
        
        menu_bar = QMenuBar(self)
        file_menu = menu_bar.addMenu('File')
        
        Viewer_menu = menu_bar.addMenu('E-Slice Viewer')
        OpenViewer_action = QAction('Open', self)        
        OpenViewer_action.triggered.connect(self.openEsliceViewer)
        Viewer_menu.addAction(OpenViewer_action)

        
        save_action = QAction('Save 2D Spectrum (.npz)', self)        
        save_action.triggered.connect(self.save_2Dspectrum_toNPZ)
        file_menu.addAction(save_action)
        
        load_action = QAction('Load 2D Spectrum (.npz)', self)        
        load_action.triggered.connect(self.Load2Dspectrum)
        file_menu.addAction(load_action)

        main_layout.setMenuBar(menu_bar)
        
        
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax1 = self.figure.add_subplot(121)  
        self.ax2 = self.figure.add_subplot(122) 
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        main_layout.addWidget(self.toolbar)
        main_layout.addWidget(self.canvas)
        self.parentwidget = parent
        
        params_groupbox = QGroupBox("Simulation Parameters")
        params_layout = QVBoxLayout()
        
        x_layout = QHBoxLayout()
        x_layout.addWidget(QLabel("X Start:"))
        slitsize_h = parent.WBSsizeH
        slitsize_v= parent.WBSsizeV
        # slitcenter_h = parent.WBScenterH
        # slitcenter_v = parent.WBScenterV
        self.xstart_edit = QLineEdit(f'{-1*slitsize_h/2:3f}')
        x_layout.addWidget(self.xstart_edit)
        
        x_layout.addWidget(QLabel("X End:"))
        self.xend_edit = QLineEdit(f'{slitsize_h/2:3f}')
        x_layout.addWidget(self.xend_edit)
        
        y_layout = QHBoxLayout()
        y_layout.addWidget(QLabel("Y Start:"))
        self.ystart_edit = QLineEdit(f'{-1*slitsize_v/2:3f}')
        y_layout.addWidget(self.ystart_edit)
        
        y_layout.addWidget(QLabel("Y End:"))
        self.yend_edit = QLineEdit(f'{slitsize_v/2:3f}')
        y_layout.addWidget(self.yend_edit)
        
        numpts_layout = QHBoxLayout()
        numpts_layout.addWidget(QLabel("Number of Points:"))
        self.numpts_edit = QLineEdit("10")
        numpts_layout.addWidget(self.numpts_edit)
        
        params_layout.addLayout(x_layout)
        params_layout.addLayout(y_layout)
        params_layout.addLayout(numpts_layout)
        
        slidergroupbox = QGroupBox()
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
        sliderlayout = QVBoxLayout()
        params_layout.addWidget(slidergroupbox)
        slidergroupbox.setLayout(sliderlayout)
        
        sliderlayout.addWidget(QLabel("ColorScale:"))
        sliderlayout.addWidget(self.vmin_slider)
        sliderlayout.addWidget(self.vmax_slider)
        
        self.plot_button = QPushButton("Calculate and Plot")
        self.plot_button.clicked.connect(self.start_calculation)
        params_layout.addWidget(self.plot_button)
        
        self.progress_bar = QProgressBar()
        params_layout.addWidget(self.progress_bar)
        
        params_groupbox.setLayout(params_layout)
        main_layout.addWidget(params_groupbox)

        
        self.x_values = None
        self.y_values = None
        self.current_step = 0
        self.total_steps = 100

        self.thread_pool = QThreadPool()
        self.thread_pool.setMaxThreadCount(6)  # # of threads in the pool
        
        self.total_pixels = 0
        self.calculated_pixels = 0
        self.finished_workers = 0

        
        self.canvas.resizeEvent = self.on_resize  #  resize handler
        
        self.image_data = None
        self.image_display = None
        self.total_pixels = 0
        self.calculated_pixels = 0
        self.edat = None
        self.xyspectrum = None
        self.xvals=None
        self.yvals = None
        self.sumspectrum = None
        
        self.canvas.mpl_connect("button_press_event", self.on_image_click)
    def save_2Dspectrum_toNPZ(self):
        save_arrays(self, self.xyspectrum, self.edat, self.xvals, self.yvals)
    def openEsliceViewer(self):
        self.esliceviewer = ArrayViewer(self.xyspectrum,self.edat,self.xvals,self.yvals)
        self.esliceviewer.show()
        
    def Load2Dspectrum(self):
        self.xyspectrum, self.edat, self.xvals, self.yvals = load_arrays(self)
        
        xmin = np.min(self.xvals)
        print(xmin)
        xmax = np.max(self.xvals)
        
        ymin = np.min(self.yvals)
        ymax = np.max(self.yvals)
        
        if self.xyspectrum is not None:
            print('Loaded')
            a,b,c = np.shape(self.xyspectrum)
            self.image_data = np.zeros((a,b))
            for x,y in itertools.product(range(a),range(b)):
                self.image_data[x,y] = float(a*b)*trapezoid(self.xyspectrum[x,y,:],self.edat)
                
            self.ax1.clear()
            self.ax2.clear()
            self.ax1.set_title("Total Intensity",fontsize=8)
            self.ax2.set_title("Local Spectrum",fontsize=8)

            self.image_display = self.ax1.imshow(self.image_data, extent=[xmin,xmax,ymin,ymax],cmap="gray", vmin=0, vmax=255)
            # self.ax1.axis("off")
            self.figure.tight_layout()
            self.canvas.draw()
            self.sumspectrum = np.sum(self.xyspectrum,(0,1))
            self.update_image_display()
                    
                    
        
    def on_image_click(self, event):

        if event.inaxes != self.ax1:
            return

        x =event.xdata
        
        y = event.ydata
        
        xindex = np.argmin(abs(self.xvals-x))
        yindex = np.argmin(abs(self.yvals-y))

        spectrum = self.xyspectrum[yindex, xindex,:]
        spectrum = spectrum/np.max(spectrum)
        

        self.ax2.clear()
        self.ax2.set_title(f"Normalized Spectrum at ({x:.3f} mm, {y:.3f} mm)",fontsize=8)
        
        self.ax2.plot(self.edat, spectrum)
        self.ax2.plot(self.edat, self.sumspectrum/np.max(self.sumspectrum))
        
        self.canvas.draw_idle()
        
    def update_image_display(self):
        if self.image_data is not None:
            maxval = np.max(self.image_data)
            self.image_display.set_clim(vmin=maxval/100*np.min([self.vmin_slider.value(),self.vmax_slider.value()]), vmax=maxval/100*np.max([self.vmin_slider.value(),self.vmax_slider.value()]))
            self.canvas.draw()        
            
    def on_resize(self, event):
        """
        ensure the axes scale proportionally
        """
        self.figure.tight_layout()  
        self.canvas.draw_idle()  
        super(FigureCanvas, self.canvas).resizeEvent(event)
    
    def start_calculation(self):
        
        if self.parentwidget.APSUcheck.checkState()>0:
            APSU=1
            offset = 1.00375
        else:
            APSU=0
            offset = 1
            

        self.starttime = time.time()
            
        #  values from input fields
        try:
            xstart = float(self.xstart_edit.text())
            xend = float(self.xend_edit.text())
            ystart = float(self.ystart_edit.text())
            yend = float(self.yend_edit.text())
            numpts = int(self.numpts_edit.text())
        except ValueError:
            print("Invalid input; please enter numeric values.")
            return
        
        undulatorvalue,period,harmonicnumber,energy,npts,emin,emax = self.parentwidget.getcurrentscreenvalues()
        if undulatorvalue==0:
            return
        ringcurrent = float(self.parentwidget.ringcurrent_edit.text())
        if self.parentwidget.useKcheckbox.checkState()==2:
            K = float(self.parentwidget.KVal_edit.text())
        else:
            offsetenergy = energy*offset
            K = self.parentwidget.KfromRingEnergyAndUndulatorPeriod(self.parentwidget.ringenergy,period,offsetenergy,harmonicnumber)
            
        UndulatorConfig = UndulatorCalcConfig(self,undulatorvalue,period,harmonicnumber,energy,npts,emin,emax,K,ringcurrent,APSU)

        #just get the e domain once, so it doesn't have to be passed a bunch
        self.edat,_ = undulatorcalc(UndulatorConfig,0,0,.001,.001)
        self.xyspectrum = np.zeros((numpts,numpts,len(self.edat)))
        
        

        self.image_data = np.zeros((numpts, numpts), dtype=np.uint16)
        

        self.ax1.clear()
        self.ax2.clear()
        self.ax1.set_title("Total Intensity")
        self.ax2.set_title("Local Spectrum")


        self.image_display = self.ax1.imshow(self.image_data, cmap="gray", vmin=0, vmax=255)
        self.ax1.axis("off")
        self.figure.tight_layout()
        
        self.total_pixels = numpts * numpts
        self.calculated_pixels = 0
        self.progress_bar.setMaximum(self.total_pixels)
        self.progress_bar.setValue(0)

        self.xvals = np.linspace(xstart,xend,numpts)
        self.yvals = np.linspace(ystart,yend,numpts)
        
        # divide image into chunks and start workers
        if self.thread_pool.maxThreadCount()>numpts:
            numthreads = numpts
        else:
            numthreads = self.thread_pool.maxThreadCount()

        self.finishedworkers=0
        chunks = np.append(np.int16(np.linspace(0,numpts-1,numthreads)),numpts)
        
        
        for i in range(numthreads):
            
            x_range = range(chunks[i],chunks[i+1])
            print(*x_range)
            y_range = range(numpts)
            
            signals = WorkerSignals()
            signals.progress.connect(self.update_image_pixel)
            signals.finished.connect(self.finish_calculation)
            
            task = PixelCalculationTask(self.xvals, self.yvals, x_range,y_range,signals, UndulatorConfig)
            self.thread_pool.start(task)

    
    def update_image_pixel(self, x, y, specpowerreturned):
        tote = trapezoid(specpowerreturned,self.edat)
        self.xyspectrum[y,x,:]=specpowerreturned
        # Update the pixel in the image data
        self.image_data[y, x] = tote*self.total_pixels
        # print(value)
        self.image_display.set_data(self.image_data)   
        
        # Update progress and display
        self.calculated_pixels += 1
        # print(self.calculated_pixels)
        self.progress_bar.setValue(self.calculated_pixels)
        
        # if np.mod(self.calculated_pixels,20)==0:
        #     self.canvas.draw()


            
    def finish_calculation(self):
        # Check if all pixels are calculated
        self.finishedworkers +=1
        self.sumspectrum = np.sum(self.xyspectrum,(0,1))
        print('Workers finished: '+str(self.finishedworkers))
        endtime = time.time()
        print(f'Duration: {endtime-self.starttime:3f}')
        if self.calculated_pixels >= self.total_pixels:
            print("Image generation finished.")
        self.canvas.draw()
        self.update_image_display()
        
def save_arrays(parent, array_3d, array_2d_1, array_2d_2, array_2d_3):
    """
    Save a 3D numpy array along with 3 2D arrays to a file selected by the user.


    """
    # Check dimensions
    if array_3d.ndim != 3 or array_2d_1.ndim != 1 or array_2d_2.ndim != 1 or array_2d_3.ndim != 1:
        QMessageBox.critical(parent, "Error", "Ensure the arrays are 3D (array_3d) and 1D (array_2d_1, array_2d_2, array_2d_3).")
        return

    # Open a file dialog to select the save location
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getSaveFileName(
        parent,
        "Save Arrays",
        "",
        "NumPy Archive (*.npz);;All Files (*)",
        options=options
    )

    if file_path:
        try:
            # Save the arrays to a .npz file
            np.savez(file_path, array_3d=array_3d, array_2d_1=array_2d_1, array_2d_2=array_2d_2, array_2d_3=array_2d_3)
            QMessageBox.information(parent, "Success", f"Arrays saved successfully to {file_path}.")
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to save arrays: {e}")


def load_arrays(parent):
    """
    Load a 3D numpy array and two 2D arrays from a file

    :param parent: The parent widget for the QFileDialog.
    :return: A tuple (array_3d, array_2d_1, array_2d_2, array_2d_3) or None if an error occurred.
    """

    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(
        parent,
        "Load Arrays",
        "",
        "NumPy Archive (*.npz);;All Files (*)",
        options=options
    )

    if file_path:
        try:
            # Load the arrays from the .npz file
            data = np.load(file_path)
            array_3d = data["array_3d"]
            array_2d_1 = data["array_2d_1"]
            array_2d_2 = data["array_2d_2"]
            array_2d_3 = data["array_2d_3"]

            # Ensure dimensions are correct
            if array_3d.ndim != 3 or array_2d_1.ndim != 1 or array_2d_2.ndim != 1:
                QMessageBox.critical(parent, "Error", "The loaded arrays do not have the correct dimensions.")
                return None, None, None, None

            QMessageBox.information(parent, "Success", f"Arrays loaded successfully from {file_path}.")
            return array_3d, array_2d_1, array_2d_2, array_2d_3
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to load arrays: {e}")
            return None, None, None, None
    return None, None, None, None

        
def save_numpy_array(parent, array):
    """
    Save a 3D numpy array to a file
    """
    
    if array.ndim != 3:# check if 3D
        QMessageBox.critical(parent, "Error", "The array must be 3D.")
        return
 
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getSaveFileName(
        parent,
        "Save Numpy Array",
        "",
        "NumPy Files (*.npy);;All Files (*)",
        options=options
    )

    if file_path:
        try: 
            np.save(file_path, array)
            QMessageBox.information(parent, "Success", f"Array saved successfully to {file_path}.")
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to save array: {e}")

def load_numpy_array(parent):
    """
    Load a 3D numpy array from a file 

    :param parent: The parent widget for the QFileDialog.
    :return: The loaded 3D numpy array or None if an error occurred.
    """
    # Open a file dialog to select the file to load
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(
        parent,
        "Load Numpy Array",
        "",
        "NumPy Files (*.npy);;All Files (*)",
        options=options
    )

    if file_path:
        try:
            # Load the numpy array from the selected file
            array = np.load(file_path)
            if array.ndim != 3:
                QMessageBox.critical(parent, "Error", "The loaded array is not 3D.")
                return None
            QMessageBox.information(parent, "Success", f"Array loaded successfully from {file_path}.")
            return array
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to load array: {e}")
            return None
    return None

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DCSSpectrum()
    window.resize(600, 200)



    window.show()
    sys.exit(app.exec_())
