import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel,QComboBox, QFormLayout)
from PyQt5.QtCore import Qt
import h5py
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt
import os
ML_dir_path = os.path.dirname(os.path.realpath(__file__))


class MultilayerConfig(QWidget):
    def __init__(self,parent_obj,startenergy=19000):
        
        super().__init__()
        self.parent_obj = parent_obj
        self.setWindowTitle("Multilayer Config")
        self.layout = QVBoxLayout()
        
        self.formlayout = QFormLayout()
        
        self.enableMLcheck = QCheckBox()
        # self.checkboxlayout.addWidget(self.enableMLcheck)
        self.formlayout.addRow(QLabel('Enable ML'), self.enableMLcheck)
        self.enableMLcheck.stateChanged.connect(self.setMLenable)

        self.layout.addLayout(self.formlayout)
        self.combo = QComboBox()
        self.combo.addItems(["3.9nm #706", "2.2nm #692"])
        self.combo.currentTextChanged.connect(self.MLchanged)
        
        self.formlayout.addRow(QLabel('Multilayer'), self.combo)
        self.setLayout(self.layout)
        
        self.MLenergyedit = QLineEdit()
        self.formlayout.addRow(QLabel('Desired Energy (eV)'), self.MLenergyedit)
        self.MLenergyedit.setText(str(startenergy))
        
        
        self.setanglebutton = QPushButton('Set Angle')
        self.plotanglebutton = QPushButton('Plot vs Angle')
        self.formlayout.addRow(self.setanglebutton,self.plotanglebutton)
        self.setanglebutton.clicked.connect(self.setangle)
        self.plotanglebutton.clicked.connect(self.plotangle)
        
        self.angleedit = QLineEdit()
        self.formlayout.addRow(QLabel('Incidence Angle (deg)'),self.angleedit)
        self.angleedit.setPlaceholderText('Incidence angle (deg)')
        
        
        
        self.MLReflectivityInterpolator = []
               
    def MLchanged(self):
            self.MLReflectivityInterpolator=[]
        
    def setMLenable(self,state):
        if state == 0:
            self.parent_obj.multilayerenabled=0
        elif state == 2:
            self.parent_obj.multilayerenabled=1
            
    def initialguessMLAngle(self,energyeV,dspacing_nm):
        return 180/np.pi*np.arcsin(1239.84198/(energyeV)/2/(dspacing_nm)) 
        
    def readMLData(self):
        
        MLname = self.combo.currentText()
        print(MLname)
        
        if MLname == "3.9nm #706": #Add more cases later
            # self.multilayerdatafile_path ='I:\Sinclair\Oasys\multilayer706_BIG_5keVto80keV_theta0p01To1deg.h5'
            
            self.multilayerdatafile_path  = os.path.join(ML_dir_path,'multilayer706_BIG_5keVto80keV_theta0p01To3deg_1.h5')#='I:\Sinclair\Oasys\multilayer706_BIG_5keVto80keV_theta0p01To3deg_1.h5'
        elif MLname == "2.2nm #692": #Add more cases later
            self.multilayerdatafile_path =os.path.join(ML_dir_path,'multilayer692_BIG_5keVto80keV_theta0p01To3deg_1.h5')
        else:
            print('Invalid ML specified')
            return
            
            
        with h5py.File(self.multilayerdatafile_path, 'r') as h5pyfile:
            MLEdat = h5pyfile['/MLayer/EnergyAngleScan/axis_x'][:]
            MLangledat = h5pyfile['/MLayer/EnergyAngleScan/axis_y'][:]
            MLimgdat = h5pyfile['/MLayer/EnergyAngleScan/image_data'][:]
        
        

        
        self.MLReflectivityInterpolator = RegularGridInterpolator((MLangledat,MLEdat), MLimgdat,
                                 bounds_error=False, fill_value=None)    
        self.datastored = MLname
    
    def setangle(self):        
        MLname = self.combo.currentText()
        # print(MLname)
        if MLname == "3.9nm #706": #Add more cases later
            dspacing = 3.93 #nm
        elif MLname== "2.2nm #692":
            dspacing = 2.2
        desiredEnergy = float(self.MLenergyedit.text())
        
        
        if not self.MLReflectivityInterpolator:
            self.readMLData()
            
        guessangle = self.initialguessMLAngle(desiredEnergy,dspacing)         
        print(guessangle)
        anglemat = np.linspace(0.9*guessangle,1.1*guessangle,200)            
        reflvalues = self.MLReflectivityInterpolator((anglemat,desiredEnergy))
        # plt.plot(anglemat,reflvalues)
        
        argm = np.argmax(reflvalues);
        besttransangle = anglemat[argm]
        self.angleedit.setText(f'{besttransangle:.4f}')
        # plt.plot(self.besttransangle,reflvalues[argm],marker='x')
        # plt.show()
        
    def plotangle(self):        
        MLname = self.combo.currentText()
        # print(MLname)
        if MLname == "3.9nm #706": #Add more cases later
            dspacing = 3.93 #nm
        elif MLname== "2.2nm #692":
            dspacing = 2.2
        desiredEnergy = float(self.MLenergyedit.text())        
        
        if not self.MLReflectivityInterpolator:
            self.readMLData()
            
        guessangle = self.initialguessMLAngle(desiredEnergy,dspacing)         
        print(guessangle)
        anglemat = np.linspace(0.9*guessangle,1.1*guessangle,200)            
        reflvalues = self.MLReflectivityInterpolator((anglemat,desiredEnergy))
        plt.plot(anglemat,reflvalues)
        
        argm = np.argmax(reflvalues);
        besttransangle = anglemat[argm]
        plt.plot(besttransangle,reflvalues[argm],marker='x')
        plt.show() 
        
    



if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MultilayerConfig(app)
    window.resize(200, 200)


    window.show()
    sys.exit(app.exec_())
