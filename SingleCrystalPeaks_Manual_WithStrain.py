import numpy as np
from ase.io import read
from ase.build import make_supercell
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget,
    QLineEdit, QLabel, QHBoxLayout, QDoubleSpinBox,
    QGridLayout, QTableWidget, QTableWidgetItem,
    QGroupBox,QPushButton,QRadioButton,QButtonGroup
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QColor,QCloseEvent
from scipy.spatial.transform import Rotation as R
import spglib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D  # necessary for 3D plots
from xoppylib.scattering_functions.xoppy_calc_f0 import xoppy_calc_f0 
import xraylib
from xoppylib.scattering_functions.f1f2_calc import f1f2_calc
from scipy.interpolate import interp1d

np.seterr(all='raise')

# Constants
hc = 12398.4198  # electron volts * angstrom
r0 = 2.81794032E-5  # angstroms
JouleIneV = 6.242e+18  # eV per Joule

def setFusionpalette(app):
    app.setStyle("Fusion")
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
    
# --- Generic Reciprocal Lattice Calculation ---
def compute_reciprocal(cell):
    cell_mat = np.array(cell)
    a = cell_mat[0]
    b = cell_mat[1]
    c = cell_mat[2]
    V = np.dot(a, np.cross(b, c))
    a_star = np.cross(b, c) / V
    b_star = np.cross(c, a) / V
    c_star = np.cross(a, b) / V
    return np.array([a_star, b_star, c_star])

def generate_q0(h, k, l, structure, wavelength, U):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    hkl = np.array([h, k, l])
    try:
        d_spacing = 1 / np.sqrt(hkl @ np.linalg.inv(G) @ hkl.T)
    except:
        return 0
    reciprocal_lattice = compute_reciprocal(structure.cell)
    hkl_vector = np.dot([h, k, l], reciprocal_lattice)
    hkl_vector = np.dot(U, hkl_vector)
    q_magnitude = 2 * np.pi / d_spacing
    q0 = hkl_vector / np.linalg.norm(hkl_vector) * q_magnitude
    return q0

def sinAngleBetweenTwoVectors(a, b):
    a = np.array(a)
    b = np.array(b)
    val = np.linalg.norm(np.cross(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    if 1.01 > val > 1:
        val = 1
    return val

def cosAngleBetweenTwoVectors(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

def compute_d_spacing(structure, h, k, l):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    hkl = np.array([h, k, l])
    try:
        d_spacing = 1 / np.sqrt(hkl @ np.linalg.inv(G) @ hkl.T)
    except np.linalg.LinAlgError:
        return 0
    return d_spacing

def compute_listof_d_spacings(structure, hkllist):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    dspacinglist = []
    for el in hkllist:
        hkl = np.array(el)
        try:
            d_spacing = 1 / np.sqrt(hkl @ np.linalg.inv(G) @ hkl.T)
            dspacinglist.append(d_spacing)
        except np.linalg.LinAlgError:
            dspacinglist.append(0)
    return dspacinglist

def compute_listof_d_spacings_and_angles(structure, hkllist, U, max2th, maxE=np.inf, minE=0):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    reciprocal_lattice = compute_reciprocal(structure.cell)
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    goodhkllist = []
    dspacinglist = []
    twothlist = []
    E_bragg_list = []
    qveclist = []
    kslist = []
    F_hkl_list = []
    kveclist = []
    elements = set(atom.symbol for atom in structure)
    f0_values = precompute_f0(xraylib, elements, x_min=0.0, x_max=8.0, num_points=100)
    
    for el in hkllist:
        h, k, l = el
        hkl = np.array(el)
        try:
            d_spacing = 1 / np.sqrt(hkl @ np.linalg.inv(G) @ hkl.T)
        except np.linalg.LinAlgError:
            d_spacing = 0
        
        qdir = generate_q_unit_vector(h, k, l, reciprocal_lattice, U)
        if d_spacing != 0 and qdir[2] != 0:
            qnorm = 2 * np.pi / d_spacing
            qvec = qnorm * qdir
            qx, qy, qz = qvec
            k0 = -(qnorm**2) / (2 * qz)
            lambda_bragg = 2 * np.pi / k0 if k0 != 0 else 0
            kvec = k0 * np.array([0, 0, 1])
            ks = kvec + qvec
            
            if lambda_bragg != 0 and k0 > 0:
                try:
                    E_bragg = hc / lambda_bragg
                    sintwoth = sinAngleBetweenTwoVectors(ks, kvec)
                    costwoth = cosAngleBetweenTwoVectors(ks, kvec)
                    twoth = np.arcsin(sintwoth)
                    if np.degrees(twoth) < max2th and costwoth > 0:
                        if minE < E_bragg < maxE:
                            F_hkl = calculate_structure_factor(structure, h, k, l, E_bragg, f0_values)
                            if F_hkl > 0.01:
                                E_bragg_list.append(E_bragg)
                                F_hkl_list.append(F_hkl)
                                goodhkllist.append(hkl)
                                twothlist.append(twoth)
                                dspacinglist.append(d_spacing)
                                qveclist.append(qvec)
                                kslist.append(ks)
                                kveclist.append(kvec)
                except Exception as e:
                    print(e, f': sintwoth = {sintwoth}')
                    pass
    return dspacinglist, twothlist, E_bragg_list, kveclist, kslist, qveclist, goodhkllist, F_hkl_list

def gen_shifted_U(U,N = 300,deltaangle=1):
    phi_list = np.linspace(-180,180,num=N,endpoint=False)
    U_list = []
    for phi in phi_list:
        U_list.append(R.from_euler('ZXZ', (phi,deltaangle,-phi),degrees=True).as_matrix() @ U)
    return U_list
        
        

def VisualizeSingle(structure, hkl, U):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    reciprocal_lattice = compute_reciprocal(structure.cell)
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    G = np.array([
        [a**2, a * b * np.cos(gamma), a * c * np.cos(beta)],
        [b * a * np.cos(gamma), b**2, b * c * np.cos(alpha)],
        [c * a * np.cos(beta), c * b * np.cos(alpha), c**2]
    ])
    h, k, l = hkl
    hkl_arr = np.array(hkl)
    try:
        d_spacing = 1 / np.sqrt(hkl_arr @ np.linalg.inv(G) @ hkl_arr.T)
    except np.linalg.LinAlgError:
        d_spacing = 0
    s = generate_q_unit_vector(h, k, l, reciprocal_lattice, U)
    beamaxis = [0, 0, 1]
    sintwoth = sinAngleBetweenTwoVectors(s, beamaxis)
    twoth = np.arcsin(sintwoth)
    sinth = np.sin(twoth / 2)
    lambda_bragg = 2 * d_spacing * sinth
    if lambda_bragg != 0:
        try:
            E_bragg_list = hc / lambda_bragg
        except:
            E_bragg_list = None

def generate_unique_hkl(structure, max_index):
    cell = (
        structure.cell,
        structure.get_scaled_positions(),
        structure.get_atomic_numbers(),
    )
    symmetry_dataset = spglib.get_symmetry_dataset(cell)
    rotations = symmetry_dataset.rotations
    identityrotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], type(rotations[0][0][0]))
    rotations = [r for r in rotations if not np.array_equal(r, identityrotation)]
    rotations.append(identityrotation)
    hkl_list = [(h, k, l) for h in range(-max_index, max_index + 1)
                for k in range(-max_index, max_index + 1)
                for l in range(-max_index, max_index + 1) if h != 0 or k != 0 or l != 0]
    unique_hkl_list = []
    for ij in range(len(hkl_list)):
        hkl = hkl_list[ij]
        if hkl != (None, None, None):
            transformed_hkls = {tuple(np.dot(rotation, hkl)) for rotation in rotations}
            matches = transformed_hkls.intersection(hkl_list)
            hkl_list = [(None, None, None) if t in transformed_hkls else t for t in hkl_list]
            hkl_list[ij] = hkl
            multiplicity = len(matches)
            refl = {'hkl': hkl, 'multiplicity': multiplicity}
            unique_hkl_list.append(refl)
    return unique_hkl_list
def solveBragg_single(hkl_vector_xtalframe,U,d_spacing,detdist,emin=5000,emax=100000):
 
    hkl_lab = np.dot(U, hkl_vector_xtalframe)
    q_magnitude = 2 * np.pi / d_spacing
    q_vec = hkl_lab / np.linalg.norm(hkl_lab) * q_magnitude
  
    k0 = -1* q_magnitude**2 / (2 * q_vec[2])
    lam = 2 * np.pi / k0
    E_bragg = 12398.4198 / lam
    if emin<E_bragg<emax:
        ks = np.array((0,0,k0)) + q_vec
        kdir = ks/k0
        peakvector = kdir / kdir[2] * detdist 
        dx = peakvector[0]
        dy = peakvector[1]
        return (dx, dy)
    else:
        return None
    
def generate_q_unit_vector(h, k, l, reciprocal_lattice, U):
    hkl_vector = np.dot([h, k, l], reciprocal_lattice)
    hkl_vector = np.dot(U, hkl_vector)
    s_unitvec = hkl_vector / np.linalg.norm(hkl_vector)
    return s_unitvec



def calculate_structure_factor(structure, h, k, l, energy, f0_values):
    F_hkl = 0.0 + 0.0j
    d_spacing = compute_d_spacing(structure, h, k, l)
    if d_spacing == 0:
        return 0
    x = 1 / (2 * d_spacing)
    for atom in structure:
        try:
            x_pos, y_pos, z_pos = atom.position / structure.get_cell().lengths()
            element = atom.symbol
            f0 = f0_values[element](x)
            atomic_number = atom.number
            f1, f2 = f1f2_calc(element, energy, theta=None, F=0, density=None, rough=None,
                                 verbose=False, material_constants_library=xraylib)
            f1_corrected = f1 - atomic_number
            f = f0 + f1_corrected + 1j * f2
            phase = np.exp(2j * np.pi * (h * x_pos + k * y_pos + l * z_pos))
            F_hkl += f * phase
        except Exception as e:
            print(e, f'Energy: {energy}')
            pass
            
    return np.float64(abs(np.squeeze(F_hkl)))

def precompute_f0(material_constants_library, elements, x_min, x_max, num_points):
    f0_values = {}
    for element in elements:
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
        x_grid = f0_data['data'][0, :]
        f0_grid = f0_data['data'][1, :]
        f0_values[element] = interp1d(x_grid, f0_grid, kind='cubic', fill_value="extrapolate")
    return f0_values
 
class StrainWindow(QMainWindow):
    def __init__(self, strain_callback, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Edit Strain Tensor")
        self.strain_callback = strain_callback
        
        self.table = QTableWidget(3, 3)
        self.table.setFixedSize(300, 150)
        for i in range(3):
            for j in range(3):
                item = QTableWidgetItem("0")
                item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(i, j, item)
                
        self.table.cellChanged.connect(self.onCellChanged)
        
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.addWidget(self.table)
        self.setCentralWidget(widget)
    
    def onCellChanged(self, row, column):
        strain = np.zeros((3,3))
        try:
            for i in range(3):
                for j in range(3):
                    item = self.table.item(i, j)
                    strain[i, j] = float(item.text())
            self.strain_callback(strain)
        except Exception as e:
            print("Error parsing strain value:", e)
 
class CrystalReflectionVisualizer(QMainWindow):
    def __init__(self, structure, max_index=5, xtal_to_det_distance_meters=0.1,
                 minimumEnergy=5000, maximumEnergy=70000, phi=0,
                 max2th=50, theta=0, psi=0):
        super().__init__()
        self.target_mount_axis = 90-28 #degrees, for laser shock
        self.currview = None
        self.qvariationplot =None
        self.doplotvariation=False
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.max_index = max_index
        self.structure = structure
        self.max2th = max2th
        self.F_hkl_list =None
        self.minimumEnergy = minimumEnergy
        self.maximumEnergy = maximumEnergy
        self.base_hkl_list = [(h, k, l) for h in range(-self.max_index, self.max_index + 1)
                              for k in range(-self.max_index, self.max_index + 1)
                              for l in range(-self.max_index, self.max_index + 1)
                              if h != 0 or k != 0 or l != 0]
        self.hkl=None
        self.axlim = 80
        self.detdistance_mm = xtal_to_det_distance_meters * 1000
        self.showgeo = False

        self.setWindowTitle("3D Crystal Reflection Visualizer")
        self.setGeometry(50, 50, 800, 1000)
        
        # Store the original cell to allow reapplication of strain
        self.original_cell = self.structure.cell.copy()
        self.current_strain = np.zeros((3,3))
        
        # Add menu bar for strain editing
        menubar = self.menuBar()
        strain_menu = menubar.addMenu("Strain")
        strain_action = strain_menu.addAction("Edit Strain Tensor")
        strain_action.triggered.connect(self.openStrainWindow)
        
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        layoutGrid = QGridLayout(self.central_widget)
                
        self.detector_fig = Figure(figsize=(5, 4))
        self.canvas_det = FigureCanvas(self.detector_fig)
        self.canvas_det.mpl_connect('pick_event', self.on_text_pick)
        self.detector_ax = self.canvas_det.figure.add_subplot(111)
        
        layoutfig = QVBoxLayout()
        layoutGrid.addLayout(layoutfig,  0, 0, 1, 2)
        layoutfig.addWidget(self.canvas_det)
        self.toolbar = NavigationToolbar2QT(self.canvas_det, self)
        layoutfig.addWidget(self.toolbar)
        
        
        # --- Diffraction Geometry view (3D) ---
        self.geobox = QWidget()
        
        self.geobox.destroyed.connect(self.showgeo_clicked)
        self.geobox.closeEvent=self.showgeoclosed
        self.layoutgeo = QHBoxLayout()
        self.layoutgeo_1 = QVBoxLayout()
        self.layoutgeo.addLayout(self.layoutgeo_1)
        
        self.geobox.setLayout(self.layoutgeo)
        self.canvas_geo = FigureCanvas(Figure(figsize=(5, 4)))
        self.ax_geo = self.canvas_geo.figure.add_subplot(111, projection='3d')
        
        self.layoutgeo_1.addWidget(self.canvas_geo)
        
        # layoutGrid.addWidget(self.geobox, 0, 2, 1, 6)
        self.axlabels=[]
        
        for i in range(5):
            self.axlabels.append(QLabel())
            self.layoutgeo_1.addWidget(self.axlabels[i])
        
        # --- Crystal Structure view (3D) ---
        self.canvas_crystal = FigureCanvas(Figure(figsize=(5, 4)))
        self.ax_crystal = self.canvas_crystal.figure.add_subplot(111, projection='3d')
        self.layoutgeo.addWidget(self.canvas_crystal)
        
        # --- Angle and parameter editors ---
        sideLayout = QVBoxLayout()
        radiolayout = QHBoxLayout()
        sideLayout.addLayout(radiolayout)
        self.b1 = QRadioButton("Gun")
        self.b1.setChecked(True)
        self.b1.toggled.connect(lambda:self.driverselected(self.b1))
        radiolayout.addWidget(self.b1)
        
        self.b2 = QRadioButton("Laser")
        self.b2.toggled.connect(lambda:self.driverselected(self.b2))
        self.group1 = QButtonGroup()
        self.group1.addButton(self.b1)
        self.group1.addButton(self.b2)
        
        radiolayout.addWidget(self.b2)
        
        
        layoutGrid.addLayout(sideLayout, 1, 0, 1, 2)
        
        def makeangleedit(label):
            hbox = QHBoxLayout()
            sideLayout.addLayout(hbox)
            hbox.addWidget(QLabel(label))
            spinner = QDoubleSpinBox()
            spinner.setMinimum(-180)
            spinner.setMaximum(180)
            spinner.setValue(0.0)
            spinner.valueChanged.connect(self.update)
            hbox.addWidget(spinner)
            return spinner

        self.phi_spinner = makeangleedit('Phi')
        self.theta_spinner = makeangleedit('Theta')
        self.psi_spinner = makeangleedit('Psi')

        def addgeneraledit(label, startvalstr):
            hbox = QHBoxLayout()
            sideLayout.addLayout(hbox)
            hbox.addWidget(QLabel(label))
            edit = QLineEdit(startvalstr)
            edit.editingFinished.connect(self.update)
            hbox.addWidget(edit)
            return edit
        
        self.showgeo_button = QPushButton('Show Geometry')
        sideLayout.addWidget(self.showgeo_button)
        self.showgeo_button.clicked.connect(self.showgeo_clicked)

        self.detdistedit = addgeneraledit('Det Dist (mm)', '100')
        self.minimumEnergyedit = addgeneraledit('minE (eV)', f'{self.minimumEnergy}')
        self.maximumEnergyedit = addgeneraledit('maxE (eV)', f'{self.maximumEnergy}')
        self.axrangeedit = addgeneraledit('Axis Range (mm)', f'{self.axlim}')
        self.TextSize_edit = addgeneraledit('TextSize', '12')
        self.update()
        startindex = 0
        Estart = self.energies[startindex]
        self.wavelength = hc / Estart * 1e-10  # 1 Angstrom
        self.hkl = self.hkl_list[0]
        self.visualize_geometry_and_crystal()
        
    def driverselected(self,b):
    
        if b.text() == "Gun":
            if b.isChecked() == True:
                self.target_mount_axis = 90-28 #degrees, for gun
 
				
        if b.text() == "Laser":
            if b.isChecked() == True:
                self.target_mount_axis = 90+38 #degrees, for laser shock
          
        self.update()
        
    def genRotation(self,phi, theta, psi):
        # rotation = R.from_euler('YZY', [phi, theta, psi], degrees=True)
        firstrotation = R.from_euler('YZY', [self.target_mount_axis, 0, 0], degrees=True)
        rotation = R.from_euler('ZXZ', [phi, theta, psi], degrees=True)
        return firstrotation.as_matrix() @ rotation.as_matrix()
        
    def showgeoclosed(self,event:QCloseEvent):
        self.showgeo_clicked()
        event.accept()
        
    def showgeo_clicked(self):
        if self.showgeo:
            self.showgeo = False
            self.showgeo_button.setText('Show Geometry')
            
            self.geobox.hide()
        else:
            self.showgeo = True
            self.showgeo_button.setText('Hide Geometry')
            self.geobox.show()

    def update(self):
        try:
            self.detdistance_mm = float(self.detdistedit.text())
            self.minimumEnergy = float(self.minimumEnergyedit.text())
            self.maximumEnergy = float(self.maximumEnergyedit.text())
            self.axlim = float(self.axrangeedit.text())
        except:
            pass
        self.doplotvariation=False
        self.qvariationplot = None
        self.phi = self.phi_spinner.value()
        self.theta = self.theta_spinner.value()
        self.psi = self.psi_spinner.value()
        self.U = self.genRotation(self.phi, self.theta, self.psi)

        res = compute_listof_d_spacings_and_angles(self.structure, self.base_hkl_list, self.U, self.max2th,
                                                     minE=self.minimumEnergy, maxE=self.maximumEnergy)
        (dspacinglist, twothlist, E_bragg_list, k0, k_prime_list,
         q_hkl_list, goodhkllist, F_hkl_list) = res

        self.k0 = k0
        self.hkl_list = goodhkllist
        self.q_hkl_list = q_hkl_list
        self.k_prime_list = k_prime_list
        self.energies = E_bragg_list
        self.F_hkl_list = F_hkl_list
        if len(E_bragg_list) > 0:
            self.energy_colors = (np.array(E_bragg_list) - min(E_bragg_list)) / (max(E_bragg_list) - min(E_bragg_list))
            msizeindicator = F_hkl_list#np.log10(F_hkl_list)
            
            self.markersizes =(np.array(msizeindicator) - min(msizeindicator)) / (max(msizeindicator) - min(msizeindicator) + 1e-9) 
            self.draw_det_reflections()
        else:
            print('No reflections')
            self.detector_ax.clear()
            self.canvas_det.draw()
        
        # refresh the 3D views after updating
        if self.showgeo:
            self.visualize_geometry_and_crystal()
            
    def on_text_pick(self, event):
        # when pick event occurs on the detector canvas
        artist = event.artist
        if hasattr(artist, '_hkl') and hasattr(artist, '_energy'):
            hkl = artist._hkl
            energy = artist._energy
            q = artist._q
            F = artist._Fhkl
            self.reflection_clicked(hkl, energy,q,F)

    def reflection_clicked(self, hkl, energy,q,F):
        self.hkl = hkl
        self.energyplotted = energy
        self.wavelength = hc / energy * 1e-10  
        self.qhkl = q
        self.Fhkl = F
        print(f"Reflection: {hkl}, Energy: {energy:.1f}, F: {F:.3f}")
        self.doplotvariation=True
        self.draw_det_reflections()
        self.visualize_geometry_and_crystal()
        
        
        # self.visualize_crystal_structure()
    # def draw_det_reflections(self):
    #     self.detector_ax.clear()
    #     self.detector_ax.set_xlabel('Detector Position (mm)')
    #     self.detector_ax.set_ylabel('Detector Position (mm)')
    #     self.detector_ax.set_xlim(-self.axlim, self.axlim)
    #     self.detector_ax.set_ylim(-self.axlim, self.axlim)
    #     self.detector_ax.invert_yaxis()
    #     colormap = cm.plasma
    #     for k0, hkl, q, k_prime, norm_color in zip(
    #         self.k0, self.hkl_list, self.q_hkl_list, self.k_prime_list, self.energy_colors
    #     ):
    #         ks_x, ks_y, ks_z = k_prime
    #         det_x = ks_x * self.detdistance_mm / ks_z
    #         det_y = ks_y * self.detdistance_mm / ks_z
    #         color = colormap(norm_color)
    #         self.detector_ax.scatter(det_x, det_y, color=color, s=25)
    #         self.detector_ax.text(det_x+2, det_y, str(hkl), color=color, fontsize=6)
    #     self.canvas_det.draw()
    
    def draw_det_reflections(self):
        # Clear the 2D detector axes and set labels/limits.
        self.detector_fig.clf()
        ax = self.detector_fig.add_subplot(111)
        ax.set_xlabel("Detector Position (mm)")
        ax.set_ylabel("Detector Position (mm)")
        ax.set_xlim(-self.axlim, self.axlim)
        ax.set_ylim(-self.axlim, self.axlim)
        ax.invert_yaxis()
        ax.plot(0,0,'x',color='black',markersize=12)
        reciprocal_lattice = compute_reciprocal(self.structure.cell)
        if self.doplotvariation:
            U_list = gen_shifted_U(self.U,N = 100,deltaangle=1) 
            self.detpoints_variations = []
        # Draw each reflection marker and add an interactive text label.
        
        for energy, hkl, q, k_prime, norm_color,markerscale,F in zip(
                self.energies, self.hkl_list, self.q_hkl_list, self.k_prime_list, self.energy_colors,self.markersizes,self.F_hkl_list):
            ks_x, ks_y, ks_z = k_prime
            det_x = ks_x * self.detdistance_mm / ks_z
            det_y = ks_y * self.detdistance_mm / ks_z
            color = cm.plasma(norm_color)
            # Draw the marker.
            try: 
                markersize_edit = float(self.TextSize_edit.text())
            except:
                markersize_edit = 12
                
            if markersize_edit==0:
                markersize = 6+3*markerscale
            else:
                markersize = markersize_edit
                
            ax.plot(det_x, det_y, 'o', color=color, markersize=markersize)
            # Add the text label with picker enabled.
            text_obj = ax.text(det_x, det_y, str(hkl), color=color, picker=True, fontsize = markersize)
            # Store custom attributes on the text object.
            text_obj._hkl = hkl
            text_obj._energy = energy
            text_obj._q = q
            text_obj._Fhkl = F
            h,k,l = hkl
            d_spacing= compute_d_spacing(self.structure, h, k, l)
            if self.doplotvariation:
                emin = float(self.minimumEnergyedit.text())
                emax = float(self.maximumEnergyedit.text())
          
                hkl_vector_xtalframe = np.dot([h,k,l], reciprocal_lattice)
                for Ushift in U_list:
                    self.detpoints_variations.append(solveBragg_single(hkl_vector_xtalframe,Ushift,d_spacing,self.detdistance_mm,emin=emin,emax=emax))
                new_detpoints_variations_list = [x for x in self.detpoints_variations if x is not None]
                
        if self.doplotvariation:
            x_val_variations, y_val_variations = zip(*new_detpoints_variations_list)
            if self.qvariationplot is not None:
                try:
                    self.qvariationplot.remove()
                except NotImplementedError: 
                    pass
            
            lines=ax.plot(x_val_variations,y_val_variations,linestyle='none', marker='o')
            self.qvariationplot = lines[-1]
                
        self.canvas_det.draw()

    def visualize_geometry_and_crystal(self):
        if self.currview is None:
            self.currview = (0, 90) #elev, azim
        else:
            self.currview = (self.ax_geo.elev, self.ax_geo.azim)
            
        evelationviewangle, azimuthalviewangle = self.currview
        
        
        if self.hkl is None:
            # If there is no reflection to use, clear the 3D axes and return.
            self.ax_geo.clear()
            self.canvas_geo.draw()
            self.ax_crystal.clear()
            self.canvas_crystal.draw()
            return
    
        # --- Diffraction Geometry View ---
        self.ax_geo.clear()
        q0 = generate_q0(self.hkl[0], self.hkl[1], self.hkl[2], self.structure, self.wavelength, self.U)
        q0_dir = q0 / np.linalg.norm(q0)
        k_in_magnitude = 2 * np.pi / self.wavelength
        k_in = np.array([0, 0, k_in_magnitude])
        k_out = k_in + q0
        k_in_dir = k_in / np.linalg.norm(k_in)
        k_out_dir = k_out / np.linalg.norm(k_out)
        
        # self.ax_geo.quiver(0, 0, 0, k_in_dir[0], k_in_dir[1], k_in_dir[2], color='b',
        #                    arrow_length_ratio=0.1, label='k_in')
        # self.ax_geo.quiver(0, 0, 0, k_out_dir[0], k_out_dir[1], k_out_dir[2], color='r',
        #                    arrow_length_ratio=0.1, label='k_out')
        # self.ax_geo.quiver(0, 0, 0, q0_dir[0], q0_dir[1], q0_dir[2], color='g',
        #                    arrow_length_ratio=0.1, label='q0')
        gunaxis = self.genRotation(0,0,0) @ np.array((0,0,1))
        crystal_100 = self.U @ np.array((1,0,0))
        crystal_010 = self.U @ np.array((0,1,0))
        crystal_001 = self.U @ np.array((0,0,1))
        origin = np.array([0, 0, 0])
        def addarrow(num,startvec,vector,color,label):
            self.axlabels[num].setText((f'<font color="{color}">{label}</font>'))
            
            self.ax_geo.quiver(*startvec, *vector,   color=color, label=label)
        
        addarrow(0,origin,k_in_dir,'blue','k_in')
        addarrow(1,origin,gunaxis,'red','gun axis')
        addarrow(2,origin,crystal_100,'yellow','crystal_100')
        addarrow(3,origin,crystal_010,'cyan','crystal_010')
        addarrow(4,origin,crystal_001,'green','crystal_001')
        
        
        
        
        
        # Plot arrows using quiver.
        # print(f"deltak/k: {(k_out_mag-k_in_magnitude)/k_in_magnitude:.4f}")
        # kscaling = 1/k_in_magnitude*self.detdistance_mm / 1000.0
        # self.ax_geo.quiver(*origin, *k_in*kscaling,   color='blue', label='k_in')
        # self.ax_geo.quiver(*origin, *k_out*kscaling,  color='red', label='k_out')
        # self.ax_geo.quiver(*k_in * kscaling, *q0 * kscaling,  color='green', label='q0')
        # self.ax_geo.set_box_aspect([0.2, 0.2, 0.2])
        
        
        # self.ax_geo.quiver(*origin, *k_in_dir,   color='blue', label='k_in')
        # self.ax_geo.quiver(*origin, *gunaxis,   color='red', label='gun axis')
        # self.ax_geo.quiver(*origin, *crystal_100,   color='yellow', label='crystal 100')
        # self.ax_geo.quiver(*origin, *crystal_010,   color='cyan', label='crystal 010')
        # self.ax_geo.quiver(*origin, *crystal_001,   color='green', label='crystal 001')        
        # self.ax_geo.legend()
        self.ax_geo.set_xlabel('X')
        self.ax_geo.set_ylabel('Y')
        self.ax_geo.set_zlabel('Z')
        z_plane = self.detdistance_mm / 1000.0
        xx, yy = np.meshgrid(np.linspace(-0.1, 0.1, 5), np.linspace(-0.1, 0.1, 5))
        zz = np.full(xx.shape, z_plane)
        self.ax_geo.plot_wireframe(xx, yy, zz, color='gray', linewidth=0.5, alpha=0.5)
        
        self.ax_geo.set_title("Diffraction Geometry")
        self.ax_geo.set_xlim(-1, 1)
        self.ax_geo.set_ylim(-1, 1)
        self.ax_geo.set_zlim(0, 2)
        self.ax_geo.view_init(elev=evelationviewangle, azim=azimuthalviewangle)
        # 
        
        self.canvas_geo.draw()
        
        # --- Crystal Structure View ---
        self.ax_crystal.clear()
        P = np.eye(3) * 3
        supercell = make_supercell(self.structure, P)
        positions = supercell.get_positions()
        rotpos = np.array([self.U @ pos for pos in positions])
        self.ax_crystal.scatter(rotpos[:,0], rotpos[:,1], rotpos[:,2],
                                 color='gray', s=20)
        scale = 10
        self.ax_crystal.quiver(0, 0, 0, *(scale * (self.U @ np.array([1, 0, 0]))),
                               color='lightskyblue', arrow_length_ratio=0.1, label='100')
        self.ax_crystal.quiver(0, 0, 0, *(scale * (self.U @ np.array([0, 1, 0]))),
                               color='mistyrose', arrow_length_ratio=0.1, label='010')
        self.ax_crystal.quiver(0, 0, 0, *(scale * (self.U @ np.array([0, 0, 1]))),
                               color='magenta', arrow_length_ratio=0.1, label='001')
        self.ax_crystal.quiver(0, 0, 0, *(scale * k_in_dir),
                               color='b', arrow_length_ratio=0.1, label='k_in')
        self.ax_crystal.quiver(0, 0, 0, *(scale * k_out_dir),
                               color='r', arrow_length_ratio=0.1, label='k_out')
        self.ax_crystal.quiver(0, 0, 0, *(q0_dir * 10),
                               color='g', arrow_length_ratio=0.1, label='q0')
        self.ax_crystal.set_title("Crystal Structure")
        self.ax_crystal.set_xlim(-20, 20)
        self.ax_crystal.set_ylim(-20, 20)
        self.ax_crystal.set_zlim(-20, 20)
        # self.ax_crystal.view_init(elev=evelationviewangle, azim=azimuthalviewangle)
        self.canvas_crystal.draw()

    def openStrainWindow(self):
        self.strain_window = StrainWindow(self.update_strain, self)
        self.strain_window.show()
    
    def update_strain(self, new_strain):
        self.current_strain = new_strain
        self.structure.cell = self.original_cell.dot(np.eye(3) + new_strain)
        self.update()
        
    def closeEvent(self, a0):
        QApplication.closeAllWindows()
        
# ----------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------
if __name__ == "__main__":
    import sys

    # Read crystal structure from CIF file
    structure = read('.\\ciffiles\\LiF_COD9008667.cif') 

    app = QApplication(sys.argv)
    setFusionpalette(app)
    main_window = CrystalReflectionVisualizer(
        structure, max_index=5, xtal_to_det_distance_meters=0.1,
        minimumEnergy=5000, maximumEnergy=70000, phi=0,
        max2th=50, theta=0, psi=0
    )
    main_window.show()
    sys.exit(app.exec_())
