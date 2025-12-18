import sys
import numpy as np
from ase.io import read
from ase.build import make_supercell
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget,
    QLineEdit, QLabel, QHBoxLayout, QDoubleSpinBox,
    QGridLayout
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPalette, QColor
from scipy.spatial.transform import Rotation as R
import spglib
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D  # needed for 3d plotting
import matplotlib.cm as cm
from scipy.interpolate import interp1d
 
from xoppylib.scattering_functions.xoppy_calc_f0 import xoppy_calc_f0 
import xraylib
from xoppylib.scattering_functions.f1f2_calc import f1f2_calc
import multiprocessing
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

def generate_q0(h, k, l, structure, wavelength, U):
    ''' returns q0 in inverse angstrom '''
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
    except Exception:
        return 0
    reciprocal_lattice = structure.cell.reciprocal()  # Reciprocal lattice vectors
    hkl_vector = np.dot([h, k, l], reciprocal_lattice)  # in crystal frame
    hkl_vector = np.dot(U, hkl_vector)  # rotated to lab frame
    q_magnitude = 2 * np.pi / d_spacing
    q0 = hkl_vector / np.linalg.norm(hkl_vector) * q_magnitude
    #
    return q0

def sinAngleBetweenTwoVectors(a, b):
    a = np.array(np.float64(a))
    b = np.array(np.float64(b))
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
    reciprocal_lattice = structure.cell.reciprocal()
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
        
        # Use the provided reciprocal lattice and U to get the unit q-vector
        qdir = generate_q_unit_vector(h, k, l, reciprocal_lattice, U)
        if d_spacing != 0 and qdir[2] != 0:
            qnorm = 2 * np.pi / d_spacing
            qvec = qnorm * qdir
            qx, qy, qz = qvec
            k0 = -(qnorm**2) / (2 * qz)
            kvec = k0 * np.array([0, 0, 1])
            ks = kvec + qvec
            
            if k0 > 0:
                # try:
                sinang = sinAngleBetweenTwoVectors(ks, kvec)
                if sinang >0:
                    E_bragg = hc / (2 * d_spacing * np.sin(np.arcsin(sinang)/2))
                    sintwoth = sinAngleBetweenTwoVectors(ks, kvec)
                    costwoth = cosAngleBetweenTwoVectors(ks, kvec)
                    twoth = np.arcsin(sintwoth)
                    if np.degrees(twoth) < max2th and costwoth > 0:
                        if minE < E_bragg < maxE:
                            F_hkl = calculate_structure_factor(structure, h, k, l, E_bragg, f0_values)
                            if F_hkl > 0.01:
                                k0mag = np.linalg.norm(k0)
                                ksmag = np.linalg.norm(ks)
                                
                                # print(hkl, f": (ks-k0)/k0 = {(ksmag-k0mag)/k0mag:.4f}")
                                # print('two-theta: ', twoth)
                                # print('d-spacing: ', d_spacing)
                                # print('q_z:', qz)
                                # print(E_bragg, hc/E_bragg)
                                # print(kvec)
                                # print(ks)
                                # print(qvec)
                                E_bragg_list.append(E_bragg)
                                F_hkl_list.append(F_hkl)
                                goodhkllist.append(hkl)
                                twothlist.append(twoth)
                                dspacinglist.append(d_spacing)
                                qveclist.append(qvec)
                                kslist.append(ks)
                                kveclist.append(kvec)
                            
                # except Exception as e:
                #     print(e, f': sintwoth = {sintwoth}')
                #     pass
    return dspacinglist, twothlist, E_bragg_list, kveclist, kslist, qveclist, goodhkllist, F_hkl_list

def VisualizeSingle(structure, hkl, U):
    a, b, c, alpha, beta, gamma = structure.cell.cellpar() 
    reciprocal_lattice = structure.cell.reciprocal()
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
    else:
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

def generate_q_unit_vector(h, k, l, reciprocal_lattice, U):
    hkl_vector = np.dot([h, k, l], reciprocal_lattice)
    hkl_vector = np.dot(U, hkl_vector)
    s_unitvec = hkl_vector / np.linalg.norm(hkl_vector)
    return s_unitvec

def genRotation(phi, theta, psi):
    rotation = R.from_euler('ZXZ', [phi, theta, psi], degrees=True)
    return rotation.as_matrix()

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

 
class CrystalReflectionVisualizer(QMainWindow):
    def __init__(self, structure, max_index=5, xtal_to_det_distance_meters=0.1,
                 minimumEnergy=18000, maximumEnergy=20000, phi=0,
                 max2th=50, theta=0, psi=0):
        super().__init__()
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.max_index = max_index
        self.structure = structure
        self.max2th = max2th
        self.minimumEnergy = minimumEnergy
        self.maximumEnergy = maximumEnergy
        self.base_hkl_list = [(h, k, l) for h in range(-self.max_index, self.max_index + 1)
                              for k in range(-self.max_index, self.max_index + 1)
                              for l in range(-self.max_index, self.max_index + 1)
                              if h != 0 or k != 0 or l != 0]
        
        self.axlim = 86
        self.detdistance_mm = xtal_to_det_distance_meters * 1000

        self.setWindowTitle("Crystal Reflection Visualizer")
        self.setGeometry(50, 50, 1600, 800)

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        
        layoutGrid = QGridLayout(self.central_widget)
        layoutDetView = QVBoxLayout()
        layoutGrid.addLayout(layoutDetView, 0, 0, 1, 2)
        
        # Create the Detector View (2D) using matplotlib
        self.detector_fig = Figure(figsize=(5, 4))
        self.detector_canvas = FigureCanvas(self.detector_fig)
        layoutDetView.addWidget(self.detector_canvas)
        # Connect the pick event for interactivity (for text labels)
        self.detector_canvas.mpl_connect('pick_event', self.on_text_pick)
        
        # Create the Diffraction Geometry view (3D) using matplotlib
        self.view1_fig = Figure(figsize=(5, 4))
        self.view1_ax = self.view1_fig.add_subplot(111, projection='3d')
        self.view1_canvas = FigureCanvas(self.view1_fig)
        layoutGrid.addWidget(self.view1_canvas, 0, 2, 1, 6)
        
        # Create the Crystal Structure view (3D) using matplotlib
        self.view2_fig = Figure(figsize=(5, 4))
        self.view2_ax = self.view2_fig.add_subplot(111, projection='3d')
        self.view2_canvas = FigureCanvas(self.view2_fig)
        layoutGrid.addWidget(self.view2_canvas, 0, 8, 1, 6)
        
        # Angle editors
        def makeangleedit(label,startval = 0.0):
            hbox = QHBoxLayout()
            layoutDetView.addLayout(hbox)
            hbox.addWidget(QLabel(label))
            spinner = QDoubleSpinBox()
            spinner.setMinimum(-180)
            spinner.setMaximum(180)
            spinner.setValue(startval)
            spinner.valueChanged.connect(self.update)
            hbox.addWidget(spinner)
            return spinner

        self.phi_spinner = makeangleedit('Phi',startval = phi)
        self.theta_spinner = makeangleedit('Theta',startval = theta)
        self.psi_spinner = makeangleedit('Psi',startval = psi)

        def addgeneraledit(label, startvalstr):
            hbox = QHBoxLayout()
            layoutDetView.addLayout(hbox)
            hbox.addWidget(QLabel(label))
            edit = QLineEdit(startvalstr)
            edit.editingFinished.connect(self.update)
            hbox.addWidget(edit)
            return edit

        self.detdistedit = addgeneraledit('Det Dist (mm)', '100')
        self.minimumEnergyedit = addgeneraledit('minE (eV)', f'{self.minimumEnergy}')
        self.maximumEnergyedit = addgeneraledit('maxE (eV)', f'{self.maximumEnergy}')
        self.axrangeedit = addgeneraledit('Axis Range (mm)', f'{self.axlim}')

        self.update()
        # Use the first good reflection for geometry visualization
        startindex = 0
        self.energyplotted = self.energies[startindex]
        self.wavelength = hc / self.energyplotted  * 1e-10  # ~1 Angstrom scale
        self.hkl = self.hkl_list[startindex]
        self.qhkl = self.q_hkl_list[startindex]
        self.visualize_geometry_and_crystal()
        self.visualize_crystal_structure()

    def update(self):
        try:
            self.detdistance_mm = float(self.detdistedit.text())
            self.minimumEnergy = float(self.minimumEnergyedit.text())
            self.maximumEnergy = float(self.maximumEnergyedit.text())
            self.axlim = float(self.axrangeedit.text())
        except:
            pass

        self.phi = self.phi_spinner.value()
        self.theta = self.theta_spinner.value()
        self.psi = self.psi_spinner.value()
        self.U = genRotation(self.phi, self.theta, self.psi)

        res = compute_listof_d_spacings_and_angles(self.structure, self.base_hkl_list, self.U, self.max2th,
                                                     minE=self.minimumEnergy, maxE=self.maximumEnergy)
        (dspacinglist, twothlist, E_bragg_list, k0, k_prime_list,
         q_hkl_list, goodhkllist, F_hkl_list) = res

        # Save computed quantities
        self.k0 = k0
        self.hkl_list = goodhkllist
        self.q_hkl_list = q_hkl_list
        self.k_prime_list = k_prime_list
        self.energies = E_bragg_list
        self.F_hkl_list = F_hkl_list

        if len(E_bragg_list) > 0:
            self.energy_colors = (np.array(E_bragg_list) - min(E_bragg_list)) / (max(E_bragg_list) - min(E_bragg_list) + 1e-9)
            
            msizeindicator = F_hkl_list#np.log10(F_hkl_list)
            
            self.markersizes =(np.array(msizeindicator) - min(msizeindicator)) / (max(msizeindicator) - min(msizeindicator) + 1e-9) 
            
            self.draw_det_reflections()
        else:
            print('No reflections')
            self.detector_fig.clf()
            self.detector_canvas.draw()

    def draw_det_reflections(self):
        # Clear the 2D detector axes and set labels/limits.
        self.detector_fig.clf()
        ax = self.detector_fig.add_subplot(111)
        ax.set_xlabel("Detector Position (mm)")
        ax.set_ylabel("Detector Position (mm)")
        ax.set_xlim(-self.axlim, self.axlim)
        ax.set_ylim(-self.axlim, self.axlim)
        ax.invert_yaxis()
        
        # Draw each reflection marker and add an interactive text label.
        for energy, hkl, q, k_prime, norm_color,markerscale,F in zip(
                self.energies, self.hkl_list, self.q_hkl_list, self.k_prime_list, self.energy_colors,self.markersizes,self.F_hkl_list):
            ks_x, ks_y, ks_z = k_prime
            det_x = ks_x * self.detdistance_mm / ks_z
            det_y = ks_y * self.detdistance_mm / ks_z
            color = cm.plasma(norm_color)
            # Draw the marker.
            ax.plot(det_x, det_y, 'o', color=color, markersize=4+1.5*markerscale)
            # Add the text label with picker enabled.
            text_obj = ax.text(det_x, det_y, str(hkl), color=color, picker=True, fontsize = 4+1.5*markerscale)
            # Store custom attributes on the text object.
            text_obj._hkl = hkl
            text_obj._energy = energy
            text_obj._q = q
            text_obj._Fhkl = F
        self.detector_canvas.draw()

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
        self.visualize_geometry_and_crystal()
        self.visualize_crystal_structure()

    def visualize_geometry_and_crystal(self):
        # print('wavelength', self.wavelength)
        # Compute q0 and beam vectors using the first reflection.
        # q0 = generate_q0(self.hkl[0], self.hkl[1], self.hkl[2], self.structure, self.wavelength*1e10, self.U)
        q = self.qhkl
        q_dir = q / np.linalg.norm(q)
        k_in_magnitude = 2 * np.pi / self.wavelength /1e10
        k_in = np.array([0, 0, k_in_magnitude])
        k_out = k_in + q
        k_in_dir = k_in / np.linalg.norm(k_in)
        k_out_mag = np.linalg.norm(k_out)
        k_out_dir = k_out / k_out_mag
        
        # print(k_in)
        # print(k_out)
        # print(q)
        
        # Save these directions for use in the crystal view.
        
        self.q0_dir = q_dir
        self.k_in_dir = k_in_dir
        self.k_out_dir = k_out_dir
        self.qdir = q_dir

        ax = self.view1_ax
        ax.cla()
        ax.set_title(f"Diffraction Geometry, hkl=[{self.hkl[0]},{self.hkl[1]},{self.hkl[2]}], {self.energyplotted:.0f} keV")
        origin = np.array([0, 0, 0])
        # Plot arrows using quiver.
        # print(f"deltak/k: {(k_out_mag-k_in_magnitude)/k_in_magnitude:.4f}")
        kscaling = 1/k_in_magnitude*self.detdistance_mm / 1000.0
        ax.quiver(*origin, *k_in*kscaling,   color='blue', label='k_in')
        ax.quiver(*origin, *k_out*kscaling,  color='red', label='k_out')
        ax.quiver(*k_in * kscaling, *q * kscaling,  color='green', label='q0')
        ax.set_box_aspect([0.2, 0.2, 0.2])
        
        # Draw a simple grid (detector plane) at z = detdistance (converted to meters).
        grid_z = self.detdistance_mm / 1000.0
        grid_size = 0.2
        grid_spacing = 0.01
        x_vals = np.arange(-grid_size, grid_size, grid_spacing)
        y_vals = np.arange(-grid_size, grid_size, grid_spacing)
        X, Y = np.meshgrid(x_vals, y_vals)
        Z = np.full_like(X, grid_z)
        ax.plot_wireframe(X, Y, Z, color='gray', linewidth=0.5)
        
        ax.set_xlim(-0.1, 0.1)
        ax.set_ylim(-0.1, 0.1)
        ax.set_zlim(-0.1, 0.1)
        
        ax.view_init(elev=-85., azim=50,roll=-140)
        ax.invert_yaxis()
        
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.legend()
        self.view1_canvas.draw()

    def visualize_crystal_structure(self):
        ax = self.view2_ax
        ax.cla()
        ax.set_title("Crystal Structure Orientation")
        P = np.eye(3) * 5
        supercell = make_supercell(self.structure, P)
        positions = supercell.get_positions()
        rotpos = np.array([self.U @ pos for pos in positions])
        centercell = np.mean(rotpos,axis=0)
        lengths = np.linalg.norm(centercell-rotpos,axis=1)
        centercell = rotpos[np.argmin(lengths)]
        ax.scatter(rotpos[:,0], rotpos[:,1], rotpos[:,2], s=5, color='gray')
        
        scale = 10
        dir_100 = scale * (self.U @ np.array([1, 0, 0]))
        dir_010 = scale * (self.U @ np.array([0, 1, 0]))
        dir_001 = scale * (self.U @ np.array([0, 0, 1]))
        origin = np.array([0, 0, 0])
        ax.quiver(*origin, *dir_100, length=1, color='lightskyblue', label='100')
        ax.quiver(*origin, *dir_010, length=1, color='orange', label='010')
        ax.quiver(*origin, *dir_001, length=1, color='magenta', label='001')
        ax.quiver(*origin, *(self.q0_dir * 10), length=1, color='green', label='q0')
        ax.quiver(*origin, *(self.k_in_dir * scale), length=1, color='blue', label='k_in')
        ax.quiver(*origin, *(self.k_out_dir * scale), length=1, color='red', label='k_out')
        
        plane_size =3
        u = np.linspace(-plane_size, plane_size, 10)
        v = np.linspace(-plane_size, plane_size, 10)
        U_plane, V_plane = np.meshgrid(u, v)
        W_plane = -(self.qdir[0] * U_plane + self.qdir[1] * V_plane) / self.qdir[2]
        
        ax.plot_surface(U_plane+centercell[0], V_plane+centercell[1], W_plane+centercell[2], alpha=0.5, color="orange", label="Plane Normal to q0")
        ax.view_init(elev=-85., azim=50,roll=-140)
        ax.invert_yaxis()
        # ax.set_box_aspect([1.0, 1.0, 1.0])
        
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
        ax.set_zlim(-20, 20)
        
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.legend()
        self.view2_canvas.draw()

    def closeEvent(self, a0):
        QApplication.closeAllWindows()

# ----------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------
if __name__ == "__main__":
    # Read crystal structure from CIF file (update the path as needed)
    ciffile = r"C:\Users\sinclair\Documents\GitHub\pyXRDImage\ciffiles\LiF_COD9008667.cif"
    structure = read(ciffile) 

    app = QApplication(sys.argv)
    setFusionpalette(app)
    main_window = CrystalReflectionVisualizer(
        structure, max_index=5, xtal_to_det_distance_meters=0.1,
        minimumEnergy=6000, maximumEnergy=37000, phi=-90,
        max2th=70, theta=142, psi=0
    )
    main_window.show()
    sys.exit(app.exec_())
