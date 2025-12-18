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
import pyqtgraph as pg
import pyqtgraph.opengl as gl
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
    reciprocal_lattice = structure.cell.reciprocal()  # Reciprocal lattice vectors
    hkl_vector = np.dot([h, k, l], reciprocal_lattice)  # in crystal frame
    hkl_vector = np.dot(U, hkl_vector)  # rotated to lab frame
    q_magnitude = 2 * np.pi / d_spacing
    q0 = hkl_vector / np.linalg.norm(hkl_vector) * q_magnitude
    return q0

def sinAngleBetweenTwoVectors(a, b):
    a = np.array(a)
    b = np.array(b)
    val = np.linalg.norm(np.cross(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    if 1.01>val>1:
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
                    print(e,f': sintwoth = {sintwoth}')
                    pass
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
    # Use provided reciprocal_lattice rather than recomputing from structure.
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

# ----------------------------------------------------------------
# Pyqtgraph-based visualization classes
# ----------------------------------------------------------------
class CrystalReflectionVisualizer(QMainWindow):
    def __init__(self, structure, max_index=5, xtal_to_det_distance_meters=0.1,
                 minimumEnergy=5000, maximumEnergy=100000, phi=0,
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
        
        self.axlim = 80
        self.detdistance_mm = xtal_to_det_distance_meters * 1000

        self.setWindowTitle("3D Crystal Reflection Visualizer")
        self.setGeometry(50, 50, 1200, 800)

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        
        layoutGrid = QGridLayout(self.central_widget)
        layoutDetView = QVBoxLayout()
        layoutGrid.addLayout(layoutDetView,0, 0, 1,2)
        
        
        
        # self.graphics_layout = pg.GraphicsLayoutWidget()
        # layoutDetView.addWidget(self.graphics_layout)
        self.plotwidget_detView = pg.PlotWidget(title="Detector View (Peak Locations)")
        layoutDetView.addWidget(self.plotwidget_detView)
        
        self.view1 = gl.GLViewWidget()
        self.view2 = gl.GLViewWidget()
        # self.view1.setBackgroundColor('w')
        # self.view2.setBackgroundColor('w')
        layoutGrid.addWidget(self.view1,0, 2, 1,6)
        layoutGrid.addWidget(self.view2,0, 8, 1,6)
        
        # Angle editors
        def makeangleedit(label):
            hbox = QHBoxLayout()
            layoutDetView.addLayout(hbox)
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
        startindex = 0
        Estart = self.energies[startindex]
        self.wavelength = hc/Estart*1e-10  # 1 Angstrom
        self.hkl = self.hkl_list[0]
        self.visualize_geometry_and_crystal()
        # self.xtalviewwidget = XRDVisualizationApp(self.U, self.hkl_list[startindex], self.structure, startwavel, self.detdistance_mm/1000)
        # layout.addWidget(self.xtalviewwidget)

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

        self.k0 = k0
        self.hkl_list = goodhkllist
        self.q_hkl_list = q_hkl_list
        self.k_prime_list = k_prime_list
        self.energies = E_bragg_list

        if len(E_bragg_list) > 0:
            self.energy_colors = self.normalize_colors(E_bragg_list)
            self.draw_det_reflections()
        else:
            print('No reflections')
            self.plotwidget_detView.clear()

    def normalize_colors(self, energies):
        energy_min = min(energies)
        energy_max = max(energies)
        return (np.array(energies) - energy_min) / (energy_max - energy_min)

    def draw_det_reflections(self):
        self.plotwidget_detView.clear()
        self.plotwidget_detView.setLabel('left', 'Detector Position (mm)')
        self.plotwidget_detView.setLabel('bottom', 'Detector Position (mm)')
        self.plotwidget_detView.setXRange(-self.axlim, self.axlim)
        self.plotwidget_detView.setYRange(-self.axlim, self.axlim)
        plasma = pg.colormap.get('plasma')
        for k0, hkl, q, k_prime, norm_color in zip(
            self.k0, self.hkl_list, self.q_hkl_list, self.k_prime_list, self.energy_colors
        ):
            ks_x, ks_y, ks_z = k_prime
            det_x = ks_x * self.detdistance_mm / ks_z
            det_y = ks_y * self.detdistance_mm / ks_z
            qcolor = plasma.map(norm_color, mode='qcolor')
            # Plot a marker at the detector hit position
            self.plotwidget_detView.plot([det_x], [det_y], pen=None, symbol='o', symbolBrush=qcolor, symbolSize=10)
            # Add a text label for the hkl
            text = pg.TextItem(text=str(hkl), color=qcolor.getRgb()[0:3], anchor=(0.5, -0.5))
            text.setPos(det_x, det_y)
            self.plotwidget_detView.addItem(text)

    def visualize_geometry_and_crystal(self):
        # Compute q0 and beam vectors
        q0 = generate_q0(self.hkl[0], self.hkl[1], self.hkl[2], self.structure, self.wavelength, self.U)
        q0_dir = q0 / np.linalg.norm(q0)
        k_in_magnitude = 2 * np.pi / self.wavelength
        k_in = np.array([0, 0, k_in_magnitude])
        k_out = k_in + q0
        k_in_dir = k_in / np.linalg.norm(k_in)
        k_out_dir = k_out / np.linalg.norm(k_out)

        # Left view: Diffraction Geometry
        axes = gl.GLAxisItem()
        axes.setSize(1, 1, 1)
        self.view1.addItem(axes)
        self.addArrow(self.view1, np.array([0, 0, 0]), k_in_dir, color=(0, 0, 255, 255), width=2)
        self.addArrow(self.view1, np.array([0, 0, 0]), k_out_dir, color=(255, 0, 0, 255), width=2)
        self.addArrow(self.view1, np.array([0, 0, 0]), q0_dir, color=(0, 255, 0, 255), width=2)
        # Add a simple grid as detector plane at z = xtal_to_det_distance
        grid = gl.GLGridItem()
        grid.setSize(0.2, 0.2)
        grid.setSpacing(0.01, 0.01)
        grid.rotate(90, 1, 0, 0)
        grid.translate(0, 0, self.detdistance_mm/1000)
        self.view1.addItem(grid)
        self.view1.setCameraPosition(distance=2)

        # Right view: Crystal Structure (supercell)
        P = np.eye(3) * 5
        supercell = make_supercell(self.structure, P)
        positions = supercell.get_positions()
        rotpos = np.array([self.U @ pos for pos in positions])
        # rotpos = np.array([1,1,0])
        sp = gl.GLScatterPlotItem(pos=rotpos, size=5, color=(0.5, 0.5, 0.5, 1), pxMode=True)
        self.view2.addItem(sp)
        scale = 10
        dir_100 = scale * (self.U @ np.array([1, 0, 0]))
        dir_010 = scale * (self.U @ np.array([0, 1, 0]))
        dir_001 = scale * (self.U @ np.array([0, 0, 1]))
        origin = np.array([0, 0, 0])
        
        self.addArrow(self.view2, origin, dir_100, color=QColor('lightskyblue'), width=2)
        self.addArrow(self.view2, origin, dir_010, color=QColor('mistyrose'), width=2)
        self.addArrow(self.view2, origin, dir_001, color=QColor('magenta'), width=4)
        self.addArrow(self.view2, origin, q0_dir * 10, color=(0, 255, 0, 255), width=2)
        self.addArrow(self.view2, origin, k_in_dir * scale, color=(0, 0, 255, 255), width=2)
        self.addArrow(self.view2, origin, k_out_dir * scale, color=(255, 0, 0, 255), width=2)

        # Create a plane normal to q0 using a mesh
        drawplane = False
        if drawplane: 
            plane_size = 1
            u = np.linspace(-plane_size, plane_size, 10)
            v = np.linspace(-plane_size, plane_size, 10)
            U_plane, V_plane = np.meshgrid(u, v)
            W_plane = -(q0_dir[0] * U_plane + q0_dir[1] * V_plane) / q0_dir[2]
            vertices = np.c_[U_plane.ravel(), V_plane.ravel(), W_plane.ravel()]
            faces = []
            nrows, ncols = U_plane.shape
            
            for i in range(nrows - 1):
                for j in range(ncols - 1):
                    idx = i * ncols + j
                    faces.append([idx, idx + 1, idx + ncols])
                    faces.append([idx + 1, idx + ncols + 1, idx + ncols])
            faces = np.array(faces)
            mesh = gl.GLMeshItem(vertexes=vertices, faces=faces, faceColors=(1, 0.5, 0, 0.5),
                                 smooth=False, drawEdges=False,shader=None)
            self.view2.addItem(mesh)
        self.view2.setCameraPosition(distance=100)

    def addArrow(self, view, start, direction, color=(255, 255, 255, 255), width=2):
        # Draw an arrow as a simple line from start to start+direction.
        end = start + direction
        pts = np.array([start, end])
        arrow = gl.GLLinePlotItem(pos=pts, color=color, width=width, antialias=True)
        view.addItem(arrow)
    def closeEvent(self, a0):
        app.closeAllWindows()
        
# class XRDVisualizationApp(QWidget):
#     def __init__(self, U, hkl, structure, wavelength, xtal_to_det_distance):
#         super().__init__()
        
#         self.U = U
#         self.hkl = hkl
#         self.structure = structure
#         self.wavelength = wavelength
#         self.detdistance_mm/1000 = xtal_to_det_distance
        
#         layout = QHBoxLayout()
#         self.setLayout(layout)

#         # Create two GLViewWidgets for side-by-side 3D views
        

#         self.visualize_geometry_and_crystal()

    

# ----------------------------------------------------------------
# Main execution
# ----------------------------------------------------------------
if __name__ == "__main__":
    import sys

    # Read crystal structure from CIF file
    structure = read('.\\ciffiles\\Ge.cif') 

  

    app = QApplication(sys.argv)
    setFusionpalette(app)
    # Launch the 2D reflection visualizer
    main_window = CrystalReflectionVisualizer(
        structure, max_index=5, xtal_to_det_distance_meters=0.1,
        minimumEnergy=5000, maximumEnergy=100000, phi=0,
        max2th=50, theta=0, psi=0
    )
    main_window.show()

    # Optionally, you can also launch the 3D visualization.
    # For example, using the first good reflection (if available):
    # (dspacinglist, twothlist, E_bragg_list, k0, k_prime_list, q_hkl_list, goodhkllist, F_hkl_list) = 
    #     compute_listof_d_spacings_and_angles(structure, some_hkl_list, genRotation(phi, theta, psi), max2th,
    #                                          minE=minimumEnergy, maxE=maximumEnergy)
    # if goodhkllist:
    #     vis3d = XRDVisualizationApp(genRotation(phi, theta, psi), goodhkllist[0], structure, wavelength, xtal_to_det_distance)
    #     vis3d.show()

    sys.exit(app.exec_())
