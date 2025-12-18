import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas,NavigationToolbar2QT
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class GeometryVisualizer(QWidget):
    def __init__(self, sample_thickness=.1, sample_angle=70, filterlist = None, conelist = None, sample_to_detector_distance=100, pixel_size=.043, num_pixels=2048,detectorxycenter=(0,0)):
        super().__init__()
        self.sample_thickness = sample_thickness
        self.sample_angle = sample_angle
        # self.filter_thickness = filter_thickness
        # self.filter_angle = filter_angle
        self.sample_to_detector_distance = sample_to_detector_distance
        self.pixel_size = pixel_size
        self.num_pixels = num_pixels
        self.filterlist = filterlist
        self.setWindowTitle('Experiment Geometry')
        self.conelist = conelist
        self.detectorxycenter = detectorxycenter
        
        self.angleforlateralextent = 70
        
        self.filterzoffset = 0
        # self.filter_lateralextent = 100
        self.sample_lateralextent = 5
        
        self.init_ui()
        self.update_plot()
    
    def init_ui(self):
        layout = QVBoxLayout()
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas)
        self.figure.tight_layout()
        self.setLayout(layout)
    
    def draw_box(self, ax, center, dimensions, angle=0,facecolor='b',alpha = 1,label=None):
        """Helper function to draw a rotated 3D box."""
        

        x, y, z = center
        dx, dy, dz = dimensions
        # Define the 8 corners of the box
        corners = np.array([
            [-dx/2, -dy/2, -dz/2],
            [ dx/2, -dy/2, -dz/2],
            [ dx/2,  dy/2, -dz/2],
            [-dx/2,  dy/2, -dz/2],
            [-dx/2, -dy/2,  dz/2],
            [ dx/2, -dy/2,  dz/2],
            [ dx/2,  dy/2,  dz/2],
            [-dx/2,  dy/2,  dz/2],
        ])
        
        # Apply rotation about the Y axis
        rotation_matrix = np.array([
            [ np.cos(np.radians(angle)), 0, np.sin(np.radians(angle))],
            [ 0, 1, 0],
            [-np.sin(np.radians(angle)), 0, np.cos(np.radians(angle))],
        ])
        rotated_corners = np.dot(corners, rotation_matrix.T) + [x, y, z]
        
        # Define the 6 faces of the box using the rotated corners
        faces = [
            [rotated_corners[j] for j in [0, 1, 5, 4]],  # Bottom face
            [rotated_corners[j] for j in [2, 3, 7, 6]],  # Top face
            [rotated_corners[j] for j in [0, 3, 7, 4]],  # Left face
            [rotated_corners[j] for j in [1, 2, 6, 5]],  # Right face
            [rotated_corners[j] for j in [0, 1, 2, 3]],  # Front face
            [rotated_corners[j] for j in [4, 5, 6, 7]],  # Back face
        ]
        if label is None:
            box = ax.add_collection3d(Poly3DCollection(faces, alpha=alpha, linewidths=1, edgecolors='k'))
        else:
            box = ax.add_collection3d(Poly3DCollection(faces, alpha=alpha, linewidths=1, edgecolors='k', label=label))
        box.set_facecolor(facecolor)

        
    # def draw_cone(self, ax, start, theta):
    #    """Draw a cone and its intersection with the detector plane."""
    #    x0, y0, z0 = start
    #    tan_theta = np.tan(np.radians(theta))
    #    r = tan_theta * (self.sample_to_detector_distance - z0)
       
    #    # Draw the cone surface
    #    phi = np.linspace(0, 2 * np.pi, 50)
    #    x_circle = r * np.cos(phi)
    #    y_circle = r * np.sin(phi)
    #    z_circle = np.full_like(x_circle, self.sample_to_detector_distance)
       
    #    # Cone edges
    #    for phi_edge in np.linspace(0, 2 * np.pi, 30):
    #        x_edge = [x0, x0 + r * np.cos(phi_edge)]
    #        y_edge = [y0, y0 + r * np.sin(phi_edge)]
    #        z_edge = [z0, self.sample_to_detector_distance]
    #        ax.plot(x_edge, y_edge, z_edge, color='purple')
       
    #    # Intersection circle on the detector plane
    #    ax.plot(x_circle, y_circle, z_circle, color='blue', label=f'2-theta = {theta}{chr(176)}')    
    
    def draw_cone(self, ax, start, theta,alpha):
        """Draw a cone and its intersection with the detector plane using plot_surface."""
        x0, y0, z0 = start
        tan_theta = np.tan(np.radians(theta))
        h = self.sample_to_detector_distance - z0  # Height of the cone
        r = tan_theta * h  # Radius of the cone at the detector plane
    
        # Create the cone surface
        phi = np.linspace(0, 2 * np.pi, 20)  # Angle around the cone
        z = np.linspace(z0, self.sample_to_detector_distance, 3)  # Height of the cone
        Z, PHI = np.meshgrid(z, phi)  # Create a grid for the surface
        R = tan_theta * (Z - z0)  # Radius at each height
        
        if np.tan(np.pi/2-np.radians(self.sample_angle))<0:
            X = np.min([x0 + R * np.cos(PHI),-1*Z * np.tan(np.pi/2-np.radians(self.sample_angle))],0)  # X-coordinates
        else:
            X = np.max([x0 + R * np.cos(PHI),-1*Z * np.tan(np.pi/2-np.radians(self.sample_angle))],0)  # X-coordinates
        Y = y0 + R * np.sin(PHI)  # Y-coordinates
        
        
    
        ax.plot_surface(X, Y, Z, color='purple', alpha=alpha)
    
        # Draw the intersection circle on the detector plane
        circle_phi = np.linspace(0, 2 * np.pi, 30)
        # x_circle = x0 + r * np.cos(circle_phi)
        min_x_at_detdistance = x0 - self.sample_to_detector_distance * np.tan(np.pi/2-np.radians(self.sample_angle))
        minxarray = np.ones_like(circle_phi)*min_x_at_detdistance
        
        if np.tan(np.pi/2-np.radians(self.sample_angle))<0:
            x_circle =np.min([x0 + r * np.cos(circle_phi),minxarray],0)  # X-coordinates
        else:
            x_circle =np.max([x0 + r * np.cos(circle_phi),minxarray],0)  # X-coordinates
            

        y_circle = y0 + r * np.sin(circle_phi)
        z_circle = np.full_like(y_circle, self.sample_to_detector_distance)
        
        
        ax.plot(x_circle, y_circle, z_circle, color='blue', label=f'2{chr(952)} = {theta:.1f}{chr(176)}')

    
    def update_plot(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111, projection='3d')
        
        # Plot the beam (Z-axis line)
        ax.plot([0, 0], [0, 0], [-50, self.sample_to_detector_distance], color='red', label="Beam")
        
        # Plot the sample
        self.draw_box(ax, center=(0, 0, 0), dimensions=(self.sample_lateralextent, self.sample_lateralextent, self.sample_thickness), angle=self.sample_angle, facecolor = 'b',alpha =1)
        
        # Plot the post-sample filters
        if self.filterlist is not None:
            # z0 = self.filterzoffset 
            # + self.sample_thickness/np.cos(np.radians(self.sample_angle))
            # + self.sample_lateralextent/np.sin(np.radians(self.sample_angle))
            for filterdict in self.filterlist:
                
                filter_thickness=filterdict['thickness']
                filter_angle=filterdict['angle']
                zpos = filterdict['zpos']
                filter_center = (filter_thickness/2*np.sin(np.radians(filter_angle)), 0, zpos+filter_thickness/2*np.cos(np.radians(filter_angle)))
                # z0 = z0 + filter_thickness + 10
                filter_lateralextent= 2*(zpos)*np.tan(np.radians(self.angleforlateralextent)) + 2*filter_thickness*np.tan(np.radians(filter_angle))
                self.draw_box(ax, center=filter_center, dimensions=(filter_lateralextent, filter_lateralextent, filter_thickness), angle=filter_angle, facecolor = 'r',alpha=0.2,label=filterdict['mat'])
                
        # Plot the XRD Cones
        
        if self.conelist is not None:
            n = len(self.conelist)
            alpha = 0.2/n
            for coneangle in self.conelist:
                xyzstart = (0,0,0)
                self.draw_cone(ax, start=xyzstart, theta=coneangle,alpha=alpha)                
        
        # Plot the detector
        x0,y0 = self.detectorxycenter
        
        detector_width = self.pixel_size * self.num_pixels
        detector_corners = np.array([
            [x0-detector_width/2, y0-detector_width/2, self.sample_to_detector_distance],
            [ x0+detector_width/2, y0-detector_width/2, self.sample_to_detector_distance],
            [ x0+detector_width/2, y0+ detector_width/2, self.sample_to_detector_distance],
            [x0+-detector_width/2, y0+ detector_width/2, self.sample_to_detector_distance],
        ])
        ax.add_collection3d(Poly3DCollection([detector_corners], alpha=0.5, linewidths=1, edgecolors='k', facecolor='m'))
        
        # Set viewing angle
        ax.view_init(elev=0, azim=90)
        
        # Set axes limits and labels
        ax.set_xlim(-100, 100)
        ax.set_ylim(-100, 100)
        ax.set_zlim(-50, 1.05*self.sample_to_detector_distance)
        # ax.set_xlabel('X (mm)')
        # ax.set_ylabel('Y (mm)')
        # ax.set_zlabel('Z (mm)')
        ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        ax.set_aspect('equal')
        self.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    filterlist = [{'mat':'LiF','thickness':0.1,'angle':28,'zpos':0.0},{'mat':'Window','thickness':1,'angle':28+90,'zpos':10.0}]
    conelist = [10,20,30]
    window = GeometryVisualizer(filterlist = filterlist, conelist = conelist)
    window.show()
    sys.exit(app.exec_())
