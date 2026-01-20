import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel, QMessageBox, QMenuBar, QAction, QMenu
)
from PyQt5.QtCore import QTimer, Qt
from xoppylib.xoppy_xraylib_util import descriptor_kind_index, nist_compound_list
import numpy as np

class filterRow(QWidget):
    def __init__(self, compound="", thickness=0.0, density=0.0,parent=None):
        super(filterRow, self).__init__(parent)
        
        # Create layout for the row
        self.layout = QHBoxLayout()

        # Checkbox for enabling/disabling
        self.checkbox = QCheckBox()
        self.layout.addWidget(self.checkbox)

        # Input fields for Compound, Thickness, and Density
        if compound:
            self.compound_edit = QLineEdit(compound)
        else:
            self.compound_edit = QLineEdit()
        self.compound_edit.setPlaceholderText("Compound")
        self.compound_edit.setContextMenuPolicy(Qt.CustomContextMenu)
        self.compound_edit.customContextMenuRequested.connect(self.show_context_menu)
        self.layout.addWidget(self.compound_edit)
        
        
        if self.isvalidvalue(thickness):
            self.thickness_edit = QLineEdit(str(thickness))
        else:
            self.thickness_edit = QLineEdit()
        self.thickness_edit.setPlaceholderText("Thickness (mm)")
        self.layout.addWidget(self.thickness_edit)

        if self.isvalidvalue(density):
            self.density_edit = QLineEdit(str(density))
        else:
            self.density_edit = QLineEdit()
        self.density_edit.setPlaceholderText("Density (g/cc)")
        self.layout.addWidget(self.density_edit)

        # Delete button
        self.delete_button = QPushButton("Delete")
        self.delete_button.clicked.connect(self.delete_row)
        self.layout.addWidget(self.delete_button)
        
        
        

        # Connect signals to update the lists in the parent
        self.compound_edit.editingFinished.connect(lambda: parent.update_lists())
        self.thickness_edit.editingFinished.connect(lambda: parent.update_lists())
        self.density_edit.editingFinished.connect(lambda: parent.update_lists())
        self.checkbox.toggled.connect(lambda: parent.update_lists())
        
        
        self.setLayout(self.layout)

    def show_context_menu(self, position):
        menu = QMenu()
        nist_menu = menu.addMenu("NIST compounds")
        
        try:
            compounds = nist_compound_list()
            for compound in compounds:
                action = QAction(compound, self)
                action.triggered.connect(lambda checked, text=compound: self.set_compound(text))
                nist_menu.addAction(action)
        except Exception as e:
            print(f"Error loading NIST compounds: {e}")

        menu.exec_(self.compound_edit.mapToGlobal(position))

    def set_compound(self, text):
        self.compound_edit.setText(text)
        self.density_edit.setText('?')
        self.parent().update_lists()
        
    def isvalidvalue(self,value):
        if value=='?':
            return 1
        try:
            if float(value)>0:
                return 1
            else:
                return 0
        except:
            return 0
                
        
    def get_values(self):
        if self.checkbox.isChecked():
            return (
                self.compound_edit.text(),
                self.thickness_edit.text(),
                self.density_edit.text()
            )
        return None
    def delete_row(self):
        self.parent().remove_row(self)
        
    def blink_red(self,widgettoblink):
        # Set the edit box to red
        widgettoblink.setStyleSheet("QLineEdit { background-color: red; }")
    
        # Create a timer to revert the color back
        QTimer.singleShot(500, lambda: self.revert_color(widgettoblink))  # Change color back after 500 ms

    def revert_color(self,widgettoblink):
        # Revert to the original style
        widgettoblink.setStyleSheet("")

class FilterManager(QWidget):
    def __init__(self,parent_obj):
        super().__init__()
        self.parent_obj = parent_obj
        self.setWindowTitle("Filter Manager")
        self.layout = QVBoxLayout()
        
        # Menu Bar
        menu_bar = QMenuBar(self)
        # For QWidget, add the menu bar to the layout directly
        self.layout.setMenuBar(menu_bar) # This line is problematic for QWidget, it should be self.layout.addWidget(menu_bar)
        # Corrected: Add menu_bar to the layout
        self.layout.addWidget(menu_bar)

        file_menu = menu_bar.addMenu('File')
        
        load_std_action = QAction('Load Standard Filters', self)
        load_std_action.triggered.connect(self.load_standard_filters)
        file_menu.addAction(load_std_action)
        
        

        button_row = QHBoxLayout()
        # Button to add new filters
        self.add_button = QPushButton("Add filter")
        self.add_button.clicked.connect(self.add_filter)
        button_row.addWidget(self.add_button)
        self.plot_button = QPushButton("Plot attenuation")
        self.plot_button.clicked.connect(self.plot_attenuation)
        button_row.addWidget(self.plot_button)
        self.layout.addLayout(button_row)

        
        # Container for rows
        self.filters_container = QVBoxLayout()
        
        # # Create layout for the row
        # self.labellayout = QHBoxLayout()        
        # checkbox = QCheckBox()
        # self.labellayout.addWidget(checkbox)
        # # checkbox.setVisible(1)        
        # self.labellayout.addWidget(QLabel('Compound'))
        
        # self.filters_container.addLayout(self.labellayout)
        self.layout.addLayout(self.filters_container)

        self.setLayout(self.layout)
        self.activeCompounds = []
        self.activeThicknesses = []
        self.activeDensities = []
        

    def load_standard_filters(self):
        # Based on user's standard set
        self.add_filter("C", 0.2, 3.53)
        self.add_filter("Be", 1.27, '?')
        self.add_filter("Kapton Polyimide Film", 0.25, '?')
        self.add_filter("He", 300.0, '?')
        self.add_filter("C0.000124N0.755268O0.231781Ar0.012827", 100.0, 0.00120479)
        
    def add_filter(self, compound="", thickness=0.0, density=0.0,enable=False):
        new_row = filterRow(compound, thickness, density, self)
        self.filters_container.addWidget(new_row)
        if enable:
            new_row.checkbox.setCheckState(2)
            self.update_lists()
            

        
    def isvalidcompound(self,compound):
        
        try:
            kind = descriptor_kind_index(compound)
        except:
            return 0         
        if not (kind==-1):
            return 1
            

    def remove_row(self, row):
        row.deleteLater()  # Safely remove the row
        self.filters_container.removeWidget(row)  # Remove from layout
        row.setParent(None)  # Clear parent
        
    def update_lists(self):
        # Clear the lists
        self.activeCompounds.clear()
        self.activeThicknesses.clear()
        self.activeDensities.clear()

        # Iterate through the rows and update the lists based on checked boxes
        for i in range(self.filters_container.count()):
            row = self.filters_container.itemAt(i).widget()
            values = row.get_values()
            badval = 0
            if values:
                compound, thickness, density = values
                if compound and thickness and density:
                    if self.isvalidvalue(thickness) and self.isvalidvalue(density) and self.isvalidcompound(compound):
                        
                        
                        self.activeCompounds.append(compound)
                        self.activeThicknesses.append(float(thickness))
                        if density=='?':
                            self.activeDensities.append(density)
                        else:
                            self.activeDensities.append(float(density))
                    else:
                        badval=1
                else:
                    badval=1
                if badval:
                    row.checkbox.setCheckState(0)
                    if not self.isvalidvalue(thickness):
                        row.blink_red(row.thickness_edit)
                    if not self.isvalidvalue(density):
                        row.blink_red(row.density_edit)  
                    if not self.isvalidcompound(compound):
                        row.blink_red(row.compound_edit)
            else:
                badval=1
            
                          
            

        # Optionally print the lists to see the updates
        print("Active Compounds:", self.activeCompounds)
        print("Active Thicknesses:", self.activeThicknesses)
        print("Active Densities:", self.activeDensities)        
        
    def isvalidvalue(self,value):
        if value=='?':
            return 1
        try:
            if float(value)>0:
                return 1
            else:
                return 0
        except:
            return 0
        

    def plot_attenuation(self):
        """Plot attenuation using the parent DCSSpectrum window."""
        required_attrs = ("XrayFilter", "plot_filter_attenuation", "subwindows")
        if any(not hasattr(self.parent_obj, attr) for attr in required_attrs):
            QMessageBox.warning(
                self,
                "Unavailable",
                "This plot is only available when the Filter Manager is opened from DCSSpectrum.",
            )
            return

        self.update_lists()
        if not self.activeCompounds:
            QMessageBox.information(self, "No Active Filters", "Enable at least one valid filter to plot attenuation.")
            return

        try:
            emin = float(self.parent_obj.Emin_edit.text())
            emax = float(self.parent_obj.Emax_edit.text())
            npts = int(float(self.parent_obj.NPts_edit.text()))
        except (AttributeError, ValueError):
            QMessageBox.warning(self, "Invalid Energy Range", "Enter valid energy range values in the main GUI.")
            return

        if npts < 2 or emax <= emin:
            QMessageBox.warning(self, "Invalid Energy Range", "Ensure Emin < Emax and # Points is at least 2.")
            return

        energy = np.linspace(emin, emax, npts)
        flat_input = np.ones_like(energy)

        try:
            filtered_energy, filtered_power = self.parent_obj.XrayFilter(energy, flat_input)
        except Exception as exc:
            QMessageBox.critical(self, "Filter Error", f"Unable to evaluate filter attenuation:\n{exc}")
            return

        self.parent_obj.plot_filter_attenuation(filtered_energy, filtered_power)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FilterManager(app)
    window.resize(600, 200)
    
    # Example: Adding some initial rows with values
    window.add_filter("Air", 1, 0.02)

    window.show()
    sys.exit(app.exec_())
