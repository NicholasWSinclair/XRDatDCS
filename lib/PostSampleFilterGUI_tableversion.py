import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel, QGroupBox, QStyle, QTableWidget, 
    QTableWidgetItem, QHeaderView, QSizePolicy
)
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QTimer
from xoppylib.xoppy_xraylib_util import descriptor_kind_index
from xoppylib.xoppy_xraylib_util import density as getdensity


class PostSampleFilterRow(QWidget):
    def __init__(self, compound="", thickness=0.0, density=0.0, angle=0.0, zpos=0.0, parent=None):
        super(PostSampleFilterRow, self).__init__(parent)
        
        # Initialize widgets
        self.checkbox = QCheckBox()
        self.checkbox.row_widget=self
        
        self.compound_edit = QLineEdit(compound)
        self.compound_edit.setPlaceholderText("Compound")
        
        
        self.thickness_edit = QLineEdit(str(thickness) if self.isvalidvalue_greaterthan0(thickness) else "")
        self.thickness_edit.setPlaceholderText("Thickness (mm)")
        
        self.density_edit = QLineEdit(str(density) if self.isvalidvalue_greaterthan0(density) else "")
        self.density_edit.setPlaceholderText("Density (g/cc)")
        
        self.angle_edit = QLineEdit(str(angle) if self.isvalidangle(angle) else "")
        self.angle_edit.setPlaceholderText("Angle Norm. to Beam (deg)")
        
        self.zpos_edit = QLineEdit(str(zpos) if self.isvalidvalue_greaterorequal0(zpos) else "")
        self.zpos_edit.setPlaceholderText("Angle Norm. to Beam (deg)")

        # Delete button (trash icon)
        Trashicon = self.style().standardIcon(QStyle.SP_DialogDiscardButton)
        self.delete_button = QPushButton(icon=Trashicon)
        self.delete_button.clicked.connect(self.delete_row)
        
        # Connect signals to update the lists in the parent
        self.compound_edit.editingFinished.connect(lambda: self.updatecompound(parent))
        self.thickness_edit.editingFinished.connect(lambda: parent.update_lists())
        self.density_edit.editingFinished.connect(lambda: parent.update_lists())
        self.angle_edit.editingFinished.connect(lambda: parent.update_lists())
        self.zpos_edit.editingFinished.connect(lambda: parent.update_lists())
        self.checkbox.toggled.connect(lambda: parent.update_lists())

        for edit in [self.compound_edit, self.thickness_edit, self.density_edit, self.angle_edit, self.zpos_edit]:
            edit.setAlignment(Qt.AlignCenter)
        # if parent is not None:
        #     parent.table.resizeColumnsToContents()
        
            
    def isvalidvalue_greaterthan0(self, value):
        try:
            return float(value) > 0
        except:
            return False

    def isvalidvalue_greaterorequal0(self, value):
        try:
            return float(value) >= 0
        except:
            return False    
    def updatecompound(self,parent):
        try:
            dens = getdensity(self.compound_edit.text())
            self.density_edit.setText(f'{dens:.3f}')
        except:
            pass
        parent.update_lists()
        

    def isvalidangle(self,value):    
        try:
            value = float(value)
            return True if float(value)>=-180 and float(value)<=180 else False
        except:
            return False 

    def get_values(self):
        if self.checkbox.isChecked():
            return (
                self.compound_edit.text(),
                self.thickness_edit.text(),
                self.density_edit.text(),
                self.angle_edit.text(),
                self.zpos_edit.text()
            )
        return None

    def delete_row(self):
        self.parent().remove_row(self)

    def blink_red(self, widgettoblink):
        widgettoblink.setStyleSheet("QLineEdit { background-color: red; }")
        QTimer.singleShot(500, lambda: self.revert_color(widgettoblink))  # Change color back after 500 ms

    def revert_color(self, widgettoblink):
        widgettoblink.setStyleSheet("")


class PostSampleFilterManager(QGroupBox):
    def __init__(self, parent_obj):
        super().__init__()
        self.setTitle('Post-Sample Filters')
        self.parent_obj = parent_obj
        
        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)

        # Button to add new filters
        self.add_button = QPushButton("Add filter")
        self.add_button.clicked.connect(lambda: self.add_filter(compound="", thickness=0.1, density=1.0, angle=0.0,zpos=0.0))
        self.layout.addWidget(self.add_button)

        # Table for managing filter rows
        self.table = QTableWidget(0, 7)  # 6 columns for checkbox, compound, thickness, density, angle, actions
        self.table.setHorizontalHeaderLabels(['', 'Compound', 'Thickness', 'Density', 'Angle', 'Z-pos','Del'])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        
        for col in [0,6]:
            self.table.setColumnWidth(col, 10)  # Set the width 
            self.table.horizontalHeader().setSectionResizeMode(col, QHeaderView.Fixed)
            
        for col in [4,5]:
            self.table.setColumnWidth(col, 75)  # Set the width 
            self.table.horizontalHeader().setSectionResizeMode(col, QHeaderView.Fixed)            
            
        for col in [2,3]:
            self.table.setColumnWidth(col, 75)  # Set the width 
            self.table.horizontalHeader().setSectionResizeMode(col, QHeaderView.Fixed)            
        
        
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        # self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.resizeColumnsToContents()
        self.layout.addWidget(self.table)

        self.setLayout(self.layout)

        # Initialize lists to store filter data
        self.activeCompounds = []
        self.activeThicknesses = []
        self.activeDensities = []
        self.activeAngles = []
        self.activezpos = []

    def add_filter(self, compound="", thickness=0.1, density=1.0, angle=0.0, zpos=0.0,checked=0):
        row_position = self.table.rowCount()
        self.table.insertRow(row_position)

        # Create a new row widget
        row = PostSampleFilterRow(compound=compound, thickness=thickness, density=density, angle=angle, zpos=zpos,parent=self)
        row.row_position=row_position

        # Add widgets to respective columns
        self.table.setCellWidget(row_position, 0, row.checkbox)
        if checked==1:
            row.checkbox.setCheckState(2)
        self.table.setCellWidget(row_position, 1, row.compound_edit)
        self.table.setCellWidget(row_position, 2, row.thickness_edit)
        self.table.setCellWidget(row_position, 3, row.density_edit)
        self.table.setCellWidget(row_position, 4, row.angle_edit)
        self.table.setCellWidget(row_position, 5, row.zpos_edit)
        self.table.setCellWidget(row_position, 6, row.delete_button)

    def remove_row(self, row):
        # row_position = self.table.indexFromItem(row.checkbox).row()
        self.table.removeRow(row.row_position)

    def update_lists(self):
        self.activeCompounds.clear()
        self.activeThicknesses.clear()
        self.activeDensities.clear()
        self.activeAngles.clear()
        self.activezpos.clear()

        # Update lists based on the rows in the table
        for row_index in range(self.table.rowCount()):
            
            checkboxwidget = self.table.cellWidget(row_index, 0)  # Get the row widget (checkbox is in column 0)
            row = checkboxwidget.row_widget
            values = row.get_values()  # Correctly call get_values() on the row widget
            if values:
                compound, thickness, density, angle,zpos = values

                if self.isvalidvalue_greaterthan0(thickness) and self.isvalidvalue_greaterthan0(density) and self.isvalidangle(angle) and self.isvalidcompound(compound) and self.isvalidvalue_greaterorequal0(zpos):
                    self.activeCompounds.append(compound)
                    self.activeThicknesses.append(float(thickness))
                    self.activeAngles.append(float(angle))
                    self.activezpos.append(float(zpos))
                    if density == '?':
                        self.activeDensities.append(density)
                    else:
                        self.activeDensities.append(float(density))
                else:
                    row.checkbox.setCheckState(0)  # Uncheck if any value is invalid
                    if not self.isvalidvalue_greaterthan0(thickness):
                        row.blink_red(row.thickness_edit)
                    if not self.isvalidvalue_greaterthan0(density):
                        row.blink_red(row.density_edit)
                    if not self.isvalidvalue_greaterorequal0(zpos):
                        row.blink_red(row.zpos_edit)                        
                    if not self.isvalidangle(angle):
                        row.blink_red(row.angle_edit)
                    if not self.isvalidcompound(compound):
                        row.blink_red(row.compound_edit)

        # Print the updated lists for debugging
        for C, T, D, A,Z in zip(self.activeCompounds, self.activeThicknesses, self.activeDensities, self.activeAngles,self.activezpos):
            if D == '?':
                print(f'{C}: Thickness: {T:.3f}, Rho: {D}, Angle: {A:.3f}, Z-Position: {Z:.3f}')
            else:
                print(f'{C}: Thickness: {T:.3f}, Rho: {D:.3f}, Angle: {A:.3f}, Z-Position: {Z:.3f}')

    def isvalidcompound(self, compound):
        try:
            kind = descriptor_kind_index(compound)
            return kind != -1
        except:
            return False
        
    def isvalidvalue_greaterthan0(self, value):
        try:
            return float(value) > 0
        except:
            return False

    def isvalidvalue_greaterorequal0(self, value):
        try:
            return float(value) >= 0
        except:
            return False

    def isvalidangle(self,value):    
        try:
            value = float(value)
            return True if float(value)>=-180 and float(value)<=180 else False
        except:
            return False 
                
            


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PostSampleFilterManager(app)
    window.resize(600, 200)

    # Example: Adding some initial rows with values
    window.add_filter("Air", 1, 0.02, 0.0)

    window.show()
    sys.exit(app.exec_())
