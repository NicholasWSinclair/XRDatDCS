import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QHBoxLayout,
    QCheckBox, QLineEdit, QLabel, QGroupBox, QStyle
)
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QTimer
from xoppylib.xoppy_xraylib_util import descriptor_kind_index

class PostSampleFilterRow(QWidget):
    def __init__(self, compound="", thickness=0.0, density=0.0, angle=0.0, parent=None):
        super(PostSampleFilterRow, self).__init__(parent)
        
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
        self.compound_edit.setToolTip("Compound")
        self.layout.addWidget(self.compound_edit)
        
        
        if self.isvalidvalue(thickness):
            self.thickness_edit = QLineEdit(str(thickness))
        else:
            self.thickness_edit = QLineEdit()
        self.thickness_edit.setPlaceholderText("Thickness (mm)")
        self.thickness_edit.setToolTip("Thickness (mm)")
        self.layout.addWidget(self.thickness_edit)

        if self.isvalidvalue(density):
            self.density_edit = QLineEdit(str(density))
        else:
            self.density_edit = QLineEdit()
        self.density_edit.setPlaceholderText("Density (g/cc)")
        self.density_edit.setToolTip("Density (g/cc)")
        self.layout.addWidget(self.density_edit)
        
        if self.isvalidangle(angle):
            self.angle_edit = QLineEdit(str(angle))
        else:
            self.angle_edit = QLineEdit()
        self.angle_edit.setPlaceholderText("Angle Norm. to Beam (deg)")
        self.angle_edit.setToolTip("Angle Norm. to Beam (deg)")
        self.layout.addWidget(self.angle_edit)

        # Delete button
        Trashicon = self.style().standardIcon(QStyle.SP_DialogDiscardButton)
        # Trashicon = self.style().standardIcon(QStyle.SP_TrashIcon)
        self.delete_button = QPushButton(icon=Trashicon)
        self.delete_button.clicked.connect(self.delete_row)
        self.layout.addWidget(self.delete_button)
        
        
        

        # Connect signals to update the lists in the parent
        self.compound_edit.editingFinished.connect(lambda: parent.update_lists())
        self.thickness_edit.editingFinished.connect(lambda: parent.update_lists())
        self.density_edit.editingFinished.connect(lambda: parent.update_lists())
        self.angle_edit.editingFinished.connect(lambda: parent.update_lists())
        self.checkbox.toggled.connect(lambda: parent.update_lists())
        
        self.setLayout(self.layout)
        
    def isvalidvalue(self,value):
        if value=='?':
            return 1
        try:
            if float(value)>=0:
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
                self.density_edit.text(),
                self.angle_edit.text()
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
    def isvalidangle(self,value):        
        try:
            if float(value)>=0:
                return 1
                
            else:
                return 0
        except:
            return 0               

class PostSampleFilterManager(QGroupBox):
    def __init__(self,parent_obj):
        super().__init__()
        self.setTitle('Post-Sample Filters')
        self.parent_obj = parent_obj
        
        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)
        

        # Button to add new filters
        self.add_button = QPushButton("Add filter")
        self.add_button.clicked.connect(self.add_filter)
        self.layout.addWidget(self.add_button)


        
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
        self.activeAngles = []
        

    def add_filter(self, compound="", thickness=0.0, density=0.0, angle=0.0):
        new_row = PostSampleFilterRow(compound, thickness, density, angle, self)
        self.filters_container.addWidget(new_row)
        
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
        self.activeAngles.clear()

        # Iterate through the rows and update the lists based on checked boxes
        for i in range(self.filters_container.count()):
            row = self.filters_container.itemAt(i).widget()
            values = row.get_values()
            badval = 0
            if values:
                compound, thickness, density,angle = values

                if self.isvalidvalue(thickness) and self.isvalidvalue(density) and self.isvalidangle(angle) and self.isvalidcompound(compound):
                    
                    
                    self.activeCompounds.append(compound)
                    self.activeThicknesses.append(float(thickness))
                    self.activeAngles.append(float(angle))
                    if density=='?':
                        self.activeDensities.append(density)
                    else:
                        self.activeDensities.append(float(density))
                else:
                    badval=1
                if badval:
                    row.checkbox.setCheckState(0)
                    if not self.isvalidvalue(thickness):
                        row.blink_red(row.thickness_edit)
                    if not self.isvalidvalue(density):
                        row.blink_red(row.density_edit)  
                    if not self.isvalidangle(angle):
                        row.blink_red(row.angle_edit)  
                    if not self.isvalidcompound(compound):
                        row.blink_red(row.compound_edit)
            else:
                badval=1
            
                          
            

        # Optionally print the lists to see the updates
        for C,T,D,A in zip(self.activeCompounds, self.activeThicknesses, self.activeDensities, self.activeAngles):
            if D=='?':
                print(f'{C}: Thickness: {T:.3f}, Rho: {D}, Angle: {A:.3f}')
            else:
                print(f'{C}: Thickness: {T:.3f}, Rho: {D:.3f}, Angle: {A:.3f}')
        
    def isvalidvalue(self,value):
        try:
            if float(value)>0:
                return 1
            else:
                return 0
        except:
            return 0
    def isvalidangle(self,value):        
        try:
            float(value)
            return 1
        except:
            return 0                  

 


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PostSampleFilterManager(app)
    window.resize(600, 200)
    
    # Example: Adding some initial rows with values
    window.add_filter("Air", 1, 0.02, 0.)

    window.show()
    sys.exit(app.exec_())
