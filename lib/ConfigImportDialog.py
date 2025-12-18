from PyQt5.QtWidgets import QDialog, QVBoxLayout, QLabel, QCheckBox, QDialogButtonBox

class ConfigImportDialog(QDialog):
    def __init__(self, available_components, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import Configuration")
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(QLabel("Found the following configuration items.\nSelect components to import:"))
        
        self.checkboxes = {}
        for comp in available_components:
            cb = QCheckBox(comp)
            cb.setChecked(True)
            self.layout.addWidget(cb)
            self.checkboxes[comp] = cb
            
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        self.layout.addWidget(self.buttons)
        
    def get_selected(self):
        return [comp for comp, cb in self.checkboxes.items() if cb.isChecked()]
