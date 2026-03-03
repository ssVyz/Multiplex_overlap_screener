from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
    QSpinBox, QCheckBox, QGroupBox, QDialogButtonBox, QLabel,
)


class SettingsDialog(QDialog):
    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumWidth(400)
        self._settings = dict(settings)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # --- Analysis Parameters ---
        analysis_group = QGroupBox("Analysis Parameters")
        form = QFormLayout()

        self.min_overlap_spin = QSpinBox()
        self.min_overlap_spin.setRange(2, 15)
        self.min_overlap_spin.setValue(self._settings.get("min_overlap", 3))
        form.addRow("Minimum overlap length:", self.min_overlap_spin)

        self.max_overlap_spin = QSpinBox()
        self.max_overlap_spin.setRange(3, 20)
        self.max_overlap_spin.setValue(self._settings.get("max_overlap", 10))
        form.addRow("Maximum overlap length:", self.max_overlap_spin)

        self.max_mm_spin = QSpinBox()
        self.max_mm_spin.setRange(0, 5)
        self.max_mm_spin.setValue(self._settings.get("max_mismatches", 1))
        form.addRow("Maximum mismatches:", self.max_mm_spin)

        self.ambiguity_check = QCheckBox("Treat ambiguous bases as matches if any variation could match")
        self.ambiguity_check.setChecked(self._settings.get("consider_ambiguity", False))
        form.addRow(self.ambiguity_check)

        analysis_group.setLayout(form)
        layout.addWidget(analysis_group)

        # --- Risk Level Thresholds ---
        risk_group = QGroupBox("Risk Level Thresholds")
        risk_form = QFormLayout()

        risk_form.addRow(QLabel("HIGH risk — overlap length >= and mismatches <="))

        self.high_min_ol_spin = QSpinBox()
        self.high_min_ol_spin.setRange(1, 20)
        self.high_min_ol_spin.setValue(self._settings.get("high_risk_min_overlap", 4))
        risk_form.addRow("  Min overlap length:", self.high_min_ol_spin)

        self.high_max_mm_spin = QSpinBox()
        self.high_max_mm_spin.setRange(0, 5)
        self.high_max_mm_spin.setValue(self._settings.get("high_risk_max_mismatches", 0))
        risk_form.addRow("  Max mismatches:", self.high_max_mm_spin)

        risk_form.addRow(QLabel("MEDIUM risk — overlap length >= and mismatches <="))

        self.med_min_ol_spin = QSpinBox()
        self.med_min_ol_spin.setRange(1, 20)
        self.med_min_ol_spin.setValue(self._settings.get("medium_risk_min_overlap", 2))
        risk_form.addRow("  Min overlap length:", self.med_min_ol_spin)

        self.med_max_mm_spin = QSpinBox()
        self.med_max_mm_spin.setRange(0, 5)
        self.med_max_mm_spin.setValue(self._settings.get("medium_risk_max_mismatches", 1))
        risk_form.addRow("  Max mismatches:", self.med_max_mm_spin)

        risk_group.setLayout(risk_form)
        layout.addWidget(risk_group)

        # --- Buttons ---
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def get_settings(self):
        """Return the edited settings dict."""
        return {
            "min_overlap": self.min_overlap_spin.value(),
            "max_overlap": self.max_overlap_spin.value(),
            "max_mismatches": self.max_mm_spin.value(),
            "consider_ambiguity": self.ambiguity_check.isChecked(),
            "high_risk_min_overlap": self.high_min_ol_spin.value(),
            "high_risk_max_mismatches": self.high_max_mm_spin.value(),
            "medium_risk_min_overlap": self.med_min_ol_spin.value(),
            "medium_risk_max_mismatches": self.med_max_mm_spin.value(),
        }
