from PySide6.QtWidgets import (
    QDialog, QVBoxLayout, QFormLayout, QComboBox,
    QSpinBox, QDoubleSpinBox, QCheckBox, QGroupBox, QDialogButtonBox, QLabel,
)

from core.analysis import MODE_OVERLAP, MODE_THERMO


class SettingsDialog(QDialog):

    MODES = [
        ("3' end overlap", MODE_OVERLAP),
        ("Thermodynamic (ΔG)", MODE_THERMO),
    ]

    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Settings")
        self.setMinimumWidth(400)
        self._settings = dict(settings)
        self._build_ui()
        self._on_mode_changed()

    def _dg_spin(self, key, default):
        """Spin box for a kcal/mol threshold (negative = more stable)."""
        spin = QDoubleSpinBox()
        spin.setRange(-40.0, 0.0)
        spin.setDecimals(1)
        spin.setSingleStep(0.5)
        spin.setSuffix(" kcal/mol")
        spin.setValue(self._settings.get(key, default))
        return spin

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # --- Screening mode ---
        mode_group = QGroupBox("Screening Mode")
        mode_form = QFormLayout()

        self.mode_combo = QComboBox()
        for label, value in self.MODES:
            self.mode_combo.addItem(label, value)
        current = self._settings.get("screen_mode", MODE_OVERLAP)
        index = self.mode_combo.findData(current)
        self.mode_combo.setCurrentIndex(index if index >= 0 else 0)
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        mode_form.addRow("Mode:", self.mode_combo)

        mode_group.setLayout(mode_form)
        layout.addWidget(mode_group)

        # --- Shared parameters ---
        shared_group = QGroupBox("Sequence Parameters")
        shared_form = QFormLayout()

        self.ambiguity_check = QCheckBox("Treat ambiguous bases as matches if any variation could match")
        self.ambiguity_check.setChecked(self._settings.get("consider_ambiguity", False))
        shared_form.addRow(self.ambiguity_check)

        self.na_conc_spin = QDoubleSpinBox()
        self.na_conc_spin.setRange(1.0, 1000.0)
        self.na_conc_spin.setDecimals(1)
        self.na_conc_spin.setSuffix(" mM")
        self.na_conc_spin.setValue(self._settings.get("na_concentration", 50.0))
        shared_form.addRow("Na⁺ concentration:", self.na_conc_spin)

        shared_group.setLayout(shared_form)
        layout.addWidget(shared_group)

        # --- Overlap mode ---
        self.overlap_group = QGroupBox("Overlap Parameters")
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

        form.addRow(QLabel("HIGH risk — overlap length >= and mismatches <="))

        self.high_min_ol_spin = QSpinBox()
        self.high_min_ol_spin.setRange(1, 20)
        self.high_min_ol_spin.setValue(self._settings.get("high_risk_min_overlap", 4))
        form.addRow("  Min overlap length:", self.high_min_ol_spin)

        self.high_max_mm_spin = QSpinBox()
        self.high_max_mm_spin.setRange(0, 5)
        self.high_max_mm_spin.setValue(self._settings.get("high_risk_max_mismatches", 0))
        form.addRow("  Max mismatches:", self.high_max_mm_spin)

        form.addRow(QLabel("MEDIUM risk — overlap length >= and mismatches <="))

        self.med_min_ol_spin = QSpinBox()
        self.med_min_ol_spin.setRange(1, 20)
        self.med_min_ol_spin.setValue(self._settings.get("medium_risk_min_overlap", 2))
        form.addRow("  Min overlap length:", self.med_min_ol_spin)

        self.med_max_mm_spin = QSpinBox()
        self.med_max_mm_spin.setRange(0, 5)
        self.med_max_mm_spin.setValue(self._settings.get("medium_risk_max_mismatches", 1))
        form.addRow("  Max mismatches:", self.med_max_mm_spin)

        self.overlap_group.setLayout(form)
        layout.addWidget(self.overlap_group)

        # --- Thermodynamic mode ---
        self.thermo_group = QGroupBox("ΔG Parameters")
        dg_form = QFormLayout()

        self.dg_temp_spin = QDoubleSpinBox()
        self.dg_temp_spin.setRange(0.0, 100.0)
        self.dg_temp_spin.setDecimals(1)
        self.dg_temp_spin.setSuffix(" °C")
        self.dg_temp_spin.setValue(self._settings.get("dg_temperature", 37.0))
        dg_form.addRow("Temperature:", self.dg_temp_spin)

        dg_form.addRow(QLabel("Report a pair if either ΔG is at or below"))

        self.report_dg_min_spin = self._dg_spin("report_dg_min", -2.0)
        dg_form.addRow("  Strongest interaction:", self.report_dg_min_spin)

        self.report_dg_ens_spin = self._dg_spin("report_dg_ens", -3.0)
        dg_form.addRow("  Total (ensemble):", self.report_dg_ens_spin)

        dg_form.addRow(QLabel("HIGH risk — ΔG at or below either of"))

        self.high_dg_3p_spin = self._dg_spin("high_risk_dg_3prime", -6.0)
        dg_form.addRow("  3' anchored:", self.high_dg_3p_spin)

        self.high_dg_any_spin = self._dg_spin("high_risk_dg_any", -10.0)
        dg_form.addRow("  Strongest interaction:", self.high_dg_any_spin)

        dg_form.addRow(QLabel("MEDIUM risk — ΔG at or below either of"))

        self.med_dg_3p_spin = self._dg_spin("medium_risk_dg_3prime", -4.0)
        dg_form.addRow("  3' anchored:", self.med_dg_3p_spin)

        self.med_dg_any_spin = self._dg_spin("medium_risk_dg_any", -7.0)
        dg_form.addRow("  Strongest interaction:", self.med_dg_any_spin)

        self.thermo_group.setLayout(dg_form)
        layout.addWidget(self.thermo_group)

        # --- Buttons ---
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _on_mode_changed(self):
        """Show only the parameter group belonging to the selected mode."""
        thermo = self.mode_combo.currentData() == MODE_THERMO
        self.overlap_group.setVisible(not thermo)
        self.thermo_group.setVisible(thermo)
        self.adjustSize()

    def get_settings(self):
        """Return the edited settings dict."""
        return {
            "screen_mode": self.mode_combo.currentData(),
            "min_overlap": self.min_overlap_spin.value(),
            "max_overlap": self.max_overlap_spin.value(),
            "max_mismatches": self.max_mm_spin.value(),
            "consider_ambiguity": self.ambiguity_check.isChecked(),
            "high_risk_min_overlap": self.high_min_ol_spin.value(),
            "high_risk_max_mismatches": self.high_max_mm_spin.value(),
            "medium_risk_min_overlap": self.med_min_ol_spin.value(),
            "medium_risk_max_mismatches": self.med_max_mm_spin.value(),
            "na_concentration": self.na_conc_spin.value(),
            "dg_temperature": self.dg_temp_spin.value(),
            "report_dg_min": self.report_dg_min_spin.value(),
            "report_dg_ens": self.report_dg_ens_spin.value(),
            "high_risk_dg_3prime": self.high_dg_3p_spin.value(),
            "high_risk_dg_any": self.high_dg_any_spin.value(),
            "medium_risk_dg_3prime": self.med_dg_3p_spin.value(),
            "medium_risk_dg_any": self.med_dg_any_spin.value(),
        }
