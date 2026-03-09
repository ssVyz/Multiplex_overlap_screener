from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QTextEdit, QFormLayout
from PySide6.QtGui import QFont


class OligoPreview(QWidget):
    """Bottom-left panel: shows details for the currently selected oligo."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        self.title_label = QLabel("Oligo Preview")
        self.title_label.setStyleSheet("font-weight: bold;")
        layout.addWidget(self.title_label)

        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)

        self.name_label = QLabel("—")
        form.addRow("Name:", self.name_label)

        self.length_label = QLabel("—")
        form.addRow("Length:", self.length_label)

        self.gc_label = QLabel("—")
        form.addRow("GC%:", self.gc_label)

        self.tm_label = QLabel("—")
        form.addRow("Tm:", self.tm_label)

        self.status_label = QLabel("—")
        form.addRow("Status:", self.status_label)

        self.mix_label = QLabel("—")
        form.addRow("Mix:", self.mix_label)

        layout.addLayout(form)

        self.seq_display = QTextEdit()
        self.seq_display.setReadOnly(True)
        self.seq_display.setFont(QFont("Courier New", 10))
        self.seq_display.setMaximumHeight(120)
        layout.addWidget(self.seq_display)

        layout.addStretch()

    def set_settings(self, settings: dict):
        """Update settings used for Tm calculation."""
        self._na_concentration = settings.get("na_concentration", 50.0)

    def show_oligo(self, oligo, mix_name=None):
        """Display details for the given Oligo object."""
        self.title_label.setText(f"Oligo Preview: {oligo.name}")
        self.name_label.setText(oligo.name)
        self.length_label.setText(f"{len(oligo.sequence)} bp")
        self.gc_label.setText(f"{oligo.gc_content():.1f}%")
        na_mM = getattr(self, "_na_concentration", 50.0)
        tm = oligo.calc_tm(na_mM)
        self.tm_label.setText(f"{tm:.1f} \u00b0C")
        self.status_label.setText("Active" if oligo.active else "Inactive")
        self.mix_label.setText(mix_name if mix_name else "Unallocated")
        self.seq_display.setPlainText(oligo.sequence)

    def clear(self):
        self.title_label.setText("Oligo Preview")
        self.name_label.setText("—")
        self.length_label.setText("—")
        self.gc_label.setText("—")
        self.tm_label.setText("—")
        self.status_label.setText("—")
        self.mix_label.setText("—")
        self.seq_display.clear()
