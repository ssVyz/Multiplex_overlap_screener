from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QTextEdit, QHeaderView, QSplitter,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont

from core.analysis import MODE_OVERLAP, MODE_THERMO


class DetailPanel(QWidget):
    """Right panel: shows interactions (top) and overlap visualization (bottom)."""

    RISK_COLORS = {
        "HIGH": QColor(255, 120, 120, 90),
        "MEDIUM": QColor(255, 200, 80, 90),
    }

    # column headers per screening mode; the first and last are shared
    COLUMNS = {
        MODE_OVERLAP: ["Partner", "Overlap", "Mismatches", "Risk"],
        MODE_THERMO: ["Partner", "ΔG 3'", "ΔG worst", "ΔG total", "Risk"],
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._interactions = []
        self._current_oligo_id = None
        self._mode = MODE_OVERLAP
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Vertical)

        # Top: interaction table
        top_widget = QWidget()
        top_layout = QVBoxLayout(top_widget)
        top_layout.setContentsMargins(0, 0, 0, 0)

        self.title_label = QLabel("Select an oligo to view interactions")
        top_layout.addWidget(self.title_label)

        self.table = QTableWidget(0, 4)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        self.table.currentCellChanged.connect(self._on_row_selected)
        self._apply_columns(MODE_OVERLAP)
        top_layout.addWidget(self.table)
        splitter.addWidget(top_widget)

        # Bottom: visualization
        bottom_widget = QWidget()
        bottom_layout = QVBoxLayout(bottom_widget)
        bottom_layout.setContentsMargins(0, 0, 0, 0)

        self.viz_label = QLabel("Overlap Visualization")
        self.viz_label.setStyleSheet("font-weight: bold;")
        bottom_layout.addWidget(self.viz_label)

        self.viz_text = QTextEdit()
        self.viz_text.setReadOnly(True)
        self.viz_text.setFont(QFont("Courier New", 10))
        bottom_layout.addWidget(self.viz_text)
        splitter.addWidget(bottom_widget)

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter)

    def _apply_columns(self, mode):
        """Set the table's columns for the given screening mode."""
        headers = self.COLUMNS[mode]
        self.table.setColumnCount(len(headers))
        self.table.setHorizontalHeaderLabels(headers)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        for col in range(1, len(headers)):
            header.setSectionResizeMode(col, QHeaderView.ResizeToContents)

    def _partner_name(self, r, oligo_id):
        if r["primer1_id"] == oligo_id and r["primer2_id"] == oligo_id:
            return f"{r['primer1_name']} (self-dimer)"
        if r["primer1_id"] == oligo_id:
            return r["primer2_name"]
        return r["primer1_name"]

    def _row_values(self, r, oligo_id):
        """Cell text for one interaction, excluding the trailing risk column."""
        if self._mode == MODE_THERMO:
            return [
                self._partner_name(r, oligo_id),
                f"{r['dg_3prime']:.1f}",
                f"{r['dg_min']:.1f}",
                f"{r['dg_ens']:.1f}",
            ]
        return [
            self._partner_name(r, oligo_id),
            str(r["overlap_length"]),
            str(r["mismatches"]),
        ]

    def show_interactions(self, oligo_id, oligo_name, interactions, mode=MODE_OVERLAP):
        """Populate the table with interactions for the given oligo."""
        self._interactions = interactions
        self._current_oligo_id = oligo_id
        if mode != self._mode:
            self._mode = mode
            self._apply_columns(mode)
        self.title_label.setText(
            f"Interactions for: {oligo_name}" if oligo_name
            else "Select an oligo to view interactions"
        )
        self.viz_label.setText(
            "Duplex Visualization" if mode == MODE_THERMO else "Overlap Visualization"
        )
        self.viz_text.clear()

        self.table.setRowCount(len(interactions))
        for row, r in enumerate(interactions):
            values = self._row_values(r, oligo_id) + [r["risk_level"]]
            color = self.RISK_COLORS.get(r["risk_level"])
            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                if col:
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                if color:
                    item.setBackground(color)
                self.table.setItem(row, col, item)

        # Auto-select the first interaction to show its visualization
        if interactions:
            self.table.setCurrentCell(0, 0)
            self._on_row_selected(0, 0, -1, -1)

    def clear(self):
        self.table.setRowCount(0)
        self.viz_text.clear()
        self._interactions = []
        self._current_oligo_id = None
        self.title_label.setText("Select an oligo to view interactions")

    def _summary_line(self, r):
        if self._mode == MODE_THERMO:
            return (
                f"ΔG 3' anchored: {r['dg_3prime']:.2f} kcal/mol  "
                f"({r['primer1_name']} 3': {r['dg_3prime_1']:.2f}, "
                f"{r['primer2_name']} 3': {r['dg_3prime_2']:.2f})\n"
                f"ΔG strongest: {r['dg_min']:.2f} kcal/mol  |  "
                f"ΔG total (ensemble): {r['dg_ens']:.2f} kcal/mol  "
                f"over {r['n_structures']} structures  |  "
                f"Risk: {r['risk_level']}"
            )
        return (
            f"Overlap: {r['overlap_length']} bp  |  "
            f"Mismatches: {r['mismatches']}  |  "
            f"Risk: {r['risk_level']}"
        )

    def _on_row_selected(self, row, _col, _prev_row, _prev_col):
        if 0 <= row < len(self._interactions):
            r = self._interactions[row]
            self.viz_text.setPlainText(
                f"{self._summary_line(r)}\n\n{r['visualization']}"
            )
