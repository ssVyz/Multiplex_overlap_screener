from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem,
    QTextEdit, QHeaderView, QSplitter,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont


class DetailPanel(QWidget):
    """Right panel: shows interactions for the selected oligo."""

    RISK_COLORS = {
        "HIGH": QColor(255, 120, 120, 90),
        "MEDIUM": QColor(255, 200, 80, 90),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._interactions = []
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.title_label = QLabel("Select an oligo to view interactions")
        layout.addWidget(self.title_label)

        splitter = QSplitter(Qt.Vertical)

        # Interaction table
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Partner", "Overlap", "Mismatches", "Risk"])
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setSelectionMode(QTableWidget.SingleSelection)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.table.currentCellChanged.connect(self._on_row_selected)
        splitter.addWidget(self.table)

        # Visualization area
        self.viz_text = QTextEdit()
        self.viz_text.setReadOnly(True)
        self.viz_text.setFont(QFont("Courier New", 10))
        splitter.addWidget(self.viz_text)

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter)

    def show_interactions(self, oligo_id, interactions):
        """Populate the table with interactions for the given oligo."""
        self._interactions = interactions
        self.title_label.setText(f"Interactions for: {oligo_id}" if oligo_id else "Select an oligo to view interactions")
        self.viz_text.clear()

        self.table.setRowCount(len(interactions))
        for row, r in enumerate(interactions):
            partner = r["primer2_id"] if r["primer1_id"] == oligo_id else r["primer1_id"]
            if r["primer1_id"] == oligo_id and r["primer2_id"] == oligo_id:
                partner = f"{oligo_id} (self-dimer)"

            items = [
                QTableWidgetItem(partner),
                QTableWidgetItem(str(r["overlap_length"])),
                QTableWidgetItem(str(r["mismatches"])),
                QTableWidgetItem(r["risk_level"]),
            ]

            color = self.RISK_COLORS.get(r["risk_level"])
            for col, item in enumerate(items):
                if color:
                    item.setBackground(color)
                self.table.setItem(row, col, item)

    def clear(self):
        self.table.setRowCount(0)
        self.viz_text.clear()
        self._interactions = []
        self.title_label.setText("Select an oligo to view interactions")

    def _on_row_selected(self, row, _col, _prev_row, _prev_col):
        if 0 <= row < len(self._interactions):
            r = self._interactions[row]
            text = (
                f"Overlap: {r['overlap_length']} bp  |  "
                f"Mismatches: {r['mismatches']}  |  "
                f"Risk: {r['risk_level']}\n\n"
                f"{r['visualization']}"
            )
            self.viz_text.setPlainText(text)
