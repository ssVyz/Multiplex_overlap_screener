from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem,
    QPushButton, QDialog, QPlainTextEdit, QLabel, QFileDialog, QMessageBox,
)
from PySide6.QtCore import Signal
from PySide6.QtGui import QColor

from Bio import SeqIO
from io import StringIO

from core.analysis import parse_fasta_text


class AddOligoDialog(QDialog):
    """Dialog for adding new oligos via FASTA text or file import."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Oligo/s")
        self.setMinimumSize(500, 400)
        self._build_ui()
        self.records = []

    def _build_ui(self):
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Paste sequences in FASTA format:"))

        self.text_edit = QPlainTextEdit()
        self.text_edit.setPlaceholderText(">Oligo1\nATCGATCGATCG\n>Oligo2\nGCTAGCTAGCTA")
        layout.addWidget(self.text_edit)

        btn_row = QHBoxLayout()

        import_btn = QPushButton("Import FASTA file...")
        import_btn.clicked.connect(self._import_fasta)
        btn_row.addWidget(import_btn)

        btn_row.addStretch()

        add_btn = QPushButton("Add")
        add_btn.setDefault(True)
        add_btn.clicked.connect(self._on_add)
        btn_row.addWidget(add_btn)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)

        layout.addLayout(btn_row)

    def _import_fasta(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA file", "",
            "FASTA files (*.fasta *.fa *.fas *.fna);;All files (*)"
        )
        if path:
            try:
                with open(path, "r", encoding="utf-8-sig") as f:
                    content = f.read()
                self.text_edit.setPlainText(content)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to read file:\n{e}")

    def _on_add(self):
        text = self.text_edit.toPlainText().strip()
        if not text:
            QMessageBox.warning(self, "Empty input", "Please enter or import FASTA sequences.")
            return

        records = parse_fasta_text(text)
        if not records:
            QMessageBox.warning(
                self, "Parse error",
                "No valid FASTA sequences found.\n\nExpected format:\n>Name\nATCGATCG..."
            )
            return

        self.records = records
        self.accept()


class OligoPanel(QWidget):
    """Left panel: oligo list with add/delete/import/export controls."""

    oligos_changed = Signal()       # emitted when the mix changes
    selection_changed = Signal(str)  # emitted with oligo id (or "" for none)

    # Risk-level background colours
    RISK_COLORS = {
        "HIGH": QColor(255, 120, 120, 90),
        "MEDIUM": QColor(255, 200, 80, 90),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        layout.addWidget(QLabel("Oligo Mix"))

        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QListWidget.SingleSelection)
        self.list_widget.currentItemChanged.connect(self._on_selection)
        layout.addWidget(self.list_widget)

        btn_row = QHBoxLayout()

        add_btn = QPushButton("New Oligo/s")
        add_btn.clicked.connect(self._on_add)
        btn_row.addWidget(add_btn)

        del_btn = QPushButton("Delete")
        del_btn.clicked.connect(self._on_delete)
        btn_row.addWidget(del_btn)

        export_btn = QPushButton("Export FASTA")
        export_btn.clicked.connect(self._on_export)
        btn_row.addWidget(export_btn)

        layout.addLayout(btn_row)

    # ---- public API ----

    def set_sequences(self, sequences):
        """Replace the full sequence list and refresh the list widget."""
        self._sequences = list(sequences)
        self._refresh_list()

    def get_sequences(self):
        return list(self._sequences)

    def update_risk_highlights(self, risk_map):
        """risk_map: dict mapping oligo_id -> highest risk level string."""
        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)
            oligo_id = item.data(1)  # stored in UserRole
            risk = risk_map.get(oligo_id, "NONE")
            color = self.RISK_COLORS.get(risk)
            if color:
                item.setBackground(color)
            else:
                item.setBackground(QColor(0, 0, 0, 0))

    # ---- private slots ----

    def _on_add(self):
        dlg = AddOligoDialog(self)
        if dlg.exec() == QDialog.Accepted and dlg.records:
            self._sequences.extend(dlg.records)
            self._refresh_list()
            self.oligos_changed.emit()

    def _on_delete(self):
        current = self.list_widget.currentRow()
        if current < 0:
            return
        self._sequences.pop(current)
        self._refresh_list()
        self.oligos_changed.emit()

    def _on_export(self):
        if not self._sequences:
            QMessageBox.information(self, "Nothing to export", "The oligo list is empty.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export FASTA", "oligos.fasta",
            "FASTA files (*.fasta);;All files (*)"
        )
        if path:
            try:
                with open(path, "w") as f:
                    SeqIO.write(self._sequences, f, "fasta")
                QMessageBox.information(self, "Exported", f"Saved {len(self._sequences)} sequences to:\n{path}")
            except Exception as e:
                QMessageBox.critical(self, "Export error", str(e))

    def _on_selection(self, current, _previous):
        if current is None:
            self.selection_changed.emit("")
        else:
            self.selection_changed.emit(current.data(1))

    # ---- helpers ----

    def _refresh_list(self):
        self.list_widget.clear()
        if not hasattr(self, "_sequences"):
            self._sequences = []
        for seq in self._sequences:
            item = QListWidgetItem(f"{seq.id}  ({len(seq.seq)} bp)")
            item.setData(1, seq.id)  # store id for lookup
            self.list_widget.addItem(item)
