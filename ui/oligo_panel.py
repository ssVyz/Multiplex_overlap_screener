from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QListWidget, QListWidgetItem,
    QPushButton, QDialog, QPlainTextEdit, QLabel, QFileDialog, QMessageBox,
    QComboBox,
)
from PySide6.QtCore import Signal, Qt
from PySide6.QtGui import QColor, QFont

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

    RISK_ORDER = {"HIGH": 3, "MEDIUM": 2, "LOW": 1, "NONE": 0}

    def __init__(self, parent=None):
        super().__init__(parent)
        self._sequences = []
        self._deactivated = set()       # oligo IDs that are deactivated
        self._risk_map = {}             # oligo_id -> risk level string
        self._display_indices = []      # maps list-widget row -> _sequences index
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        layout.addWidget(QLabel("Oligo Mix"))

        # Sort / filter controls
        controls_row = QHBoxLayout()

        controls_row.addWidget(QLabel("Sort:"))
        self.sort_combo = QComboBox()
        self.sort_combo.addItems(["Original Order", "Name (A→Z)", "Name (Z→A)", "Risk (High first)"])
        self.sort_combo.currentIndexChanged.connect(self._on_sort_filter_changed)
        controls_row.addWidget(self.sort_combo)

        controls_row.addWidget(QLabel("Filter:"))
        self.filter_combo = QComboBox()
        self.filter_combo.addItems(["All", "Active Only", "High Risk", "Medium+ Risk"])
        self.filter_combo.currentIndexChanged.connect(self._on_sort_filter_changed)
        controls_row.addWidget(self.filter_combo)

        layout.addLayout(controls_row)

        # Oligo list
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QListWidget.ExtendedSelection)
        self.list_widget.currentItemChanged.connect(self._on_selection)
        self.list_widget.itemSelectionChanged.connect(self._on_selection_changed)
        layout.addWidget(self.list_widget)

        # Action buttons — row 1
        btn_row1 = QHBoxLayout()

        add_btn = QPushButton("New Oligo/s")
        add_btn.clicked.connect(self._on_add)
        btn_row1.addWidget(add_btn)

        del_btn = QPushButton("Delete Selected")
        del_btn.clicked.connect(self._on_delete)
        btn_row1.addWidget(del_btn)

        self.deactivate_btn = QPushButton("Deactivate")
        self.deactivate_btn.clicked.connect(self._on_toggle_active)
        btn_row1.addWidget(self.deactivate_btn)

        layout.addLayout(btn_row1)

        # Action buttons — row 2
        btn_row2 = QHBoxLayout()

        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(self._on_select_all)
        btn_row2.addWidget(select_all_btn)

        export_btn = QPushButton("Export FASTA")
        export_btn.clicked.connect(self._on_export)
        btn_row2.addWidget(export_btn)

        layout.addLayout(btn_row2)

    # ---- public API ----

    def set_sequences(self, sequences):
        """Replace the full sequence list and refresh the list widget."""
        self._sequences = list(sequences)
        self._deactivated.clear()
        self._risk_map.clear()
        self._refresh_list()

    def get_sequences(self):
        """Return all sequences (active + deactivated). Used for FASTA export."""
        return list(self._sequences)

    def get_active_sequences(self):
        """Return only active (non-deactivated) sequences. Used for analysis."""
        return [s for s in self._sequences if s.id not in self._deactivated]

    def update_risk_highlights(self, risk_map):
        """risk_map: dict mapping oligo_id -> highest risk level string."""
        self._risk_map = dict(risk_map)
        self._refresh_list()

    # ---- private slots ----

    def _on_add(self):
        dlg = AddOligoDialog(self)
        if dlg.exec() == QDialog.Accepted and dlg.records:
            self._sequences.extend(dlg.records)
            self._refresh_list()
            self.oligos_changed.emit()

    def _on_delete(self):
        selected_items = self.list_widget.selectedItems()
        if not selected_items:
            return

        count = len(selected_items)
        if count > 1:
            reply = QMessageBox.question(
                self, "Confirm Delete",
                f"Delete {count} selected oligo(s)?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return

        # Map selected list-widget rows back to _sequences indices
        rows = sorted(
            [self.list_widget.row(item) for item in selected_items],
            reverse=True,
        )
        for row in rows:
            if 0 <= row < len(self._display_indices):
                seq_idx = self._display_indices[row]
                oligo_id = self._sequences[seq_idx].id
                self._deactivated.discard(oligo_id)

        # Remove from _sequences in reverse order of original index to avoid shifting
        orig_indices = sorted(
            [self._display_indices[r] for r in rows if 0 <= r < len(self._display_indices)],
            reverse=True,
        )
        for idx in orig_indices:
            self._sequences.pop(idx)

        self._refresh_list()
        self.oligos_changed.emit()

    def _on_toggle_active(self):
        selected_items = self.list_widget.selectedItems()
        if not selected_items:
            return

        # Collect oligo IDs for the selected items
        selected_ids = []
        for item in selected_items:
            row = self.list_widget.row(item)
            if 0 <= row < len(self._display_indices):
                seq_idx = self._display_indices[row]
                selected_ids.append(self._sequences[seq_idx].id)

        # If all selected are deactivated, activate them; otherwise deactivate all
        all_deactivated = all(oid in self._deactivated for oid in selected_ids)
        if all_deactivated:
            for oid in selected_ids:
                self._deactivated.discard(oid)
        else:
            for oid in selected_ids:
                self._deactivated.add(oid)

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
            self.selection_changed.emit(current.data(Qt.UserRole))

    def _on_selection_changed(self):
        """Update the deactivate/activate button text based on current selection."""
        selected_items = self.list_widget.selectedItems()
        if not selected_items:
            self.deactivate_btn.setText("Deactivate")
            return

        selected_ids = []
        for item in selected_items:
            row = self.list_widget.row(item)
            if 0 <= row < len(self._display_indices):
                seq_idx = self._display_indices[row]
                selected_ids.append(self._sequences[seq_idx].id)

        if selected_ids and all(oid in self._deactivated for oid in selected_ids):
            self.deactivate_btn.setText("Activate")
        else:
            self.deactivate_btn.setText("Deactivate")

    def _on_select_all(self):
        self.list_widget.selectAll()

    def _on_sort_filter_changed(self):
        self._refresh_list()

    # ---- helpers ----

    def _refresh_list(self):
        self.list_widget.clear()
        self._display_indices = []

        if not self._sequences:
            return

        # Build indexed list: (original_index, sequence)
        indexed = list(enumerate(self._sequences))

        # --- Apply filter ---
        filter_mode = self.filter_combo.currentText()
        if filter_mode == "Active Only":
            indexed = [(i, s) for i, s in indexed if s.id not in self._deactivated]
        elif filter_mode == "High Risk":
            indexed = [(i, s) for i, s in indexed
                       if self._risk_map.get(s.id) == "HIGH"]
        elif filter_mode == "Medium+ Risk":
            indexed = [(i, s) for i, s in indexed
                       if self._risk_map.get(s.id) in ("HIGH", "MEDIUM")]

        # --- Apply sort ---
        sort_mode = self.sort_combo.currentText()
        if sort_mode == "Name (A→Z)":
            indexed.sort(key=lambda x: x[1].id.lower())
        elif sort_mode == "Name (Z→A)":
            indexed.sort(key=lambda x: x[1].id.lower(), reverse=True)
        elif sort_mode == "Risk (High first)":
            indexed.sort(
                key=lambda x: self.RISK_ORDER.get(self._risk_map.get(x[1].id, "NONE"), 0),
                reverse=True,
            )

        # --- Populate list widget ---
        for orig_idx, seq in indexed:
            is_deactivated = seq.id in self._deactivated
            label = f"{seq.id}  ({len(seq.seq)} bp)"
            if is_deactivated:
                label += "  [inactive]"

            item = QListWidgetItem(label)
            item.setData(Qt.UserRole, seq.id)

            if is_deactivated:
                # Italic + gray for deactivated oligos
                font = item.font()
                font.setItalic(True)
                item.setFont(font)
                item.setForeground(QColor(140, 140, 140))
            else:
                # Apply risk color only for active oligos
                risk = self._risk_map.get(seq.id, "NONE")
                color = self.RISK_COLORS.get(risk)
                if color:
                    item.setBackground(color)

            self.list_widget.addItem(item)
            self._display_indices.append(orig_idx)
