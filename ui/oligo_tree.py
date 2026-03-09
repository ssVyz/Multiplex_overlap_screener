from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTreeWidget, QTreeWidgetItem,
    QDialog, QPlainTextEdit, QLabel, QFileDialog, QMessageBox, QMenu,
    QInputDialog, QLineEdit,
)
from PySide6.QtCore import Signal, Qt
from PySide6.QtGui import QColor, QFont

from core.analysis import parse_fasta_text
from core.models import Oligo, Mix, Project

# Custom data roles
ROLE_TYPE = Qt.ItemDataRole.UserRole          # "unallocated", "mix", "oligo"
ROLE_ID = Qt.ItemDataRole.UserRole + 1        # mix_id or oligo_id


class EditOligoDialog(QDialog):
    """Dialog for editing an oligo's name and sequence."""

    def __init__(self, oligo, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Edit Oligo")
        self.setMinimumSize(400, 300)
        self._oligo = oligo
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel("Name:"))
        self.name_edit = QLineEdit(self._oligo.name)
        layout.addWidget(self.name_edit)

        layout.addWidget(QLabel("Sequence:"))
        self.seq_edit = QPlainTextEdit(self._oligo.sequence)
        self.seq_edit.setFont(QFont("Courier New", 10))
        layout.addWidget(self.seq_edit)

        btn_row = QHBoxLayout()
        btn_row.addStretch()
        ok_btn = QPushButton("OK")
        ok_btn.setDefault(True)
        ok_btn.clicked.connect(self._on_ok)
        btn_row.addWidget(ok_btn)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

    def _on_ok(self):
        name = self.name_edit.text().strip()
        seq = self.seq_edit.toPlainText().strip().replace("\n", "").replace(" ", "")
        if not name:
            QMessageBox.warning(self, "Invalid", "Name cannot be empty.")
            return
        if not seq:
            QMessageBox.warning(self, "Invalid", "Sequence cannot be empty.")
            return
        self.result_name = name
        self.result_sequence = seq.upper()
        self.accept()


class ImportOligoDialog(QDialog):
    """Dialog for importing oligos via FASTA text or file."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Import Oligos")
        self.setMinimumSize(500, 400)
        self.records = []
        self._build_ui()

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
                    self.text_edit.setPlainText(f.read())
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


class OligoTree(QWidget):
    """Top-left panel: tree view with mixes and oligos, drag-drop, context menus."""

    oligo_selected = Signal(str)    # oligo UUID
    project_changed = Signal()      # data changed, re-run analysis

    RISK_COLORS = {
        "HIGH": QColor(255, 120, 120, 90),
        "MEDIUM": QColor(255, 200, 80, 90),
    }

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self._project = project
        self._risk_map = {}
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        btn_row = QHBoxLayout()
        import_btn = QPushButton("Import Oligos")
        import_btn.clicked.connect(self._on_import)
        btn_row.addWidget(import_btn)

        new_mix_btn = QPushButton("New Mix")
        new_mix_btn.clicked.connect(self._on_new_mix)
        btn_row.addWidget(new_mix_btn)
        layout.addLayout(btn_row)

        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.tree.setDragEnabled(True)
        self.tree.setAcceptDrops(True)
        self.tree.setDragDropMode(QTreeWidget.DragDrop)
        self.tree.setDefaultDropAction(Qt.MoveAction)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._show_context_menu)
        self.tree.currentItemChanged.connect(self._on_current_changed)

        # Override drop to update model
        self.tree.dropEvent = self._handle_drop

        layout.addWidget(self.tree)

    # ---- public API ----

    def set_project(self, project: Project):
        self._project = project
        self._risk_map = {}
        self.refresh()

    def update_risk_highlights(self, risk_map: dict):
        self._risk_map = dict(risk_map)
        self.refresh()

    def refresh(self):
        """Rebuild the tree from the project model, preserving selection."""
        # Save current selection
        current = self.tree.currentItem()
        selected_id = current.data(0, ROLE_ID) if current else None
        selected_type = current.data(0, ROLE_TYPE) if current else None

        self.tree.clear()

        # Unallocated section
        unalloc_item = QTreeWidgetItem(self.tree, ["Unallocated"])
        unalloc_item.setData(0, ROLE_TYPE, "unallocated")
        unalloc_item.setData(0, ROLE_ID, None)
        unalloc_item.setFlags(
            Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDropEnabled
        )
        font = unalloc_item.font(0)
        font.setBold(True)
        unalloc_item.setFont(0, font)

        for oligo in self._project.get_unallocated_oligos():
            self._add_oligo_item(unalloc_item, oligo)
        unalloc_item.setExpanded(True)

        # Mix sections
        for mix in self._project.mixes:
            mix_item = QTreeWidgetItem(self.tree, [f"{mix.name} ({len(self._project.get_oligos_in_mix(mix.id))})"])
            mix_item.setData(0, ROLE_TYPE, "mix")
            mix_item.setData(0, ROLE_ID, mix.id)
            mix_item.setFlags(
                Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDropEnabled
            )
            font = mix_item.font(0)
            font.setBold(True)
            mix_item.setFont(0, font)

            for oligo in self._project.get_oligos_in_mix(mix.id):
                self._add_oligo_item(mix_item, oligo)
            mix_item.setExpanded(True)

        # Restore selection
        self._restore_selection(selected_type, selected_id)

    def get_selected_oligo_ids(self) -> list[str]:
        """Return UUIDs of all selected oligo items."""
        ids = []
        for item in self.tree.selectedItems():
            if item.data(0, ROLE_TYPE) == "oligo":
                ids.append(item.data(0, ROLE_ID))
        return ids

    # ---- tree helpers ----

    def _add_oligo_item(self, parent, oligo):
        label = f"{oligo.name}  ({len(oligo.sequence)} bp)"
        if not oligo.active:
            label += "  [inactive]"
        item = QTreeWidgetItem(parent, [label])
        item.setData(0, ROLE_TYPE, "oligo")
        item.setData(0, ROLE_ID, oligo.id)
        item.setFlags(
            Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsDragEnabled
        )

        if not oligo.active:
            f = item.font(0)
            f.setItalic(True)
            item.setFont(0, f)
            item.setForeground(0, QColor(140, 140, 140))
        else:
            risk = self._risk_map.get(oligo.id, "NONE")
            color = self.RISK_COLORS.get(risk)
            if color:
                item.setBackground(0, color)

    def _restore_selection(self, sel_type, sel_id):
        if sel_type is None:
            return
        it = QTreeWidgetItem()
        iterator = self._iterate_all_items()
        for item in iterator:
            if item.data(0, ROLE_TYPE) == sel_type and item.data(0, ROLE_ID) == sel_id:
                self.tree.setCurrentItem(item)
                return

    def _iterate_all_items(self):
        """Yield all items in the tree."""
        for i in range(self.tree.topLevelItemCount()):
            top = self.tree.topLevelItem(i)
            yield top
            for j in range(top.childCount()):
                yield top.child(j)

    # ---- drag and drop ----

    def _handle_drop(self, event):
        """Handle drop event: move oligos to the target mix/unallocated."""
        target_item = self.tree.itemAt(event.position().toPoint())
        if not target_item:
            event.ignore()
            return

        target_type = target_item.data(0, ROLE_TYPE)

        # If dropped on an oligo, use its parent
        if target_type == "oligo":
            target_item = target_item.parent()
            if not target_item:
                event.ignore()
                return
            target_type = target_item.data(0, ROLE_TYPE)

        if target_type == "unallocated":
            target_mix_id = None
        elif target_type == "mix":
            target_mix_id = target_item.data(0, ROLE_ID)
        else:
            event.ignore()
            return

        moved = False
        for item in self.tree.selectedItems():
            if item.data(0, ROLE_TYPE) == "oligo":
                oligo_id = item.data(0, ROLE_ID)
                oligo = self._project.get_oligo_by_id(oligo_id)
                if oligo and oligo.mix_id != target_mix_id:
                    oligo.mix_id = target_mix_id
                    moved = True

        if moved:
            self.refresh()
            self.project_changed.emit()

        event.accept()

    # ---- context menus ----

    def _show_context_menu(self, pos):
        item = self.tree.itemAt(pos)
        menu = QMenu(self)

        if item is None:
            menu.addAction("Import Oligos", self._on_import)
            menu.addAction("New Mix", self._on_new_mix)
            menu.exec(self.tree.viewport().mapToGlobal(pos))
            return

        item_type = item.data(0, ROLE_TYPE)

        if item_type == "oligo":
            selected_oligo_ids = self.get_selected_oligo_ids()
            multi = len(selected_oligo_ids) > 1

            if not multi:
                menu.addAction("Edit Oligo", lambda: self._on_edit_oligo(item))
                menu.addAction("Duplicate", lambda: self._on_duplicate_oligo(item))
                menu.addAction("Rename", lambda: self._on_rename_oligo(item))
                menu.addSeparator()

            # Activate / Deactivate
            oligo = self._project.get_oligo_by_id(item.data(0, ROLE_ID))
            if multi:
                menu.addAction("Activate Selected", lambda: self._on_set_active(True))
                menu.addAction("Deactivate Selected", lambda: self._on_set_active(False))
            else:
                if oligo and oligo.active:
                    menu.addAction("Deactivate", lambda: self._on_toggle_active(item))
                else:
                    menu.addAction("Activate", lambda: self._on_toggle_active(item))

            menu.addSeparator()

            # Move to mix submenu
            move_menu = menu.addMenu("Move to...")
            move_menu.addAction("Unallocated", lambda: self._on_move_selected(None))
            for mix in self._project.mixes:
                move_menu.addAction(mix.name, lambda mid=mix.id: self._on_move_selected(mid))

            menu.addSeparator()
            if multi:
                menu.addAction("New Mix from Selection", self._on_new_mix_from_selection)
                menu.addAction(f"Delete {len(selected_oligo_ids)} Oligos", self._on_delete_selected_oligos)
            else:
                menu.addAction("Delete", lambda: self._on_delete_oligo(item))

        elif item_type == "mix":
            menu.addAction("Rename Mix", lambda: self._on_rename_mix(item))
            menu.addAction("Export Mix as FASTA", lambda: self._on_export_mix(item))
            menu.addSeparator()
            menu.addAction("Delete Mix", lambda: self._on_delete_mix(item))

        elif item_type == "unallocated":
            menu.addAction("Import Oligos", self._on_import)
            menu.addAction("New Mix", self._on_new_mix)

        menu.exec(self.tree.viewport().mapToGlobal(pos))

    # ---- actions ----

    def _on_import(self):
        dlg = ImportOligoDialog(self)
        if dlg.exec() == QDialog.Accepted and dlg.records:
            for rec in dlg.records:
                oligo = Oligo.from_seq_record(rec, mix_id=None)
                self._project.add_oligo(oligo)
            self.refresh()
            self.project_changed.emit()

    def _on_new_mix(self):
        name, ok = QInputDialog.getText(self, "New Mix", "Mix name:")
        if ok and name.strip():
            self._project.add_mix(Mix(name=name.strip()))
            self.refresh()
            self.project_changed.emit()

    def _on_new_mix_from_selection(self):
        oligo_ids = self.get_selected_oligo_ids()
        if not oligo_ids:
            return
        name, ok = QInputDialog.getText(self, "New Mix from Selection", "Mix name:")
        if ok and name.strip():
            mix = Mix(name=name.strip())
            self._project.add_mix(mix)
            for oid in oligo_ids:
                oligo = self._project.get_oligo_by_id(oid)
                if oligo:
                    oligo.mix_id = mix.id
            self.refresh()
            self.project_changed.emit()

    def _on_edit_oligo(self, item):
        oligo_id = item.data(0, ROLE_ID)
        oligo = self._project.get_oligo_by_id(oligo_id)
        if not oligo:
            return
        dlg = EditOligoDialog(oligo, self)
        if dlg.exec() == QDialog.Accepted:
            oligo.name = dlg.result_name
            oligo.sequence = dlg.result_sequence
            self.refresh()
            self.project_changed.emit()

    def _on_duplicate_oligo(self, item):
        oligo_id = item.data(0, ROLE_ID)
        dup = self._project.duplicate_oligo(oligo_id)
        if dup:
            self.refresh()
            self.project_changed.emit()

    def _on_rename_oligo(self, item):
        oligo_id = item.data(0, ROLE_ID)
        oligo = self._project.get_oligo_by_id(oligo_id)
        if not oligo:
            return
        name, ok = QInputDialog.getText(self, "Rename Oligo", "New name:", text=oligo.name)
        if ok and name.strip():
            oligo.name = name.strip()
            self.refresh()
            self.project_changed.emit()

    def _on_toggle_active(self, item):
        oligo_id = item.data(0, ROLE_ID)
        oligo = self._project.get_oligo_by_id(oligo_id)
        if oligo:
            oligo.active = not oligo.active
            self.refresh()
            self.project_changed.emit()

    def _on_set_active(self, active: bool):
        for oid in self.get_selected_oligo_ids():
            oligo = self._project.get_oligo_by_id(oid)
            if oligo:
                oligo.active = active
        self.refresh()
        self.project_changed.emit()

    def _on_move_selected(self, target_mix_id):
        for oid in self.get_selected_oligo_ids():
            oligo = self._project.get_oligo_by_id(oid)
            if oligo:
                oligo.mix_id = target_mix_id
        self.refresh()
        self.project_changed.emit()

    def _on_delete_oligo(self, item):
        oligo_id = item.data(0, ROLE_ID)
        oligo = self._project.get_oligo_by_id(oligo_id)
        if not oligo:
            return
        reply = QMessageBox.question(
            self, "Delete Oligo", f"Delete '{oligo.name}'?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply == QMessageBox.Yes:
            self._project.remove_oligo(oligo_id)
            self.refresh()
            self.project_changed.emit()

    def _on_delete_selected_oligos(self):
        ids = self.get_selected_oligo_ids()
        if not ids:
            return
        reply = QMessageBox.question(
            self, "Delete Oligos", f"Delete {len(ids)} selected oligo(s)?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply == QMessageBox.Yes:
            for oid in ids:
                self._project.remove_oligo(oid)
            self.refresh()
            self.project_changed.emit()

    def _on_rename_mix(self, item):
        mix_id = item.data(0, ROLE_ID)
        mix = self._project.get_mix_by_id(mix_id)
        if not mix:
            return
        name, ok = QInputDialog.getText(self, "Rename Mix", "New name:", text=mix.name)
        if ok and name.strip():
            mix.name = name.strip()
            self.refresh()

    def _on_delete_mix(self, item):
        mix_id = item.data(0, ROLE_ID)
        mix = self._project.get_mix_by_id(mix_id)
        if not mix:
            return
        reply = QMessageBox.question(
            self, "Delete Mix",
            f"Delete mix '{mix.name}'?\nOligos will be moved to Unallocated.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply == QMessageBox.Yes:
            self._project.remove_mix(mix_id)
            self.refresh()
            self.project_changed.emit()

    def _on_export_mix(self, item):
        mix_id = item.data(0, ROLE_ID)
        mix = self._project.get_mix_by_id(mix_id)
        if not mix:
            return
        oligos = self._project.get_oligos_in_mix(mix_id)
        if not oligos:
            QMessageBox.information(self, "Empty Mix", "This mix has no oligos to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Mix as FASTA", f"{mix.name}.fasta",
            "FASTA files (*.fasta);;All files (*)"
        )
        if path:
            try:
                count = self._project.export_fasta(path, mix_id)
                QMessageBox.information(self, "Exported", f"Saved {count} sequences to:\n{path}")
            except Exception as e:
                QMessageBox.critical(self, "Export error", str(e))

    def _on_current_changed(self, current, _previous):
        if current is None:
            self.oligo_selected.emit("")
            return
        if current.data(0, ROLE_TYPE) == "oligo":
            self.oligo_selected.emit(current.data(0, ROLE_ID))
        else:
            self.oligo_selected.emit("")
