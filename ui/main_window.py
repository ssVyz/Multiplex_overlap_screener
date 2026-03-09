from PySide6.QtWidgets import (
    QMainWindow, QSplitter, QFileDialog, QMessageBox,
)
from PySide6.QtCore import Qt

from core.settings import load_settings, save_settings
from core.analysis import analyze_mix, get_interactions_for_oligo, get_max_risk_for_oligo
from core.models import Project
from ui.oligo_tree import OligoTree
from ui.oligo_preview import OligoPreview
from ui.detail_panel import DetailPanel
from ui.settings_dialog import SettingsDialog


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Multiplex Assay Overlap Screener")
        self.resize(1200, 800)

        self._settings = load_settings()
        self._project = Project()
        self._project_path = None
        self._analysis_cache = {}  # mix_id -> results list

        self._build_menu()
        self._build_ui()
        self._connect_signals()

    def _build_menu(self):
        menu_bar = self.menuBar()

        file_menu = menu_bar.addMenu("File")
        file_menu.addAction("New Project", self._on_new_project)
        file_menu.addAction("Open Project...", self._on_open_project)
        file_menu.addAction("Save Project", self._on_save_project)
        file_menu.addAction("Save Project As...", self._on_save_project_as)
        file_menu.addSeparator()
        file_menu.addAction("Export All as FASTA...", self._on_export_all)

        settings_action = menu_bar.addAction("Settings")
        settings_action.triggered.connect(self._open_settings)

    def _build_ui(self):
        # Main horizontal splitter: left | right
        main_splitter = QSplitter(Qt.Horizontal)

        # Left vertical splitter: tree (top) | preview (bottom)
        left_splitter = QSplitter(Qt.Vertical)
        self.oligo_tree = OligoTree(self._project)
        left_splitter.addWidget(self.oligo_tree)

        self.oligo_preview = OligoPreview()
        left_splitter.addWidget(self.oligo_preview)

        left_splitter.setStretchFactor(0, 3)
        left_splitter.setStretchFactor(1, 1)
        main_splitter.addWidget(left_splitter)

        # Right side: detail panel (interactions table + visualization)
        self.detail_panel = DetailPanel()
        main_splitter.addWidget(self.detail_panel)

        main_splitter.setStretchFactor(0, 1)
        main_splitter.setStretchFactor(1, 2)

        self.setCentralWidget(main_splitter)

    def _connect_signals(self):
        self.oligo_tree.oligo_selected.connect(self._on_oligo_selected)
        self.oligo_tree.project_changed.connect(self._on_project_changed)

    # ---- settings ----

    def _open_settings(self):
        dlg = SettingsDialog(self._settings, self)
        if dlg.exec():
            self._settings = dlg.get_settings()
            save_settings(self._settings)
            self._analysis_cache.clear()
            self._refresh_current_view()

    # ---- project management ----

    def _on_new_project(self):
        if self._project.oligos or self._project.mixes:
            reply = QMessageBox.question(
                self, "New Project",
                "Discard current project and start fresh?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return
        self._project = Project()
        self._project_path = None
        self._analysis_cache.clear()
        self.oligo_tree.set_project(self._project)
        self.oligo_preview.clear()
        self.detail_panel.clear()
        self.setWindowTitle("Multiplex Assay Overlap Screener")

    def _on_open_project(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Project", "",
            "Project files (*.json);;All files (*)"
        )
        if path:
            try:
                self._project = Project.load(path)
                self._project_path = path
                self._analysis_cache.clear()
                self.oligo_tree.set_project(self._project)
                self.oligo_preview.clear()
                self.detail_panel.clear()
                self.setWindowTitle(f"Overlap Screener — {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load project:\n{e}")

    def _on_save_project(self):
        if hasattr(self, '_project_path') and self._project_path:
            try:
                self._project.save(self._project_path)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save:\n{e}")
        else:
            self._on_save_project_as()

    def _on_save_project_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Project As", "project.json",
            "Project files (*.json);;All files (*)"
        )
        if path:
            try:
                self._project.save(path)
                self._project_path = path
                self.setWindowTitle(f"Overlap Screener — {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save:\n{e}")

    def _on_export_all(self):
        if not self._project.oligos:
            QMessageBox.information(self, "Nothing to export", "No oligos in the project.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export All as FASTA", "all_oligos.fasta",
            "FASTA files (*.fasta);;All files (*)"
        )
        if path:
            try:
                count = self._project.export_fasta(path, "ALL")
                QMessageBox.information(self, "Exported", f"Saved {count} sequences to:\n{path}")
            except Exception as e:
                QMessageBox.critical(self, "Export error", str(e))

    # ---- analysis ----

    def _on_project_changed(self):
        """Called when oligos are added/removed/moved/edited."""
        self._analysis_cache.clear()
        self._refresh_current_view()

    def _run_analysis_for_mix(self, mix_id):
        """Run (or return cached) analysis for a specific mix."""
        if mix_id in self._analysis_cache:
            return self._analysis_cache[mix_id]
        active_oligos = self._project.get_active_oligos_in_mix(mix_id)
        results = analyze_mix(active_oligos, self._settings)
        self._analysis_cache[mix_id] = results
        return results

    def _refresh_current_view(self):
        """Re-run analysis and update risk highlights for all mixes."""
        # Build combined risk map across all mixes
        risk_map = {}
        for mix in self._project.mixes:
            results = self._run_analysis_for_mix(mix.id)
            for oligo in self._project.get_active_oligos_in_mix(mix.id):
                risk_map[oligo.id] = get_max_risk_for_oligo(oligo.id, results)

        self.oligo_tree.update_risk_highlights(risk_map)

        # Refresh detail panel for currently selected oligo
        current = self.oligo_tree.tree.currentItem()
        if current and current.data(0, Qt.ItemDataRole.UserRole) == "oligo":
            oligo_id = current.data(0, Qt.ItemDataRole.UserRole + 1)
            self._show_oligo_details(oligo_id)
        else:
            self.detail_panel.clear()

    def _on_oligo_selected(self, oligo_id):
        if not oligo_id:
            self.oligo_preview.clear()
            self.detail_panel.clear()
            return
        self._show_oligo_details(oligo_id)

    def _show_oligo_details(self, oligo_id):
        oligo = self._project.get_oligo_by_id(oligo_id)
        if not oligo:
            self.oligo_preview.clear()
            self.detail_panel.clear()
            return

        # Preview panel
        mix_name = None
        if oligo.mix_id:
            mix = self._project.get_mix_by_id(oligo.mix_id)
            mix_name = mix.name if mix else None
        self.oligo_preview.show_oligo(oligo, mix_name)

        # Detail panel — only show interactions if oligo is in a mix
        if oligo.mix_id:
            results = self._run_analysis_for_mix(oligo.mix_id)
            interactions = get_interactions_for_oligo(oligo.id, results)
            self.detail_panel.show_interactions(oligo.id, oligo.name, interactions)
        else:
            self.detail_panel.clear()
