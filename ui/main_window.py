from PySide6.QtWidgets import QMainWindow, QSplitter, QMenuBar
from PySide6.QtCore import Qt

from core.settings import load_settings, save_settings
from core.analysis import analyze_mix, get_interactions_for_oligo, get_max_risk_for_oligo
from ui.oligo_panel import OligoPanel
from ui.detail_panel import DetailPanel
from ui.settings_dialog import SettingsDialog


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Multiplex Assay Overlap Screener")
        self.resize(1100, 700)

        self._settings = load_settings()
        self._analysis_results = []

        self._build_menu()
        self._build_ui()
        self._connect_signals()

    def _build_menu(self):
        menu_bar = self.menuBar()
        settings_action = menu_bar.addAction("Settings")
        settings_action.triggered.connect(self._open_settings)

    def _build_ui(self):
        splitter = QSplitter(Qt.Horizontal)

        self.oligo_panel = OligoPanel()
        self.oligo_panel.set_sequences([])
        splitter.addWidget(self.oligo_panel)

        self.detail_panel = DetailPanel()
        splitter.addWidget(self.detail_panel)

        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)

        self.setCentralWidget(splitter)

    def _connect_signals(self):
        self.oligo_panel.oligos_changed.connect(self._run_analysis)
        self.oligo_panel.selection_changed.connect(self._on_oligo_selected)

    def _open_settings(self):
        dlg = SettingsDialog(self._settings, self)
        if dlg.exec():
            self._settings = dlg.get_settings()
            save_settings(self._settings)
            self._run_analysis()

    def _run_analysis(self):
        sequences = self.oligo_panel.get_sequences()
        self._analysis_results = analyze_mix(sequences, self._settings)

        # Update risk highlights on the oligo list
        risk_map = {}
        for seq in sequences:
            risk_map[seq.id] = get_max_risk_for_oligo(seq.id, self._analysis_results)
        self.oligo_panel.update_risk_highlights(risk_map)

        # Refresh detail panel if an oligo is selected
        current = self.oligo_panel.list_widget.currentItem()
        if current:
            self._on_oligo_selected(current.data(1))
        else:
            self.detail_panel.clear()

    def _on_oligo_selected(self, oligo_id):
        if not oligo_id:
            self.detail_panel.clear()
            return
        interactions = get_interactions_for_oligo(oligo_id, self._analysis_results)
        self.detail_panel.show_interactions(oligo_id, interactions)
