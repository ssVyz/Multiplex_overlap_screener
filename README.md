# Multiplex Assay Overlap Screener

A desktop tool for screening oligonucleotides in multiplex assay designs for primer-dimer interactions. It performs pairwise comparison of oligo sequences within user-defined mixes using either of two screening modes:

- **3' end overlap** — flags interactions by overlap length and mismatch count.
- **Thermodynamic (ΔG)** — scores every binding register by nearest-neighbor free energy, reporting the strongest single interaction, the total (ensemble) interaction, and 3' anchored stability.

The mode is selected in Settings and applies to the whole project.

## Requirements

- Python 3.12+
- PySide6
- Biopython

## Setup

This project uses [uv](https://docs.astral.sh/uv/). Dependencies and the Python version are declared in `pyproject.toml` / `.python-version`.

```
uv sync
```

This creates a `.venv` and installs everything pinned in `uv.lock`. uv will fetch Python 3.12 automatically if it isn't already installed.

## Running

```
uv run python main.py
```

Or, with the venv activated (`.venv\Scripts\activate` on Windows, `source .venv/bin/activate` elsewhere):

```
python main.py
```

## Usage

### Importing oligos

Use **Import Oligos** (button or right-click menu) to add oligonucleotides. Input is FASTA format, either pasted as text or loaded from a `.fasta` file. Imported oligos are placed in the **Unallocated** section and are not analyzed until assigned to a mix.

### Mixes

Mixes are flat groups (no nesting). Each mix acts as an independent analysis unit — overlap checks only run between oligos within the same mix.

- Create a mix via the **New Mix** button, right-click menu, or by selecting multiple oligos and choosing **New Mix from Selection**.
- Move oligos between mixes by drag-and-drop, or via right-click > **Move to...**.
- An oligo belongs to exactly one mix (or is unallocated). Moving it removes it from its previous location.

### Oligo operations (right-click)

- **Edit Oligo** — change name and sequence
- **Duplicate** — creates a copy in the same location with a `_N` suffix
- **Rename** — change the display name
- **Activate / Deactivate** — inactive oligos are excluded from analysis but kept in the project
- **Delete** — removes the oligo entirely

### Mix operations (right-click)

- **Rename Mix**
- **Export Mix as FASTA**
- **Delete Mix** — oligos are moved to Unallocated, not deleted

### Analysis

When an oligo inside a mix is selected, the right-side panels show:

- **Interaction partners** — a table of all pairwise interactions involving the selected oligo. Columns depend on the active mode: overlap length and mismatch count in overlap mode, the three ΔG values in thermodynamic mode.
- **Visualization** — an ASCII alignment of the two primers showing the interacting region. Click a row in the interaction table to view it.

Risk levels (HIGH, MEDIUM, LOW) are determined by configurable thresholds in **Settings**. Oligos in the tree are highlighted by their highest risk interaction. The title bar shows which mode is active.

### Projects

- **File > Save Project / Save As** — saves the full state (all oligos, mixes, active/inactive status) to a JSON file.
- **File > Open Project** — loads a previously saved project.
- **File > Export All as FASTA** — exports every oligo in the project regardless of mix assignment.

### Settings

Accessible from the menu bar. Configurable parameters:

Shared:

| Parameter | Description | Default |
|---|---|---|
| Mode | `3' end overlap` or `Thermodynamic (ΔG)` | 3' end overlap |
| Consider ambiguity | Treat IUPAC ambiguous bases as matches if any variation could pair | Off |
| Na⁺ concentration | Salt concentration, used for Tm and ΔG salt correction | 50 mM |

Overlap mode:

| Parameter | Description | Default |
|---|---|---|
| Minimum overlap length | Shortest overlap to check | 3 |
| Maximum overlap length | Longest overlap to check | 10 |
| Maximum mismatches | Max allowed mismatches in an overlap | 1 |
| High risk thresholds | Min overlap and max mismatches for HIGH | 4 / 0 |
| Medium risk thresholds | Min overlap and max mismatches for MEDIUM | 2 / 1 |

Thermodynamic mode (all ΔG in kcal/mol; more negative = more stable):

| Parameter | Description | Default |
|---|---|---|
| Temperature | Temperature at which ΔG is evaluated | 37 °C |
| Report thresholds | A pair is listed if *either* its strongest or its total ΔG is at or below these | -2.0 / -3.0 |
| High risk thresholds | 3' anchored ΔG, or strongest ΔG, at or below these | -6.0 / -10.0 |
| Medium risk thresholds | 3' anchored ΔG, or strongest ΔG, at or below these | -4.0 / -7.0 |

Settings are stored in `settings.json` in the project root. Keys absent from an
existing file fall back to defaults, so older settings files keep working.

## Project structure

```
main.py                  Entry point
settings.json            Persisted analysis settings
core/
    models.py            Data model: Oligo, Mix, Project (serialization, FASTA export)
    analysis.py          Mode dispatch, overlap analysis, risk classification, visualization
    thermo.py            Nearest-neighbor free-energy scoring of duplexes
    sequtils.py          IUPAC codes and base comparison shared by the above
    settings.py          Settings load/save
ui/
    main_window.py       Main window layout (4-panel) and coordination
    oligo_tree.py        Tree widget with mixes/oligos, drag-drop, context menus, import dialog
    oligo_preview.py     Selected oligo details (name, length, GC%, sequence)
    detail_panel.py      Interaction table and overlap ASCII visualization
    settings_dialog.py   Settings editor dialog
```

## How the analysis works

### 3' end overlap mode

For each pair of oligos in a mix (including self-pairs), the tool checks if the 3' end of one oligo is complementary to the 3' end of the other. It does this for every overlap length between the configured min and max:

1. Extract the last N bases from oligo A.
2. Extract the last N bases from oligo B and reverse-complement them.
3. Count mismatches between the two fragments (optionally treating IUPAC ambiguity codes as matches).
4. If mismatches are within the configured limit, record the interaction with its risk level.

Results are sorted by overlap length (longest first), then by mismatch count. Risk levels are assigned based on the overlap/mismatch thresholds in settings.

IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are supported. When ambiguity consideration is enabled, two bases are considered a match if any of their possible expansions overlap.

### Thermodynamic (ΔG) mode

This mode scores duplex stability rather than counting matches. For each pair it
scans every gapless register (every relative offset of the two strands) and finds
the most stable contiguous helix in each. A helix must begin and end on a
Watson-Crick pair but may carry mismatches internally.

Energies come from the SantaLucia/Allawi nearest-neighbor parameters shipped with
Biopython — stacking, internal mismatches, terminal mismatches and dangling ends
— plus helix initiation penalties, the SantaLucia (1998) salt correction from the
Na⁺ setting, and the symmetry correction for self-complementary structures. No
external thermodynamics dependency is required.

Three numbers are reported per pair, one row per pair rather than one per overlap
length:

| Column | Meaning |
|---|---|
| **ΔG 3'** | Stability of the most stable helix that actually pairs an oligo's 3' terminal base, taking the worse of the two directions. Structures that leave both 3' ends unpaired cannot be extended by polymerase, so this is usually the best predictor of primer-dimer. |
| **ΔG worst** | The single most stable structure found — the classic worst-case dimer. |
| **ΔG total** | Ensemble free energy over all registers, `-RT ln Σ exp(-ΔGᵢ/RT)`. Always at least as negative as ΔG worst; the gap between them indicates how many competing structures exist. |

ΔG total is a Boltzmann sum, not an arithmetic one. Adding the registers' energies
would scale with how many registers happen to be scanned, so long oligos would
score worse for purely geometric reasons and many weak junk interactions would
outweigh one real one. The ensemble sum is dominated by the strongest structure
and weak interactions contribute almost nothing, while several comparable
structures do pull it down.

Ambiguity is handled by scoring a degenerate position as its most stable possible
pairing, i.e. worst case. With ambiguity consideration off, degenerate positions
do not pair at all.

**Limitations.** Bulges, internal loops and multi-helix structures within a single
register are not modelled, so structures that depend on a bulge score slightly
less negative than primer3 or IDT OligoAnalyzer would report. Hairpins
(intramolecular folding) are not evaluated. Absolute values will not match primer3
exactly; validated against the published SantaLucia duplex CGTTGA the model gives
-5.41 vs -5.35 kcal/mol, the difference being ΔH/ΔS rounding in the tables.
