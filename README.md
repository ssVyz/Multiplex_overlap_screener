# Multiplex Assay Overlap Screener

A desktop tool for checking 3' end overlaps between oligonucleotides in multiplex assay designs. It performs pairwise comparison of oligo sequences within user-defined mixes and flags potential primer-dimer interactions by overlap length, mismatch count, and risk level.

## Requirements

- Python 3.12+
- PySide6
- Biopython

Install dependencies:

```
pip install PySide6 biopython
```

## Running

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

- **Interaction partners** — a table of all pairwise 3' overlaps involving the selected oligo, with overlap length, mismatch count, and risk level.
- **Overlap visualization** — an ASCII alignment of the two primers showing the 3' overlap region. Click a row in the interaction table to view it.

Risk levels (HIGH, MEDIUM, LOW) are determined by configurable thresholds in **Settings**. Oligos in the tree are highlighted by their highest risk interaction.

### Projects

- **File > Save Project / Save As** — saves the full state (all oligos, mixes, active/inactive status) to a JSON file.
- **File > Open Project** — loads a previously saved project.
- **File > Export All as FASTA** — exports every oligo in the project regardless of mix assignment.

### Settings

Accessible from the menu bar. Configurable parameters:

| Parameter | Description | Default |
|---|---|---|
| Minimum overlap length | Shortest overlap to check | 3 |
| Maximum overlap length | Longest overlap to check | 10 |
| Maximum mismatches | Max allowed mismatches in an overlap | 1 |
| Consider ambiguity | Treat IUPAC ambiguous bases as matches if any variation could pair | Off |
| High risk thresholds | Min overlap and max mismatches for HIGH | 4 / 0 |
| Medium risk thresholds | Min overlap and max mismatches for MEDIUM | 2 / 1 |

Settings are stored in `settings.json` in the project root.

## Project structure

```
main.py                  Entry point
settings.json            Persisted analysis settings
core/
    models.py            Data model: Oligo, Mix, Project (serialization, FASTA export)
    analysis.py          Pairwise 3' overlap analysis, risk classification, visualization
    settings.py          Settings load/save
ui/
    main_window.py       Main window layout (4-panel) and coordination
    oligo_tree.py        Tree widget with mixes/oligos, drag-drop, context menus, import dialog
    oligo_preview.py     Selected oligo details (name, length, GC%, sequence)
    detail_panel.py      Interaction table and overlap ASCII visualization
    settings_dialog.py   Settings editor dialog
```

## How the analysis works

For each pair of oligos in a mix (including self-pairs), the tool checks if the 3' end of one oligo is complementary to the 3' end of the other. It does this for every overlap length between the configured min and max:

1. Extract the last N bases from oligo A.
2. Extract the last N bases from oligo B and reverse-complement them.
3. Count mismatches between the two fragments (optionally treating IUPAC ambiguity codes as matches).
4. If mismatches are within the configured limit, record the interaction with its risk level.

Results are sorted by overlap length (longest first), then by mismatch count. Risk levels are assigned based on the overlap/mismatch thresholds in settings.

IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are supported. When ambiguity consideration is enabled, two bases are considered a match if any of their possible expansions overlap.
