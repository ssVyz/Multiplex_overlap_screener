from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

from core.sequtils import (  # noqa: F401
    IUPAC_CODES, bases_could_match, bases_could_pair, duplex_bases_pair,
)
from core.thermo import model_from_settings

# Screening modes. OVERLAP is the original 3' end overlap scan; THERMO scores
# every gapless register by nearest-neighbor free energy (see core/thermo.py).
MODE_OVERLAP = "overlap"
MODE_THERMO = "thermo"


def count_mismatches(seq1, seq2, consider_ambiguity=False):
    """Count mismatches between two equal-length sequences."""
    if consider_ambiguity:
        return sum(1 for a, b in zip(seq1, seq2) if not bases_could_match(a, b))
    else:
        return sum(1 for a, b in zip(seq1, seq2) if a != b)


def get_risk_level(overlap_length, mismatches, settings):
    """Determine risk level using configurable thresholds from settings."""
    high_min_ol = settings.get("high_risk_min_overlap", 4)
    high_max_mm = settings.get("high_risk_max_mismatches", 0)
    med_min_ol = settings.get("medium_risk_min_overlap", 2)
    med_max_mm = settings.get("medium_risk_max_mismatches", 1)

    if overlap_length >= high_min_ol and mismatches <= high_max_mm:
        return "HIGH"
    elif overlap_length >= med_min_ol and mismatches <= med_max_mm:
        return "MEDIUM"
    else:
        return "LOW"


def get_thermo_risk_level(dg_3prime, dg_min, settings):
    """Determine risk level from free energies.

    3' anchored stability is checked first because only 3' paired structures can
    be extended by polymerase; a very stable internal duplex still counts, via
    the separate dg_min thresholds.
    """
    high_3p = settings.get("high_risk_dg_3prime", -6.0)
    high_any = settings.get("high_risk_dg_any", -10.0)
    med_3p = settings.get("medium_risk_dg_3prime", -4.0)
    med_any = settings.get("medium_risk_dg_any", -7.0)

    if dg_3prime <= high_3p or dg_min <= high_any:
        return "HIGH"
    elif dg_3prime <= med_3p or dg_min <= med_any:
        return "MEDIUM"
    else:
        return "LOW"


def get_last_n_bases(sequence, n):
    """Get last n bases from 3' end of a sequence string."""
    return sequence[-n:].upper()


def get_last_n_bases_rc(sequence, n):
    """Get last n bases from 3' end and return reverse complement."""
    trimmed = Seq(sequence[-n:].upper())
    return str(trimmed.reverse_complement())


def visualize_duplex(oligo1, oligo2, offset, helix_start, helix_end,
                     consider_ambiguity=False):
    """ASCII view of two oligos held in duplex at a given register.

    oligo2 is drawn reversed so it reads 3'->5' left to right. `offset` places
    it against oligo1: base k of the reversed oligo2 sits under base k + offset
    of oligo1. `helix_start`/`helix_end` mark the paired span in oligo1
    coordinates, and only that span is annotated with pairing bars.
    """
    seq1 = oligo1.sequence
    seq2_rev = oligo2.sequence[::-1]

    pad = max(0, -offset)  # shift right when oligo2 overhangs oligo1's 5' end

    lines = [f"5'-{' ' * pad}{seq1}-3'  ({oligo1.name})"]

    match_line = "   " + " " * (pad + helix_start)
    for col in range(helix_start, helix_end + 1):
        k = col - offset
        bottom = seq2_rev[k] if 0 <= k < len(seq2_rev) else None
        paired = bottom is not None and duplex_bases_pair(
            seq1[col], bottom, consider_ambiguity)
        match_line += "|" if paired else " "
    lines.append(match_line)

    lines.append(f"3'-{' ' * (pad + offset)}{seq2_rev}-5'  ({oligo2.name})")

    return "\n".join(lines)


def visualize_overlap(oligo1, oligo2, overlap_length, consider_ambiguity=False):
    """Create ASCII visualization of the 3' overlap between two oligos.

    oligo1/oligo2 must have .sequence and .name attributes.
    """
    # A 3' overlap is the register where both 3' ends meet, i.e. the reversed
    # oligo2 starts where oligo1's terminal `overlap_length` bases begin.
    offset = len(oligo1.sequence) - overlap_length
    return visualize_duplex(
        oligo1, oligo2, offset, offset, len(oligo1.sequence) - 1,
        consider_ambiguity,
    )


def analyze_mix(oligos, settings):
    """Run pairwise analysis on a list of Oligo objects, dispatching on mode.

    Each oligo must have .id, .name, and .sequence attributes. Both modes return
    a list of dicts sharing the keys primer1_id, primer2_id, primer1_name,
    primer2_name, risk_level and visualization; the remaining keys differ.
    """
    if settings.get("screen_mode", MODE_OVERLAP) == MODE_THERMO:
        return analyze_mix_thermo(oligos, settings)
    return analyze_mix_overlap(oligos, settings)


def analyze_mix_overlap(oligos, settings):
    """Run pairwise 3' overlap analysis on a list of Oligo objects.

    Returns a list of result dicts, each containing:
        overlap_length, mismatches, primer1_id, primer2_id,
        primer1_name, primer2_name, risk_level, visualization
    """
    if len(oligos) < 2:
        return []

    min_overlap = settings.get("min_overlap", 3)
    max_overlap = settings.get("max_overlap", 10)
    max_mismatches = settings.get("max_mismatches", 1)
    consider_ambiguity = settings.get("consider_ambiguity", False)

    results = []

    for overlap_length in range(max_overlap, min_overlap - 1, -1):
        for mm in range(max_mismatches + 1):
            for i in range(len(oligos)):
                for j in range(i, len(oligos)):
                    oligo_i = oligos[i]
                    oligo_j = oligos[j]

                    if len(oligo_i.sequence) < overlap_length or len(oligo_j.sequence) < overlap_length:
                        continue

                    primer1_3end = get_last_n_bases(oligo_i.sequence, overlap_length)
                    primer2_3end_rc = get_last_n_bases_rc(oligo_j.sequence, overlap_length)

                    actual_mm = count_mismatches(primer1_3end, primer2_3end_rc, consider_ambiguity)

                    if actual_mm == mm:
                        risk = get_risk_level(overlap_length, mm, settings)
                        vis = visualize_overlap(oligo_i, oligo_j, overlap_length, consider_ambiguity)

                        results.append({
                            "overlap_length": overlap_length,
                            "mismatches": mm,
                            "primer1_id": oligo_i.id,
                            "primer2_id": oligo_j.id,
                            "primer1_name": oligo_i.name,
                            "primer2_name": oligo_j.name,
                            "risk_level": risk,
                            "visualization": vis,
                        })

    return results


def analyze_mix_thermo(oligos, settings):
    """Score every oligo pair by nearest-neighbor free energy.

    Unlike the overlap mode this emits one row per pair (self-pairs included),
    each carrying dg_min, dg_ens and dg_3prime. Pairs whose interaction is
    weaker than both reporting thresholds are omitted.

    Returns a list of result dicts, each containing:
        dg_min, dg_ens, dg_3prime, dg_3prime_1, dg_3prime_2, n_structures,
        primer1_id, primer2_id, primer1_name, primer2_name, risk_level,
        visualization
    """
    if len(oligos) < 2:
        return []

    model = model_from_settings(settings)
    consider_ambiguity = settings.get("consider_ambiguity", False)
    report_dg_min = settings.get("report_dg_min", -2.0)
    report_dg_ens = settings.get("report_dg_ens", -3.0)

    results = []

    for i in range(len(oligos)):
        for j in range(i, len(oligos)):
            oligo_i = oligos[i]
            oligo_j = oligos[j]

            r = model.pair_interaction(oligo_i.sequence, oligo_j.sequence)
            if r["n_structures"] == 0:
                continue
            # either filter alone is enough to surface the pair
            if r["dg_min"] > report_dg_min and r["dg_ens"] > report_dg_ens:
                continue

            risk = get_thermo_risk_level(r["dg_3prime"], r["dg_min"], settings)
            vis = visualize_duplex(
                oligo_i, oligo_j, r["offset"], r["helix_start"], r["helix_end"],
                consider_ambiguity,
            )

            results.append({
                "dg_min": r["dg_min"],
                "dg_ens": r["dg_ens"],
                "dg_3prime": r["dg_3prime"],
                "dg_3prime_1": r["dg_3prime_1"],
                "dg_3prime_2": r["dg_3prime_2"],
                "n_structures": r["n_structures"],
                "primer1_id": oligo_i.id,
                "primer2_id": oligo_j.id,
                "primer1_name": oligo_i.name,
                "primer2_name": oligo_j.name,
                "risk_level": risk,
                "visualization": vis,
            })

    results.sort(key=lambda r: (r["dg_3prime"], r["dg_min"]))
    return results


def get_interactions_for_oligo(oligo_id, results):
    """Return all results involving a given oligo (by UUID)."""
    return [r for r in results if r["primer1_id"] == oligo_id or r["primer2_id"] == oligo_id]


def get_max_risk_for_oligo(oligo_id, results):
    """Return the highest risk level for an oligo across all its interactions."""
    interactions = get_interactions_for_oligo(oligo_id, results)
    if not interactions:
        return "NONE"
    risk_order = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    max_risk = max(interactions, key=lambda r: risk_order.get(r["risk_level"], 0))
    return max_risk["risk_level"]


def parse_fasta_text(text):
    """Parse FASTA-formatted text and return list of SeqRecords."""
    text = text.replace("\xa0", " ")
    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
        return records
    except Exception:
        return []
