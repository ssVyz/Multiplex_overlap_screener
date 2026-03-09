from Bio.Seq import Seq
from Bio import SeqIO
from io import StringIO

# IUPAC ambiguity codes — maps each code to its possible bases
IUPAC_CODES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def bases_could_match(base1, base2):
    """Check if two bases could match considering IUPAC ambiguity codes."""
    bases1 = IUPAC_CODES.get(base1.upper(), {base1.upper()})
    bases2 = IUPAC_CODES.get(base2.upper(), {base2.upper()})
    return bool(bases1 & bases2)


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


def get_last_n_bases(sequence, n):
    """Get last n bases from 3' end of a sequence string."""
    return sequence[-n:].upper()


def get_last_n_bases_rc(sequence, n):
    """Get last n bases from 3' end and return reverse complement."""
    trimmed = Seq(sequence[-n:].upper())
    return str(trimmed.reverse_complement())


def visualize_overlap(oligo1, oligo2, overlap_length, consider_ambiguity=False):
    """Create ASCII visualization of the 3' overlap between two oligos.

    oligo1/oligo2 must have .sequence and .name attributes.
    """
    primer1_full = oligo1.sequence
    primer2_full = oligo2.sequence

    p2_offset = len(primer1_full) - overlap_length

    primer1_3end = get_last_n_bases(oligo1.sequence, overlap_length)
    primer2_3end_rc = get_last_n_bases_rc(oligo2.sequence, overlap_length)

    lines = []
    lines.append(f"5'-{primer1_full}-3'  ({oligo1.name})")

    match_line = "   " + " " * p2_offset
    for a, b in zip(primer1_3end, primer2_3end_rc):
        if consider_ambiguity:
            match_line += "|" if bases_could_match(a, b) else " "
        else:
            match_line += "|" if a == b else " "
    lines.append(match_line)

    primer2_reversed = primer2_full[::-1]
    lines.append(f"3'-{' ' * p2_offset}{primer2_reversed}-5'  ({oligo2.name})")

    return "\n".join(lines)


def analyze_mix(oligos, settings):
    """Run pairwise 3' overlap analysis on a list of Oligo objects.

    Each oligo must have .id, .name, and .sequence attributes.

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
