"""Shared sequence primitives: IUPAC ambiguity codes and base comparison.

Kept separate from analysis/thermo so both can depend on it without a cycle.
"""

from functools import lru_cache

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

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


@lru_cache(maxsize=None)
def expand_base(base):
    """Return the set of concrete bases an IUPAC code can represent."""
    return frozenset(IUPAC_CODES.get(base.upper(), {base.upper()}))


@lru_cache(maxsize=None)
def bases_could_match(base1, base2):
    """Check if two bases could be identical considering IUPAC ambiguity codes."""
    return bool(expand_base(base1) & expand_base(base2))


@lru_cache(maxsize=None)
def bases_could_pair(top, bottom):
    """Check if two bases sitting opposite each other could form a Watson-Crick pair.

    Unlike bases_could_match this compares a base against its *complement*, for
    use where the two strands are held in duplex orientation rather than one
    being reverse-complemented first.
    """
    comps = {COMPLEMENT[b] for b in expand_base(top) if b in COMPLEMENT}
    return bool(comps & expand_base(bottom))


# Complement of each IUPAC code, matching Bio.Seq.reverse_complement
COMPLEMENT_CODE = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


@lru_cache(maxsize=None)
def duplex_bases_pair(top, bottom, consider_ambiguity=False):
    """Check whether two bases held opposite each other count as paired.

    With ambiguity off the codes must be exact complements, which is what
    comparing a sequence against a reverse complement does. With it on, any
    shared expansion counts.
    """
    if consider_ambiguity:
        return bases_could_pair(top, bottom)
    return COMPLEMENT_CODE.get(bottom.upper()) == top.upper()
