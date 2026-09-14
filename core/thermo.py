"""Nearest-neighbor free-energy scoring for oligo dimer interactions.

Uses the SantaLucia/Allawi parameter set shipped with Biopython
(``DNA_NN3`` stacking, ``DNA_IMM1`` internal mismatches, ``DNA_TMM1`` terminal
mismatches, ``DNA_DE1`` dangling ends), so no extra dependency is needed.

Model and its limits
--------------------
For a pair of oligos every gapless register (relative offset) is scanned. Within
a register the single most stable contiguous helix is found; a helix must begin
and end on a Watson-Crick pair but may carry mismatches internally. Helix
termini get the standard initiation penalty plus a dangling-end or terminal-
mismatch contribution from the flanking column. Energies are salt-corrected per
SantaLucia (1998).

Deliberately *not* modelled: bulges, internal loops and multi-helix structures
within one register. Values are therefore slightly less negative than primer3's
for structures that rely on a bulge, and will not match it exactly.

Three numbers are reported per pair:

``dg_min``
    The most stable single structure found - the classic worst-case dimer.
``dg_ens``
    Ensemble free energy over all registers, ``-RT ln sum(exp(-dG_i/RT))``.
    Always <= dg_min; the gap indicates how many competing structures exist.
    This is the meaningful "total" - an arithmetic sum would instead scale with
    how many registers happen to be scanned.
``dg_3prime``
    Stability of the most stable helix that actually pairs an oligo's 3'
    terminal base. Structures leaving both 3' ends unpaired cannot be extended
    by polymerase, so this is the number that best predicts primer-dimer.
"""

import math
from functools import lru_cache

from Bio.SeqUtils import MeltingTemp as mt

from core.sequtils import duplex_bases_pair, expand_base

R = 1.9872e-3  # gas constant, kcal/(mol*K)

MIN_HELIX_BP = 2


def _dg(entry, temp_k):
    """Convert a (dH kcal/mol, dS cal/mol*K) table entry to dG at temp_k."""
    dh, ds = entry
    return dh - temp_k * ds / 1000.0


class ThermoModel:
    """Temperature- and salt-specific energy tables plus duplex scoring."""

    def __init__(self, temp_c=37.0, na_mM=50.0, consider_ambiguity=False):
        self.temp_c = temp_c
        self.temp_k = temp_c + 273.15
        self.consider_ambiguity = consider_ambiguity
        self.rt = R * self.temp_k

        # SantaLucia (1998) salt correction, expressed per phosphate on dG
        na_M = max(na_mM, 1e-3) / 1000.0
        self.salt_per_phosphate = -self.temp_k * 0.368 * math.log(na_M) / 1000.0

        self._stack = {}
        for key, entry in mt.DNA_NN3.items():
            if key.startswith("init") or key == "sym":
                continue
            self._stack[key] = _dg(entry, self.temp_k)
        for key, entry in mt.DNA_IMM1.items():
            self._stack.setdefault(key, _dg(entry, self.temp_k))

        self._tmm = {k: _dg(v, self.temp_k) for k, v in mt.DNA_TMM1.items()}
        self._de = {k: _dg(v, self.temp_k) for k, v in mt.DNA_DE1.items()}

        self._init_at = _dg(mt.DNA_NN3["init_A/T"], self.temp_k)
        self._init_gc = _dg(mt.DNA_NN3["init_G/C"], self.temp_k)
        self.dg_sym = _dg(mt.DNA_NN3["sym"], self.temp_k)

        # caches are per-model so they never outlive a settings change
        self._stack_dg = lru_cache(maxsize=None)(self._stack_dg_uncached)
        self._end_dg = lru_cache(maxsize=None)(self._end_dg_uncached)
        self._init_dg = lru_cache(maxsize=None)(self._init_dg_uncached)

    # ---- table lookup ----

    def _raw(self, table, t1, t2, b1, b2):
        """Look up a column-pair key, trying the equivalent flipped orientation.

        Reversing the whole "XY/AB" key string views the same duplex from the
        other strand, which is how the one-orientation tables are extended.
        """
        key = f"{t1}{t2}/{b1}{b2}"
        if key in table:
            return table[key]
        flipped = key[::-1]
        if flipped in table:
            return table[flipped]
        return None

    def _resolve(self, table, t1, t2, b1, b2, default=0.0):
        """Table lookup that expands IUPAC codes to their most stable option."""
        value = self._raw(table, t1, t2, b1, b2)
        if value is not None:
            return value
        if not self.consider_ambiguity:
            return default
        best = None
        for a1 in self._options(t1):
            for a2 in self._options(t2):
                for c1 in self._options(b1):
                    for c2 in self._options(b2):
                        value = self._raw(table, a1, a2, c1, c2)
                        if value is not None and (best is None or value < best):
                            best = value
        return default if best is None else best

    @staticmethod
    def _options(base):
        """Concrete bases to try for a column, passing the gap marker through."""
        return "." if base == "." else expand_base(base)

    def _stack_dg_uncached(self, t1, t2, b1, b2):
        """Stacking energy between two adjacent duplex columns."""
        return self._resolve(self._stack, t1, t2, b1, b2)

    def _init_dg_uncached(self, top, bottom):
        """Helix initiation penalty for a terminal base pair."""
        if self.consider_ambiguity and expand_base(top) - {"A", "T"}:
            return self._init_gc
        return self._init_at if top in ("A", "T") else self._init_gc

    def _end_dg_uncached(self, t_in, t_out, b_in, b_out):
        """Contribution of the column flanking a helix terminus.

        ``*_in`` is the terminal base pair of the helix, ``*_out`` the flanking
        column: both bases present is a terminal mismatch, one present a
        dangling end, neither nothing.
        """
        if t_out == "." and b_out == ".":
            return 0.0
        if t_out != "." and b_out != ".":
            return self._resolve(self._tmm, t_in, t_out, b_in, b_out)
        return self._resolve(self._de, t_in, t_out, b_in, b_out)

    # ---- duplex scoring ----

    def _register(self, top, bot_rev, off, self_pair):
        """Score one gapless register of two strands.

        ``top`` is oligo A 5'->3'. ``bot_rev`` is oligo B reversed, so it reads
        3'->5' left to right and column c of A faces ``bot_rev[c - off]``.
        Returns (dg, helix_start, helix_end, dg_3p_a, dg_3p_b) in A coordinates,
        or None when the register holds no stabilizing helix.
        """
        la, lb = len(top), len(bot_rev)
        start = max(0, off)
        stop = min(la, off + lb)  # exclusive
        if stop - start < MIN_HELIX_BP:
            return None

        def base_top(c):
            return top[c] if 0 <= c < la else "."

        def base_bot(c):
            k = c - off
            return bot_rev[k] if 0 <= k < lb else "."

        cols = range(start, stop)
        anchors = [c for c in cols
                   if duplex_bases_pair(top[c], base_bot(c), self.consider_ambiguity)]
        if len(anchors) < MIN_HELIX_BP:
            return None

        # running stack sum from `start`, so sum over [i, j] == prefix[j] - prefix[i]
        prefix = {start: 0.0}
        total = 0.0
        for c in range(start, stop - 1):
            total += self._stack_dg(top[c], top[c + 1], base_bot(c), base_bot(c + 1))
            prefix[c + 1] = total

        salt = self.salt_per_phosphate
        # terminal costs are split so each depends on a single endpoint; the
        # salt term is linear in helix length and splits the same way
        cost_i, cost_j = {}, {}
        for c in anchors:
            t, b = top[c], base_bot(c)
            init = self._init_dg(t, b)
            cost_i[c] = (-prefix[c] + init - salt * c
                         + self._end_dg(t, base_top(c - 1), b, base_bot(c - 1)))
            cost_j[c] = (prefix[c] + init + salt * c
                         + self._end_dg(t, base_top(c + 1), b, base_bot(c + 1)))

        def with_symmetry(value, i, j):
            """Self-complementary structures carry the NN symmetry penalty."""
            if self_pair and i + j == la - 1 + off:
                return value + self.dg_sym
            return value

        # best helix overall: sweep j, keeping the best opening i < j
        best = best_i = best_j = None
        run_i = run_i_at = None
        a3 = None
        for j in anchors:
            if run_i is not None:
                value = with_symmetry(run_i + cost_j[j], run_i_at, j)
                if best is None or value < best:
                    best, best_i, best_j = value, run_i_at, j
                if j == la - 1:
                    a3 = value
            if run_i is None or cost_i[j] < run_i:
                run_i, run_i_at = cost_i[j], j

        if best is None or best >= 0:
            return None

        # constrained to helices starting at B's 3' base, which sits at column off
        b3 = None
        if off in cost_i:
            tail = [(cost_i[off] + cost_j[j], j) for j in anchors if j > off]
            if tail:
                value, j = min(tail)
                b3 = with_symmetry(value, off, j)

        return (best, best_i, best_j,
                a3 if a3 is not None and a3 < 0 else 0.0,
                b3 if b3 is not None and b3 < 0 else 0.0)

    def pair_interaction(self, seq_a, seq_b):
        """Score every gapless register of two oligos.

        Returns a dict with dg_min, dg_ens, dg_3prime_1 / dg_3prime_2 (oligo A's
        and oligo B's 3' end respectively), dg_3prime (the worse of the two) and
        the offset/span of the best helix for visualization.
        """
        a = seq_a.upper()
        b_rev = seq_b.upper()[::-1]
        la, lb = len(a), len(b_rev)
        self_pair = a == seq_b.upper()

        empty = {
            "dg_min": 0.0, "dg_ens": 0.0, "dg_3prime": 0.0,
            "dg_3prime_1": 0.0, "dg_3prime_2": 0.0,
            "offset": 0, "helix_start": 0, "helix_end": 0, "n_structures": 0,
        }
        if la < MIN_HELIX_BP or lb < MIN_HELIX_BP:
            return empty

        energies = []
        best = None
        dg3_a = 0.0
        dg3_b = 0.0
        for off in range(-lb + 1, la):
            scored = self._register(a, b_rev, off, self_pair)
            if scored is None:
                continue
            dg, helix_start, helix_end, a3, b3 = scored
            energies.append(dg)
            if best is None or dg < best[0]:
                best = (dg, off, helix_start, helix_end)
            dg3_a = min(dg3_a, a3)
            dg3_b = min(dg3_b, b3)

        if not energies:
            return empty

        z = sum(math.exp(-dg / self.rt) for dg in energies)
        dg_ens = -self.rt * math.log(z) if z > 0 else 0.0

        return {
            "dg_min": best[0],
            "dg_ens": dg_ens,
            "dg_3prime": min(dg3_a, dg3_b),
            "dg_3prime_1": dg3_a,
            "dg_3prime_2": dg3_b,
            "offset": best[1],
            "helix_start": best[2],
            "helix_end": best[3],
            "n_structures": len(energies),
        }


def model_from_settings(settings):
    """Build a ThermoModel from the app settings dict."""
    return ThermoModel(
        temp_c=settings.get("dg_temperature", 37.0),
        na_mM=settings.get("na_concentration", 50.0),
        consider_ambiguity=settings.get("consider_ambiguity", False),
    )
