import json
import uuid
from dataclasses import dataclass, field

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


@dataclass
class Oligo:
    name: str
    sequence: str
    id: str = field(default_factory=lambda: str(uuid.uuid4()))
    active: bool = True
    mix_id: str | None = None

    def to_seq_record(self) -> SeqRecord:
        return SeqRecord(Seq(self.sequence), id=self.name, description="")

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "name": self.name,
            "sequence": self.sequence,
            "active": self.active,
            "mix_id": self.mix_id,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Oligo":
        return cls(**data)

    @classmethod
    def from_seq_record(cls, record: SeqRecord, mix_id: str | None = None) -> "Oligo":
        return cls(name=record.id, sequence=str(record.seq), mix_id=mix_id)

    def gc_content(self) -> float:
        seq = self.sequence.upper()
        if not seq:
            return 0.0
        gc = sum(1 for b in seq if b in ("G", "C"))
        return gc / len(seq) * 100


@dataclass
class Mix:
    name: str
    id: str = field(default_factory=lambda: str(uuid.uuid4()))

    def to_dict(self) -> dict:
        return {"id": self.id, "name": self.name}

    @classmethod
    def from_dict(cls, data: dict) -> "Mix":
        return cls(**data)


class Project:
    def __init__(self):
        self.oligos: list[Oligo] = []
        self.mixes: list[Mix] = []

    def get_oligos_in_mix(self, mix_id: str) -> list[Oligo]:
        return [o for o in self.oligos if o.mix_id == mix_id]

    def get_active_oligos_in_mix(self, mix_id: str) -> list[Oligo]:
        return [o for o in self.oligos if o.mix_id == mix_id and o.active]

    def get_unallocated_oligos(self) -> list[Oligo]:
        return [o for o in self.oligos if o.mix_id is None]

    def get_oligo_by_id(self, oligo_id: str) -> Oligo | None:
        for o in self.oligos:
            if o.id == oligo_id:
                return o
        return None

    def get_mix_by_id(self, mix_id: str) -> Mix | None:
        for m in self.mixes:
            if m.id == mix_id:
                return m
        return None

    def add_oligo(self, oligo: Oligo):
        self.oligos.append(oligo)

    def remove_oligo(self, oligo_id: str):
        self.oligos = [o for o in self.oligos if o.id != oligo_id]

    def add_mix(self, mix: Mix):
        self.mixes.append(mix)

    def remove_mix(self, mix_id: str):
        for o in self.oligos:
            if o.mix_id == mix_id:
                o.mix_id = None
        self.mixes = [m for m in self.mixes if m.id != mix_id]

    def duplicate_oligo(self, oligo_id: str) -> Oligo | None:
        original = self.get_oligo_by_id(oligo_id)
        if not original:
            return None
        base_name = original.name
        existing_names = {o.name for o in self.oligos}
        count = 1
        while f"{base_name}_{count}" in existing_names:
            count += 1
        dup = Oligo(
            name=f"{base_name}_{count}",
            sequence=original.sequence,
            active=original.active,
            mix_id=original.mix_id,
        )
        self.oligos.append(dup)
        return dup

    def to_dict(self) -> dict:
        return {
            "mixes": [m.to_dict() for m in self.mixes],
            "oligos": [o.to_dict() for o in self.oligos],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Project":
        proj = cls()
        proj.mixes = [Mix.from_dict(m) for m in data.get("mixes", [])]
        proj.oligos = [Oligo.from_dict(o) for o in data.get("oligos", [])]
        return proj

    def save(self, path: str):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load(cls, path: str) -> "Project":
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)

    def export_fasta(self, path: str, mix_id: str | None = "ALL"):
        """Export oligos as FASTA. mix_id=None for unallocated, 'ALL' for everything, or a specific mix id."""
        if mix_id == "ALL":
            oligos = self.oligos
        elif mix_id is None:
            oligos = self.get_unallocated_oligos()
        else:
            oligos = self.get_oligos_in_mix(mix_id)
        records = [o.to_seq_record() for o in oligos]
        with open(path, "w") as f:
            SeqIO.write(records, f, "fasta")
        return len(records)
