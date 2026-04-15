"""Input and output helpers for DNA motif finding."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from dna_motif.algorithms import normalize_sequences


def load_sequences(path: str | Path) -> tuple[str, ...]:
    """Load plain-text or FASTA DNA sequences from disk."""

    input_path = Path(path)
    if not input_path.exists():
        raise FileNotFoundError(f"DNA input file does not exist: {input_path}")
    if not input_path.is_file():
        raise ValueError(f"DNA input path is not a file: {input_path}")

    lines = input_path.read_text(encoding="utf-8").splitlines()
    return parse_sequences(lines)


def parse_sequences(lines: Iterable[str]) -> tuple[str, ...]:
    """Parse either one-sequence-per-line text or FASTA content."""

    raw_lines = [line.strip() for line in lines if line.strip()]
    if not raw_lines:
        raise ValueError("input does not contain any DNA sequences")

    if any(line.startswith(">") for line in raw_lines):
        sequences: list[str] = []
        current: list[str] = []
        for line in raw_lines:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
                continue
            current.append(line)
        if current:
            sequences.append("".join(current))
        return normalize_sequences(sequences)

    return normalize_sequences(raw_lines)
