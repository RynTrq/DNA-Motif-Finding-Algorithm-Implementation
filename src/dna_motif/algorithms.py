"""Production-ready implementations of classic DNA motif finding algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from random import Random
from typing import Iterable, Sequence

DNA_ALPHABET = "ACGT"


@dataclass(frozen=True)
class MotifResult:
    """Shared result shape for deterministic motif search algorithms."""

    algorithm: str
    k: int
    motifs: tuple[str, ...]
    positions: tuple[int, ...]
    consensus: str
    score: int
    distance: int | None = None
    pattern: str | None = None
    states_evaluated: int | None = None


@dataclass(frozen=True)
class GibbsResult(MotifResult):
    """Result shape for Gibbs sampling with convergence history."""

    score_history: tuple[int, ...] = ()


def normalize_sequences(sequences: Iterable[str]) -> tuple[str, ...]:
    """Normalize and validate DNA sequences.

    Empty lines are ignored so input files can contain spacing. The algorithms
    intentionally reject ambiguous bases because the scoring/profile logic is
    defined over the canonical DNA alphabet only.
    """

    normalized = tuple(sequence.strip().upper() for sequence in sequences if sequence.strip())
    if not normalized:
        raise ValueError("at least one DNA sequence is required")

    invalid = {
        base
        for sequence in normalized
        for base in sequence
        if base not in DNA_ALPHABET
    }
    if invalid:
        bad = ", ".join(sorted(invalid))
        raise ValueError(f"DNA sequences may only contain A, C, G, and T; found: {bad}")

    return normalized


def validate_k(sequences: Sequence[str], k: int) -> None:
    if k <= 0:
        raise ValueError("motif length k must be positive")
    shortest = min(len(sequence) for sequence in sequences)
    if k > shortest:
        raise ValueError(
            f"motif length k={k} exceeds the shortest sequence length ({shortest})"
        )


def hamming_distance(left: str, right: str) -> int:
    """Return the Hamming distance between equal-length strings."""

    if len(left) != len(right):
        raise ValueError("hamming distance requires strings of equal length")
    return sum(a != b for a, b in zip(left.upper(), right.upper()))


def consensus(motifs: Sequence[str]) -> str:
    """Return the deterministic consensus sequence for a set of motifs.

    Ties are resolved alphabetically by the canonical alphabet order A, C, G, T.
    This keeps outputs stable across runs and Python versions.
    """

    if not motifs:
        raise ValueError("at least one motif is required")

    motif_length = len(motifs[0])
    if any(len(motif) != motif_length for motif in motifs):
        raise ValueError("all motifs must have the same length")

    result: list[str] = []
    for column in zip(*(motif.upper() for motif in motifs)):
        counts = {base: column.count(base) for base in DNA_ALPHABET}
        result.append(max(DNA_ALPHABET, key=lambda base: (counts[base], -DNA_ALPHABET.index(base))))
    return "".join(result)


def motif_score(motifs: Sequence[str]) -> int:
    """Return total consensus support across all motif columns."""

    if not motifs:
        return 0
    motif_length = len(motifs[0])
    if any(len(motif) != motif_length for motif in motifs):
        raise ValueError("all motifs must have the same length")
    return sum(
        max(column.count(base) for base in DNA_ALPHABET)
        for column in zip(*(motif.upper() for motif in motifs))
    )


def motifs_from_positions(sequences: Sequence[str], positions: Sequence[int], k: int) -> tuple[str, ...]:
    """Extract k-mers from each sequence at the supplied start positions."""

    if len(sequences) != len(positions):
        raise ValueError("positions must contain exactly one start index per sequence")
    return tuple(sequence[position : position + k] for sequence, position in zip(sequences, positions))


def score_positions(sequences: Sequence[str], positions: Sequence[int], k: int) -> int:
    """Score one candidate motif alignment."""

    return motif_score(motifs_from_positions(sequences, positions, k))


def brute_force_search(
    sequences: Iterable[str],
    k: int,
    *,
    max_states: int | None = 10_000_000,
) -> MotifResult:
    """Exhaustively search all motif start-position combinations.

    This is exact for the "best aligned motif instances" formulation, but its
    state space is the product of possible starts in every sequence. The
    ``max_states`` guard prevents accidental multi-hour runs from large inputs.
    """

    dna = normalize_sequences(sequences)
    validate_k(dna, k)
    ranges = [range(len(sequence) - k + 1) for sequence in dna]
    states = prod(len(candidate_range) for candidate_range in ranges)
    if max_states is not None and states > max_states:
        raise ValueError(
            f"brute force would evaluate {states:,} states; raise max_states to run it"
        )

    best_positions: tuple[int, ...] | None = None
    best_score = -1
    evaluated = 0
    for positions in product(*ranges):
        evaluated += 1
        score = score_positions(dna, positions, k)
        if score > best_score:
            best_score = score
            best_positions = tuple(positions)

    assert best_positions is not None
    best_motifs = motifs_from_positions(dna, best_positions, k)
    return MotifResult(
        algorithm="brute_force",
        k=k,
        motifs=best_motifs,
        positions=best_positions,
        consensus=consensus(best_motifs),
        score=best_score,
        states_evaluated=evaluated,
    )


def total_hamming_distance(
    sequences: Iterable[str],
    pattern: str,
) -> tuple[int, tuple[int, ...], tuple[str, ...]]:
    """Return distance from a candidate pattern to its closest k-mer in each sequence."""

    dna = normalize_sequences(sequences)
    candidate = pattern.strip().upper()
    if not candidate:
        raise ValueError("pattern must not be empty")
    validate_k(dna, len(candidate))

    total = 0
    positions: list[int] = []
    motifs: list[str] = []
    for sequence in dna:
        best_distance = len(candidate) + 1
        best_position = 0
        best_motif = sequence[: len(candidate)]
        for position in range(len(sequence) - len(candidate) + 1):
            motif = sequence[position : position + len(candidate)]
            distance = hamming_distance(motif, candidate)
            if distance < best_distance:
                best_distance = distance
                best_position = position
                best_motif = motif
        total += best_distance
        positions.append(best_position)
        motifs.append(best_motif)
    return total, tuple(positions), tuple(motifs)


def median_string_search(sequences: Iterable[str], k: int) -> MotifResult:
    """Find the k-mer minimizing total Hamming distance to all DNA sequences."""

    dna = normalize_sequences(sequences)
    validate_k(dna, k)

    best_pattern: str | None = None
    best_distance = k * len(dna) + 1
    best_positions: tuple[int, ...] = ()
    best_motifs: tuple[str, ...] = ()
    evaluated = 0

    for letters in product(DNA_ALPHABET, repeat=k):
        evaluated += 1
        pattern = "".join(letters)
        distance, positions, motifs = total_hamming_distance(dna, pattern)
        if distance < best_distance:
            best_pattern = pattern
            best_distance = distance
            best_positions = positions
            best_motifs = motifs

    assert best_pattern is not None
    return MotifResult(
        algorithm="median_string",
        k=k,
        motifs=best_motifs,
        positions=best_positions,
        consensus=consensus(best_motifs),
        score=motif_score(best_motifs),
        distance=best_distance,
        pattern=best_pattern,
        states_evaluated=evaluated,
    )


def _profile_from_motifs(motifs: Sequence[str], k: int, pseudocount: float) -> list[dict[str, float]]:
    denominator = len(motifs) + pseudocount * len(DNA_ALPHABET)
    return [
        {
            base: (column.count(base) + pseudocount) / denominator
            for base in DNA_ALPHABET
        }
        for column in zip(*motifs)
    ] if motifs else [{base: 1.0 / len(DNA_ALPHABET) for base in DNA_ALPHABET} for _ in range(k)]


def _kmer_probability(kmer: str, profile: Sequence[dict[str, float]]) -> float:
    return prod(profile[index][base] for index, base in enumerate(kmer))


def _weighted_choice_index(weights: Sequence[float], rng: Random) -> int:
    total = sum(weights)
    if total <= 0:
        return rng.randrange(len(weights))

    threshold = rng.random() * total
    cumulative = 0.0
    for index, weight in enumerate(weights):
        cumulative += weight
        if cumulative >= threshold:
            return index
    return len(weights) - 1


def gibbs_sampler(
    sequences: Iterable[str],
    k: int,
    *,
    iterations: int = 1_000,
    restarts: int = 20,
    seed: int | None = None,
    pseudocount: float = 1.0,
) -> GibbsResult:
    """Approximate motif search using profile-based Gibbs sampling."""

    dna = normalize_sequences(sequences)
    validate_k(dna, k)
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    if restarts <= 0:
        raise ValueError("restarts must be positive")
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")

    rng = Random(seed)
    starts_by_sequence = [range(len(sequence) - k + 1) for sequence in dna]
    best_positions: tuple[int, ...] | None = None
    best_motifs: tuple[str, ...] = ()
    best_score = -1
    score_history: list[int] = []

    for _ in range(restarts):
        positions = [rng.choice(list(starts)) for starts in starts_by_sequence]
        motifs = list(motifs_from_positions(dna, positions, k))

        for _ in range(iterations):
            held_out = rng.randrange(len(dna))
            profile = _profile_from_motifs(
                [motif for index, motif in enumerate(motifs) if index != held_out],
                k,
                pseudocount,
            )
            kmers = [
                dna[held_out][position : position + k]
                for position in starts_by_sequence[held_out]
            ]
            weights = [_kmer_probability(kmer, profile) for kmer in kmers]
            sampled_index = _weighted_choice_index(weights, rng)
            positions[held_out] = sampled_index
            motifs[held_out] = kmers[sampled_index]

            current_score = motif_score(motifs)
            score_history.append(current_score)
            if current_score > best_score:
                best_score = current_score
                best_positions = tuple(positions)
                best_motifs = tuple(motifs)

    assert best_positions is not None
    return GibbsResult(
        algorithm="gibbs_sampler",
        k=k,
        motifs=best_motifs,
        positions=best_positions,
        consensus=consensus(best_motifs),
        score=best_score,
        states_evaluated=iterations * restarts,
        score_history=tuple(score_history),
    )
