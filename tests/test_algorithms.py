import pytest

from dna_motif.algorithms import (
    brute_force_search,
    consensus,
    gibbs_sampler,
    hamming_distance,
    median_string_search,
    motif_score,
    normalize_sequences,
    total_hamming_distance,
)


SAMPLE_SEQUENCES = (
    "gcggaagagggcactagcccatgtgagagggcaaggacca",
    "atctttctcttaaaaataacataattcagggccaggatgt",
    "gtcacgagctttatcctacagatgatgaatgcaaatcagc",
    "taaaagataatatcgaccctagcgtggcgggcaaggtgct",
)


def test_normalize_sequences_rejects_empty_and_ambiguous_input():
    assert normalize_sequences([" acgt ", "", "TGCA"]) == ("ACGT", "TGCA")

    with pytest.raises(ValueError, match="at least one"):
        normalize_sequences(["", "   "])

    with pytest.raises(ValueError, match="found: N"):
        normalize_sequences(["ACNT"])


def test_hamming_distance_validates_lengths():
    assert hamming_distance("GATTACA", "GACTATA") == 2

    with pytest.raises(ValueError, match="equal length"):
        hamming_distance("AAA", "AA")


def test_consensus_and_score_are_deterministic_with_ties():
    motifs = ("ACG", "ACC", "ATG", "AAG")

    assert consensus(motifs) == "ACG"
    assert motif_score(motifs) == 9

    with pytest.raises(ValueError, match="same length"):
        consensus(("AAA", "AA"))


def test_total_hamming_distance_returns_closest_positions_and_motifs():
    distance, positions, motifs = total_hamming_distance(
        ["AAACCC", "GGACCT"],
        "ACC",
    )

    assert distance == 0
    assert positions == (2, 2)
    assert motifs == ("ACC", "ACC")


def test_brute_force_search_finds_known_sample_alignment():
    result = brute_force_search(SAMPLE_SEQUENCES, 6)

    assert result.algorithm == "brute_force"
    assert result.motifs == ("GGGCAA", "GGGCCA", "ATGCAA", "GGGCAA")
    assert result.positions == (28, 28, 28, 28)
    assert result.consensus == "GGGCAA"
    assert result.score == 21
    assert result.states_evaluated == 1_500_625


def test_brute_force_guard_prevents_accidental_explosion():
    with pytest.raises(ValueError, match="would evaluate"):
        brute_force_search(SAMPLE_SEQUENCES, 6, max_states=100)


def test_median_string_search_finds_known_sample_pattern():
    result = median_string_search(SAMPLE_SEQUENCES, 6)

    assert result.algorithm == "median_string"
    assert result.pattern == "GCAAGG"
    assert result.distance == 3
    assert result.positions == (30, 30, 30, 30)
    assert result.motifs == ("GCAAGG", "GCCAGG", "GCAAAT", "GCAAGG")
    assert result.score == 21
    assert result.states_evaluated == 4_096


def test_gibbs_sampler_is_reproducible_with_seed_and_tracks_history():
    first = gibbs_sampler(SAMPLE_SEQUENCES, 6, iterations=50, restarts=5, seed=7)
    second = gibbs_sampler(SAMPLE_SEQUENCES, 6, iterations=50, restarts=5, seed=7)

    assert first == second
    assert first.algorithm == "gibbs_sampler"
    assert len(first.score_history) == 250
    assert first.score == max(first.score_history)
    assert len(first.motifs) == len(SAMPLE_SEQUENCES)
    assert all(len(motif) == 6 for motif in first.motifs)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"iterations": 0}, "iterations"),
        ({"restarts": 0}, "restarts"),
        ({"pseudocount": 0}, "pseudocount"),
    ],
)
def test_gibbs_sampler_validates_configuration(kwargs, message):
    with pytest.raises(ValueError, match=message):
        gibbs_sampler(SAMPLE_SEQUENCES, 6, **kwargs)
