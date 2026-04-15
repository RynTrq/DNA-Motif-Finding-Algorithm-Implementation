import json

import pytest

from dna_motif.cli import main
from dna_motif.io import load_sequences, parse_sequences


def test_parse_plain_text_sequences():
    assert parse_sequences(["ACGT", "", "tgca"]) == ("ACGT", "TGCA")


def test_parse_fasta_sequences_with_wrapped_records():
    assert parse_sequences([">seq1", "AC", "GT", ">seq2", "TGCA"]) == ("ACGT", "TGCA")


def test_load_sequences_reports_missing_file(tmp_path):
    missing = tmp_path / "missing.txt"

    with pytest.raises(FileNotFoundError, match="does not exist"):
        load_sequences(missing)


def test_cli_emits_json_for_median_search(tmp_path, capsys):
    input_file = tmp_path / "dna.txt"
    input_file.write_text("AAACCC\nGGACCT\n", encoding="utf-8")

    exit_code = main([str(input_file), "--algorithm", "median", "-k", "3", "--json"])

    assert exit_code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["algorithm"] == "median_string"
    assert payload["pattern"] == "ACC"
    assert payload["distance"] == 0


def test_cli_returns_usage_error_for_bad_input(tmp_path, capsys):
    input_file = tmp_path / "dna.txt"
    input_file.write_text("ACNT\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc:
        main([str(input_file), "-k", "3"])

    assert exc.value.code == 2
    assert "found: N" in capsys.readouterr().err


def test_cli_omits_gibbs_history_unless_requested(tmp_path, capsys):
    input_file = tmp_path / "dna.txt"
    input_file.write_text("AAACCC\nGGACCT\n", encoding="utf-8")

    exit_code = main(
        [
            str(input_file),
            "--algorithm",
            "gibbs",
            "-k",
            "3",
            "--iterations",
            "3",
            "--restarts",
            "2",
            "--seed",
            "1",
            "--json",
        ]
    )
    without_history = json.loads(capsys.readouterr().out)

    exit_code_with_history = main(
        [
            str(input_file),
            "--algorithm",
            "gibbs",
            "-k",
            "3",
            "--iterations",
            "3",
            "--restarts",
            "2",
            "--seed",
            "1",
            "--json",
            "--include-history",
        ]
    )
    with_history = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert exit_code_with_history == 0
    assert "score_history" not in without_history
    assert len(with_history["score_history"]) == 6
