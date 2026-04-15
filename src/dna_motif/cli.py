"""Command-line interface for DNA motif finding."""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict

from dna_motif.algorithms import brute_force_search, gibbs_sampler, median_string_search
from dna_motif.io import load_sequences


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dna-motif",
        description="Find DNA motifs with brute force, median string, or Gibbs sampling.",
    )
    parser.add_argument("input", help="Path to a plain-text or FASTA DNA sequence file.")
    parser.add_argument("-k", "--motif-length", type=int, default=6, help="Motif length.")
    parser.add_argument(
        "-a",
        "--algorithm",
        choices=("brute-force", "median", "gibbs"),
        default="gibbs",
        help="Algorithm to run.",
    )
    parser.add_argument("--iterations", type=int, default=1_000, help="Gibbs iterations per restart.")
    parser.add_argument("--restarts", type=int, default=20, help="Gibbs restart count.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible Gibbs runs.")
    parser.add_argument(
        "--include-history",
        action="store_true",
        help="Include the full Gibbs score history in JSON output.",
    )
    parser.add_argument(
        "--max-states",
        type=int,
        default=10_000_000,
        help="Safety cap for brute-force start-position combinations.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    sequences = load_sequences(args.input)
    start = time.perf_counter()

    if args.algorithm == "brute-force":
        result = brute_force_search(sequences, args.motif_length, max_states=args.max_states)
    elif args.algorithm == "median":
        result = median_string_search(sequences, args.motif_length)
    else:
        result = gibbs_sampler(
            sequences,
            args.motif_length,
            iterations=args.iterations,
            restarts=args.restarts,
            seed=args.seed,
        )

    elapsed_seconds = time.perf_counter() - start
    payload = asdict(result)
    if not args.include_history:
        payload.pop("score_history", None)
    payload["elapsed_seconds"] = elapsed_seconds
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        payload = run(args)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    print(f"Algorithm: {payload['algorithm']}")
    print(f"Motif length: {payload['k']}")
    if payload.get("pattern"):
        print(f"Median pattern: {payload['pattern']}")
    print(f"Consensus: {payload['consensus']}")
    print(f"Motifs: {', '.join(payload['motifs'])}")
    print(f"Positions: {', '.join(str(position) for position in payload['positions'])}")
    print(f"Score: {payload['score']}")
    if payload.get("distance") is not None:
        print(f"Distance: {payload['distance']}")
    print(f"States evaluated: {payload['states_evaluated']}")
    print(f"Elapsed seconds: {payload['elapsed_seconds']:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
