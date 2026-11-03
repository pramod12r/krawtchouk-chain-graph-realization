#!/usr/bin/env python3
"""Independent verifier for a bounded full-path Bellman certificate.

The checker does NOT import the solver. It enumerates c | q_i (rather than
b | q_i), derives m = n*q_i/c^2, and checks the certificate independently.
"""
from __future__ import annotations
import argparse
import json
from collections import defaultdict
from math import isqrt
from pathlib import Path


def divisors_independent(value: int):
    """Independent divisor enumeration using paired factors."""
    left, right = [], []
    k = 1
    while k * k <= value:
        quotient, remainder = divmod(value, k)
        if remainder == 0:
            left.append(k)
            if quotient != k:
                right.append(quotient)
        k += 1
    return left + right[::-1]


def checker_successors(i: int, D: int, n: int):
    """Checker-side transition enumeration from c divisors."""
    q = (i + 1) * (D - i)
    for c in divisors_independent(q):
        if c > n:
            continue
        # From balance and bc=q: m = n*b/c = n*q/c^2.
        numerator = n * q
        denominator = c * c
        m, rem = divmod(numerator, denominator)
        if rem:
            continue
        b, rem_b = divmod(q, c)
        if rem_b:
            continue
        if b > m:
            continue
        yield b, c, m


def check_witness(D, witness):
    n = witness["layer_sizes"]
    b = witness["forward_degrees"]
    c = witness["backward_degrees"]
    assert len(n) == D + 1 and len(b) == D and len(c) == D
    assert n[0] == 1 and n[-1] == 1
    for i in range(D):
        q = (i + 1) * (D - i)
        assert all(isinstance(x, int) and x > 0 for x in (n[i], n[i+1], b[i], c[i]))
        assert b[i] * c[i] == q
        assert n[i] * b[i] == n[i+1] * c[i]
        assert b[i] <= n[i+1] and c[i] <= n[i]


def verify(path: Path):
    cert = json.loads(path.read_text(encoding="utf-8"))
    D, U = cert["dimension"], cert["bound_U"]
    states = cert["states"]
    transitions = cert["transitions"]

    # Build layer maps and reject duplicate/malformed state entries.
    labels = [dict() for _ in range(D + 1)]
    state_records = {}
    for s in states:
        i, n, cost = s["layer"], s["n"], s["cost"]
        assert 0 <= i <= D and n > 0 and cost > 0
        assert n not in labels[i], "Duplicate state."
        labels[i][n] = cost
        state_records[(i, n)] = s
    assert labels[0] == {1: 1}, "Initial state must be exactly (0,1) with cost 1."

    # Independently reconstruct every permissible transition from every
    # listed source state, bounded by U.
    expected_transitions = set()
    expected_best = [dict() for _ in range(D + 1)]
    expected_best[0][1] = 1
    for i in range(D):
        for n, cost in labels[i].items():
            for b, c, m in checker_successors(i, D, n):
                candidate = cost + m
                if candidate <= U:
                    expected_transitions.add((i, n, b, c, m, candidate))
                    old = expected_best[i+1].get(m)
                    if old is None or candidate < old:
                        expected_best[i+1][m] = candidate

    actual_transitions = set()
    for t in transitions:
        item = (t["layer"], t["from"], t["b"], t["c"], t["to"], t["candidate_cost"])
        assert item not in actual_transitions, "Duplicate transition."
        actual_transitions.add(item)
    assert actual_transitions == expected_transitions, (
        f"Transition mismatch: expected {len(expected_transitions)}, "
        f"found {len(actual_transitions)}."
    )

    # The listed states must be exactly the states reached by the bounded
    # recurrence, with the exact minimum prefix costs.
    for i in range(D + 1):
        assert labels[i] == expected_best[i], f"State/cost mismatch at layer {i}."

    # Check every stored predecessor is a valid attaining transition.
    for i in range(1, D + 1):
        for n, cost in labels[i].items():
            p = state_records[(i, n)]["predecessor"]
            assert p is not None
            prev, b, c = p["from"], p["b"], p["c"]
            assert prev in labels[i-1]
            assert (i-1, prev, b, c, n, labels[i-1][prev] + n) in actual_transitions
            assert labels[i-1][prev] + n == cost
    assert state_records[(0, 1)]["predecessor"] is None

    # Terminal optimality within the bound, and witness attainability.
    assert labels[D].get(1) == U, "Terminal singleton does not have cost U."
    witness = cert["witness"]
    check_witness(D, witness)
    assert witness["order"] == U
    assert sum(witness["layer_sizes"]) == U
    assert witness["layer_sizes"][-1] == 1

    # The bounded Bellman recurrence is complete by induction over layers:
    # all source states are checked, and every allowed successor <= U is
    # included above. Therefore no complete path of order < U is omitted.
    assert cert["state_count"] == len(states)
    assert cert["transition_count"] == len(transitions)
    print(f"PASS: D={D}, M(D)={U:,}")
    print(f"PASS: {len(states)} states and {len(transitions)} transitions verified.")
    print("PASS: all degree products, balance equations and capacity constraints.")
    print("PASS: bounded-search completeness and feasible optimal witness.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("certificate", type=Path)
    args = parser.parse_args()
    verify(args.certificate)


if __name__ == "__main__":
    main()
