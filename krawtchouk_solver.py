#!/usr/bin/env python3
"""Exact-integer solver and bounded Bellman-certificate generator.

Implements Proposition 10 and Theorem 11 of the supplied manuscript.
No realizing graph is constructed.
"""
from __future__ import annotations
import argparse
import json
from math import isqrt
from pathlib import Path


def positive_divisors(x: int) -> list[int]:
    """Return all positive divisors of x, in increasing order."""
    small, large = [], []
    for d in range(1, isqrt(x) + 1):
        if x % d == 0:
            small.append(d)
            if d * d != x:
                large.append(x // d)
    return small + large[::-1]


def transitions_from_b(i: int, D: int, n: int):
    """Generator-side transitions: enumerate b | q_i, c=q_i/b."""
    q = (i + 1) * (D - i)
    for b in positive_divisors(q):
        c = q // b
        if c > n:                       # capacity c <= n
            continue
        numerator = n * b
        if numerator % c:
            continue                    # m = nb/c must be integral
        m = numerator // c
        if b > m:                       # capacity b <= m
            continue
        yield b, c, m


def solve_even_half_path(D: int):
    """Theorem 11: exact half-path DP for even D."""
    if D < 2 or D % 2:
        raise ValueError("The half-path reduction requires an even D >= 2.")
    h = D // 2
    # cost excludes the middle layer; each pre-middle layer is counted twice.
    costs = {1: 0}
    predecessors = [{1: None}]
    for i in range(h):
        next_costs, next_pred = {}, {}
        for n, old_cost in costs.items():
            for b, c, m in transitions_from_b(i, D, n):
                new_cost = old_cost + 2 * n
                if new_cost < next_costs.get(m, 10**100):
                    next_costs[m] = new_cost
                    next_pred[m] = (n, b, c)
        costs = next_costs
        predecessors.append(next_pred)

    total, middle = min((cost + n, n) for n, cost in costs.items())
    first_half = [middle]
    half_degrees = []
    current = middle
    for layer in range(h, 0, -1):
        previous, b, c = predecessors[layer][current]
        half_degrees.append((b, c))
        first_half.append(previous)
        current = previous
    first_half.reverse()
    half_degrees.reverse()

    # Reflect the first half; the central layer is included only once.
    layers = first_half + first_half[-2::-1]
    forward, backward = [], []
    for b, c in half_degrees:
        forward.append(b)
        backward.append(c)
    for b, c in reversed(half_degrees):
        # Reversal swaps forward and backward degrees.
        forward.append(c)
        backward.append(b)

    assert len(layers) == D + 1
    assert sum(layers) == total
    validate_witness(D, layers, forward, backward)
    return {
        "D": D, "minimum_order": total, "layer_sizes": layers,
        "forward_degrees": forward, "backward_degrees": backward,
        "middle_states": len(costs),
    }


def validate_witness(D, n, b, c):
    if len(n) != D + 1 or len(b) != D or len(c) != D:
        raise ValueError("Witness has incorrect list lengths.")
    if n[0] != 1 or n[D] != 1:
        raise ValueError("Both endpoint layers must be singletons.")
    for i in range(D):
        q = (i + 1) * (D - i)
        if min(n[i], n[i+1], b[i], c[i]) < 1:
            raise ValueError(f"Nonpositive value at edge {i}.")
        if b[i] * c[i] != q:
            raise ValueError(f"Quotient product fails at edge {i}.")
        if n[i] * b[i] != n[i+1] * c[i]:
            raise ValueError(f"Balance fails at edge {i}.")
        if b[i] > n[i+1] or c[i] > n[i]:
            raise ValueError(f"Simple-graph capacity fails at edge {i}.")


def build_full_certificate(D: int, U: int):
    """Full directed-path Bellman search, retaining prefixes of cost <= U."""
    labels = [{1: 1}]
    pred = [{1: None}]
    transitions = []

    for i in range(D):
        next_labels, next_pred = {}, {}
        layer_transitions = []
        for n, prefix_cost in labels[-1].items():
            for b, c, m in transitions_from_b(i, D, n):
                new_cost = prefix_cost + m
                if new_cost <= U:
                    layer_transitions.append({
                        "layer": i, "from": n, "b": b, "c": c,
                        "to": m, "candidate_cost": new_cost
                    })
                    if new_cost < next_labels.get(m, 10**100):
                        next_labels[m] = new_cost
                        next_pred[m] = {
                            "from": n, "b": b, "c": c
                        }
        transitions.extend(layer_transitions)
        labels.append(next_labels)
        pred.append(next_pred)

    states = []
    for i, layer in enumerate(labels):
        for n, cost in sorted(layer.items()):
            states.append({
                "layer": i, "n": n, "cost": cost,
                "predecessor": pred[i][n]
            })

    if 1 not in labels[D]:
        raise RuntimeError(f"No terminal singleton state found within bound U={U}.")

    # Recover a minimum-cost full path.
    layer_sizes = [None] * (D + 1)
    forward, backward = [None] * D, [None] * D
    current = 1
    layer_sizes[D] = current
    for i in range(D, 0, -1):
        p = pred[i][current]
        if p is None:
            raise RuntimeError("Broken predecessor chain.")
        layer_sizes[i-1] = p["from"]
        forward[i-1], backward[i-1] = p["b"], p["c"]
        current = p["from"]

    validate_witness(D, layer_sizes, forward, backward)
    if sum(layer_sizes) != labels[D][1]:
        raise RuntimeError("Recovered witness cost mismatch.")

    return {
        "dimension": D,
        "bound_U": U,
        "states": states,
        "transitions": transitions,
        "witness": {
            "layer_sizes": layer_sizes,
            "forward_degrees": forward,
            "backward_degrees": backward,
            "order": sum(layer_sizes)
        },
        "terminal_cost": labels[D][1],
        "state_count": len(states),
        "transition_count": len(transitions),
        "layer_state_counts": [len(x) for x in labels],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=32)
    parser.add_argument("--certificate", type=Path, default=None)
    args = parser.parse_args()

    result = solve_even_half_path(args.dimension)
    print(f"D={result['D']}")
    print(f"M(D)={result['minimum_order']:,}")
    print("Layer sizes:", result["layer_sizes"])
    print("Forward degrees:", result["forward_degrees"])
    print("Backward degrees:", result["backward_degrees"])
    print(f"Reachable middle-layer states: {result['middle_states']:,}")

    if args.certificate:
        cert = build_full_certificate(args.dimension, result["minimum_order"])
        args.certificate.write_text(json.dumps(cert, indent=2), encoding="utf-8")
        print(f"Certificate written to {args.certificate}")
        print(f"Full-path states: {cert['state_count']}")
        print(f"Permissible transitions: {cert['transition_count']}")


if __name__ == "__main__":
    main()
