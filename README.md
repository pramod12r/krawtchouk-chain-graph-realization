# Krawtchouk chain exact-integer Python package

This package implements the mathematical specification in the supplied manuscript
"Graph Order Lower Bounds for One-Cell Equitable Realizations of the Krawtchouk Chain".

## Contents

- `krawtchouk_solver.py`: exact integer half-path dynamic programming (Theorem 11),
  witness construction, and bounded full-path Bellman certificate generation
  (Proposition 10).
- `krawtchouk_verifier.py`: independent certificate checker. It does not import
  the solver; it enumerates divisors using the backward degree `c` and independently
  reconstructs all bounded transitions and Bellman labels.
- `certificate_D16.json`, `certificate_D32.json`: generated state, cost,
  predecessor, transition and witness data.

## Requirements

Python 3.9 or later. No third-party packages are required.

## Reproduce the results

```bash
python krawtchouk_solver.py --dimension 16 --certificate certificate_D16.json
python krawtchouk_verifier.py certificate_D16.json

python krawtchouk_solver.py --dimension 32 --certificate certificate_D32.json
python krawtchouk_verifier.py certificate_D32.json
```

Expected key results:

- `D=16`, `M(16)=8,874`
- `D=32`, `M(32)=216,096,168`
- D=32 full-path certificate: 270 states and 817 permissible transitions.
- D=32 half-path DP: 6,144 reachable middle-layer sizes.

All calculations use Python integers. The enormous realizing graph is never
constructed.

## Mathematical conditions checked

For each edge between layers `i` and `i+1`, with
`q_i = (i+1)*(D-i)`:

1. `b_i*c_(i+1) = q_i`
2. `n_i*b_i = n_(i+1)*c_(i+1)`
3. `b_i <= n_(i+1)`
4. `c_(i+1) <= n_i`

The endpoints satisfy `n_0 = n_D = 1`. These are the conditions in
Proposition 3 of the manuscript and are necessary and sufficient for the
specified one-cell simple-graph realization.

## Certificate interpretation

The full-path certificate retains the minimum prefix cost for each layer-size
state, subject to the bound `U` equal to the feasible witness order. The
checker independently regenerates all permissible transitions from every
listed state and confirms that no bounded transition is omitted. It then checks
the exact Bellman labels, predecessor chains, terminal singleton cost, and
the feasible witness. The manuscript's arithmetic proof in Theorem 16 remains
a separate proof of the two target minima; the computational certificate is
an independent finite verification route.

## Scope

The result is the minimum within the one-cell equitable Krawtchouk realization
class. It is not a minimum-order result for all perfect-state-transfer graphs.
