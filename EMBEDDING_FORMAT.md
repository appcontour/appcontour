# Embedding file format

Embedding files have the following form

```text
embedding { 0: (6, 7, 8), 1: (5, 8, 7), 2: (4, 6, 8, 8), 3: (4, 5, 7, 6), 4: odd(2, 5, 3, 6), 5: (1, 7, 3, 4), 6: odd(0, 2, 4, 3), 7: (0, 3, 5, 1), 8: (0, 1, 2, 2) }
```

The vertex labels, and neighbor labels must be integer literals, and
each adjacency row must be a tuple. The `odd` specification means that the overpass is
the one connecting the second and fourth vertices in the cyclic list `(x, x, x, x)`, otherwise
it is the one connecting the first and third vertices.

The over/under information can also be given as an integer following the `embedding` keyword, as in

```text
embedding:20 { 0: (6, 7, 8), 1: (5, 8, 7), 2: (4, 6, 8, 8), 3: (4, 5, 7, 6), 4: (2, 5, 3, 6), 5: (1, 7, 3, 4), 6: (0, 2, 4, 3), 7: (0, 3, 5, 1), 8: (0, 1, 2, 2) }
```

See the relevant section below.

Commas in the cyclic lists are optional.

## Underlying plane graph

The dictionary describes a combinatorial embedding of the underlying plane
graph. Vertices are numbered contiguously from `0`. The two trivalent vertices
are `0` and `1`, and the crossing vertices are the quadrivalent vertices
`2, 3, ...`.

For each vertex, the tuple lists the incident half-edges in counterclockwise
cyclic order. Repeated neighbor labels represent parallel half-edges.
For example, in

```text
2: (4, 6, 8, 8)
```

vertex `2` is incident to one half-edge ending at `4`, one at `6`, and two
distinct half-edges ending at `8`, in that counterclockwise order.

Parallel edges have no individual identifiers in this format. Their
occurrences must form cyclic blocks at both endpoints and are paired in reverse
order, which is the noncrossing pairing for counterclockwise adjacency lists.
This convention lets the checker identify the two darts of each edge
unambiguously.

## Crossing choice

The integer after `embedding:`, if given, is the over/under choice. If the diagram has
`c` crossings, interpret its binary expansion using `c` bits, padded with
leading zeroes if necessary. The bits are read from right to left:

```text
bit d = floor(choice / 2^d) mod 2
```

is the bit for crossing vertex `k + d` where `k` is the number of trivalent vertices.
Equivalently, the rightmost binary
digit is for vertex `k`, the next digit to the left is for vertex `k+1`, and so
on. The choice must satisfy `0 <= choice < 2^c`.

At a crossing vertex with counterclockwise adjacency tuple

```text
v: (a, b, c, d)
```

the bit determines which opposite pair is the overpass:

| bit | overpassing arc | underpassing arc |
| --- | --- | --- |
| `0` | `a-v-c` | `b-v-d` |
| `1` | `b-v-d` | `a-v-c` |

The trivalent vertices have no crossing bit.

For the example above, `20` is binary `10100`. With 7 crossing vertices
`2, ..., 8`, this is padded to `0010100`; reading from right to left gives
bits

```text
vertex: 2 3 4 5 6 7 8
bit:    0 0 1 0 1 0 0
```

Thus vertex `2`, whose tuple is `(4, 6, 8, 8)`, has bit `0`, so the arc
`4-2-8` through the first occurrence of `8` (at zero-based tuple index `2`) is
overpassing there. Vertex `4`, whose tuple is `(2, 5, 3, 6)`, has bit `1`, so
the arc `5-4-6` is overpassing there.

## Generic syntax

embedding[:c] [-]{0: (x, x, x), 1: (x, x, x), ..., k: (x, x, x), <k+1>: [e][o](x, x, x, x), <k+2>: [e][o](x, x, x, x),...}

where the optional `-` sign before `{` means that the cyclic lists are clockwise instead of counter-clockwise.
Indication `e` (or `even`) is optional in case the overcrossing connects the first to the third vertex in the cyclic list,
Indication `o` (or `odd`) means that the overcrossing connects the second to the fourth vertex of the cyclic list.

The overpass information can be given either locally for each quadrivalent vertex using prefix `o` if the overpass
involves the second and fourth entry in the cyclic list or globally using the `c` integer value.

The commas in the cyclic lists can be omitted.

The `choice` value can be overridden with the `appcontour` option `--choice <value>`

## The Hopf link conundrum

There is a situation where the overpass indication at each crossing is ambiguous.  This is exemplified in the following
example (the Hopf link)

```text
embedding {0 : (1 1 1 1), 1: (0 0 0 0)}
```

Here there is no way to establish which half-arc starts the cyclic list, resulting in uncertainty regarding the
overpass choice at the two vertices. Upon analysis one realizes that there are exactly two possible interpretations
(up to exchange of the two vertices): one gives the Hopf link description, the other is a situation where
a direct Reidemeister II move can be performed leading to a split link.

Conventionally we use a value of 0 or 3 for the choice value (both vertices are even or both are odd) to indicate
the case of the Hopf link whereas the values 1 or 2 correspond to the `fake` Hopf choice of overcrossings.

Note that this situation can appear inside a more complex embedding description when a split component of the
embedding involving two vertices link to each other.

## Presence of one or more `split` closed loops

A `split` closed loop is a close S1 with no vertices. The number of such trivial components can be added to an
embedding description via an additional `+<n>` before the closing brace, where <n> is the number of such trivial
components.
[Not implemented yet]
