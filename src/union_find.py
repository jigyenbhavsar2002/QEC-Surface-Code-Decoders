"""Union-find decoder (Delfosse & Nickerson, 2021), unweighted-growth version.

Works on the decoding graph of a detector error model:
  nodes = detectors (+ one virtual boundary node), edges = error mechanisms.

Algorithm
  1. Every defect (fired detector) starts as its own odd cluster.
  2. Repeat until no odd cluster remains: every odd cluster that does not touch the
     boundary grows by half an edge in every direction. Edges that become fully grown
     merge the clusters they connect (union-find with small-to-large merging).
     A cluster is "neutral" when its defect count is even or it contains the boundary.
  3. Peeling: build a spanning forest of each cluster (rooted at the boundary if the
     cluster touches it) and peel leaves, adding an edge to the correction whenever
     the leaf still carries a defect. The correction's logical effect is the XOR of
     the observable masks of the chosen edges.

Edge weights are ignored during growth (all edges have equal length). This is what
makes the decoder nearly linear-time, and also what costs it accuracy vs MWPM.
"""
import numpy as np
import pymatching
import stim


class UnionFindDecoder:
    def __init__(self, matching: pymatching.Matching):
        self.num_detectors = matching.num_detectors
        self.num_obs = matching.num_fault_ids
        self.B = self.num_detectors  # virtual boundary node
        self.eu, self.ev, self.emask = [], [], []
        self.adj = [[] for _ in range(self.num_detectors + 1)]
        for u, v, attr in matching.edges():
            u = self.B if u is None else u
            v = self.B if v is None else v
            mask = 0
            for f in attr.get("fault_ids", ()):
                mask |= 1 << f
            e = len(self.eu)
            self.eu.append(u); self.ev.append(v); self.emask.append(mask)
            self.adj[u].append((e, v))
            self.adj[v].append((e, u))

    @classmethod
    def from_circuit(cls, circuit: stim.Circuit):
        dem = circuit.detector_error_model(decompose_errors=True)
        return cls(pymatching.Matching.from_detector_error_model(dem))

    # ------------------------------------------------------------------ decoding
    def decode_mask(self, defects) -> int:
        """Return the predicted logical-observable flips as a bitmask."""
        defects = [int(x) for x in defects]
        if not defects:
            return 0
        B, adj, eu, ev, emask = self.B, self.adj, self.eu, self.ev, self.emask
        parent, members, parity, has_b = {}, {}, {}, {}
        roots = set(defects)

        def add(n):
            parent[n] = n; members[n] = [n]; parity[n] = 0; has_b[n] = (n == B)

        def find(n):
            while parent[n] != n:
                parent[n] = parent[parent[n]]
                n = parent[n]
            return n

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra == rb:
                return
            if len(members[ra]) < len(members[rb]):
                ra, rb = rb, ra
            parent[rb] = ra
            members[ra].extend(members.pop(rb))
            parity[ra] ^= parity[rb]
            has_b[ra] = has_b[ra] or has_b[rb]
            roots.discard(rb)

        for n in defects:
            add(n); parity[n] = 1

        # ---- 1. cluster growth
        growth, full_edges = {}, []
        while True:
            active = [r for r in roots if parity[r] and not has_b[r]]
            if not active:
                break
            newly_full, grew = [], False
            for r in active:
                for n in members[r]:
                    for e, _ in adj[n]:
                        g = growth.get(e, 0)
                        if g < 2:
                            growth[e] = g + 1
                            grew = True
                            if g == 1:
                                newly_full.append(e)
            if not grew:  # isolated odd cluster with nowhere to grow (malformed graph)
                break
            for e in newly_full:
                u, v = eu[e], ev[e]
                for x in (u, v):
                    if x not in parent:
                        add(x); roots.add(x)
                full_edges.append(e)
                union(u, v)

        # ---- 2. peeling on a spanning forest of the fully grown edges
        fadj = {}
        for e in full_edges:
            fadj.setdefault(eu[e], []).append((e, ev[e]))
            fadj.setdefault(ev[e], []).append((e, eu[e]))
        syndrome = set(defects)
        visited, correction = set(), 0
        starts = ([B] if B in fadj else []) + [n for n in fadj if n != B]
        for s in starts:
            if s in visited:
                continue
            visited.add(s)
            order, par = [s], {}
            i = 0
            while i < len(order):
                n = order[i]; i += 1
                for e, m in fadj[n]:
                    if m not in visited:
                        visited.add(m); par[m] = (n, e); order.append(m)
            for n in reversed(order[1:]):
                if n in syndrome:
                    p, e = par[n]
                    correction ^= emask[e]
                    syndrome.discard(n)
                    if p != B:
                        syndrome ^= {p}
        return correction

    def decode_batch(self, detections: np.ndarray) -> np.ndarray:
        """detections: (shots, num_detectors) bool -> (shots, num_obs) bool."""
        out = np.zeros((detections.shape[0], self.num_obs), dtype=bool)
        for i, row in enumerate(detections):
            m = self.decode_mask(np.flatnonzero(row))
            for k in range(self.num_obs):
                out[i, k] = (m >> k) & 1
        return out
