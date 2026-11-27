# lp_metrics: fused vs unfused (4 threads for the parallel column)

| dim | vertices | unfused ms | fused 1T ms | fused 4T ms | fused/unfused speedup | thread speedup | bit-identical across threads | fused complexity miss | unfused complexity miss |
|---|---|---|---|---|---|---|---|---|---|
| 2D | 200000 | 25.8 | 26.5 | 13.7 | 0.97x | 1.94x | yes | 2.8e-15 | 1.1e-8 |
| 3D | 200000 | 196.1 | 134.4 | 40.5 | 1.46x | 3.32x | yes | 2.2e-16 | 1.1e-5 |

* 2D: 200000 eigen-decompositions alone take 4.0 ms = 15% of the fused 1-thread time
* 3D: 200000 eigen-decompositions alone take 88.1 ms = 66% of the fused 1-thread time
