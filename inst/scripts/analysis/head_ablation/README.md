# Attention head vs average pooling (2026-09-30)

**Asked:** on the peptide-end C set, is the attention global head better than
plain average pooling ("gap")?

**Found:** no. Over 10 seeds each the two are indistinguishable.
- PR AUC p = 0.80.
- The two retrieval metrics disagree in direction, and per peptide it is 16-13.
- gap is 383 parameters smaller and more stable across seeds (ROC sd 0.005 vs 0.013).

`10_9l_head_ablation.R`
