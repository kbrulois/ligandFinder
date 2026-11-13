"""After rebuilding aa_onehot_msacons.parquet with more MSAs: only the
precursors that previously had NO alignment (all-zero channel) may have changed.
Anything else means the scaling or the join moved, and the scan would present the
trunks a channel they were not trained on.  Exit 1 refuses the scan."""
import os
import sys

import numpy as np
import pyarrow.parquet as pq

A2 = os.path.expanduser("~/AF2_analysis")
old = pq.read_table(f"{A2}/lf_plm/aa_onehot_msacons_trainbuild.parquet",
                    columns=["accession", "plm_cons_msa"])
new = pq.read_table(f"{A2}/lf_plm/aa_onehot_msacons.parquet",
                    columns=["accession", "plm_cons_msa"])
o = old.column("plm_cons_msa").to_numpy()
n = new.column("plm_cons_msa").to_numpy()
if o.shape != n.shape:
    sys.exit(f"row count changed: {o.shape} -> {n.shape}")
acc = old.column("accession").to_numpy(zero_copy_only=False).astype(str)
chg = o != n
genes = sorted(set(acc[chg].tolist()))
bad = [g for g in genes if o[acc == g].max(initial=0.0) != 0.0]
print(f"  rows changed: {int(chg.sum()):,} in {len(genes)} precursors")
if bad:
    sys.exit(f"  {len(bad)} precursors changed that already HAD a conservation value, "
             f"e.g. {bad[:5]} -- refusing to scan")
print("  every changed precursor was all-zero before: OK")
