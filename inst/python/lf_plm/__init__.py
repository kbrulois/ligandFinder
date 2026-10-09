"""lf_plm -- per-residue protein language model embeddings as lf_dcnn channels.

Two steps, both resumable and runnable with no R in the loop:

    python -m lf_plm embed  --sequences seqs.parquet --out plm_raw.h5
    python -m lf_plm reduce --h5 plm_raw.h5 --k 32 --out plm_residues.parquet

``embed`` writes one ``(L, D)`` float16 array per accession (ProtT5-XL-U50
encoder by default, 1024-d; ESM C as an alternative).  ``reduce`` fits one
PCA over a residue sample of the whole set and writes a long table
``(accession, index, AA, plm_01..plm_k)`` scaled to ``[0, 1]`` -- the join keys
``9.2_add_contact_data.R`` / ``R/plm_features.R`` use to slot the columns in
beside the 26 hand-built channels.

The input is a parquet with ``accession`` and ``sequence_uni`` columns (what R
exports from ``secretome``), not FASTA: the pipeline never has FASTA.
"""

from .embed import BACKENDS, embed_sequences
from .reduce import fit_pca, transform_to_table

__all__ = ["BACKENDS", "embed_sequences", "fit_pca", "transform_to_table"]
__version__ = "0.1.0"
