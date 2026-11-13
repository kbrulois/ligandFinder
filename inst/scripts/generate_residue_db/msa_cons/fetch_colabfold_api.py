"""Fetch ColabFold MMseqs2 MSAs from the public API for precursors the LiteFold
human-proteome set does not hold (immunoglobulin segments, UniProt sequence
revisions, ...), so every precursor's alignment comes from the same search:
MMseqs2 against UniRef30 + the ColabFold environmental DB, `mode=env`.

The result tar holds `uniref.a3m` and `bfd.mgnify30.metaeuk30.smag30.a3m`; they
are concatenated the way ColabFold itself does (env hits appended after the
UniRef hits, its duplicate query record dropped), which is the layout of the
LiteFold files.

Usage: python fetch_colabfold_api.py <api_needed.tsv> <outdir> [--workers 3]
Writes <outdir>/<accession>.a3m.gz; finished files are skipped, so re-running
resumes. Polite by design: few concurrent tickets, backoff on RATELIMIT.
"""
import argparse, csv, gzip, io, json, os, sys, tarfile, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

HOST = "https://api.colabfold.com"
UA = {"User-Agent": "ligandFinder-msa-conservation/1.0"}


def _req(url, data=None, timeout=120, tries=6):
    """One request, retried on any network error. A bare urlopen here killed a
    380/388 run on a single `[Errno 60] Operation timed out` while polling."""
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    for k in range(tries):
        try:
            r = urllib.request.Request(url, data=body, headers=UA)
            with urllib.request.urlopen(r, timeout=timeout) as f:
                return f.read()
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"  {type(e).__name__} on {url.rsplit('/', 1)[-1]}; retry {k + 1}/{tries - 1}",
                  file=sys.stderr, flush=True)
            time.sleep(30 * (k + 1))


def submit(seq):
    while True:
        try:
            out = json.loads(_req(f"{HOST}/ticket/msa", {"q": f">101\n{seq}\n", "mode": "env"}))
        except Exception as e:  # network hiccup
            print(f"  submit error {e}; retry in 60s", file=sys.stderr); time.sleep(60); continue
        if out.get("status") in ("RATELIMIT", "UNKNOWN"):
            time.sleep(120); continue
        return out


def run_one(acc, seq, outdir):
    dest = os.path.join(outdir, f"{acc}.a3m.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return acc, "skip"
    t = submit(seq)
    waited = 0
    while t.get("status") in ("PENDING", "RUNNING"):
        time.sleep(15); waited += 15
        if waited > 7200:
            return acc, "FAILED ticket still running after 2 h"
        try:
            t = json.loads(_req(f"{HOST}/ticket/{t['id']}"))
        except Exception as e:
            return acc, f"FAILED polling: {type(e).__name__}"
    if t.get("status") != "COMPLETE":
        return acc, f"FAILED status={t.get('status')}"
    raw = _req(f"{HOST}/result/download/{t['id']}", timeout=1800)
    tf = tarfile.open(fileobj=io.BytesIO(raw))
    names = {os.path.basename(m.name): m for m in tf.getmembers()}
    parts = []
    for i, nm in enumerate(["uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m"]):
        if nm not in names:
            continue
        txt = tf.extractfile(names[nm]).read().decode().replace("\x00", "")
        recs = txt.strip("\n").split("\n")
        if i > 0 and recs and recs[0].startswith(">101"):
            recs = recs[2:]  # the env file repeats the query
        parts.append("\n".join(recs))
    a3m = "\n".join(p for p in parts if p) + "\n"
    # the query must come back unchanged, or the alignment columns are not our residues
    q = a3m.split("\n")[1]
    if q != seq:
        return acc, "FAILED query mismatch"
    tmp = dest + ".tmp"
    with gzip.open(tmp, "wt") as g:
        g.write(a3m)
    os.replace(tmp, dest)
    return acc, f"ok {a3m.count(chr(62))} seqs"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv"); ap.add_argument("outdir")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    with open(a.tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    rows.sort(key=lambda r: len(r["seq"]))  # short first: most of the set finishes early
    n_fail = 0
    with ThreadPoolExecutor(a.workers) as ex:
        futs = {ex.submit(run_one, r["accession"], r["seq"], a.outdir): r["accession"] for r in rows}
        for i, fu in enumerate(as_completed(futs), 1):
            try:
                acc, msg = fu.result()
            except Exception as e:           # never let one precursor end the run
                acc, msg = futs[fu], f"FAILED {type(e).__name__}: {e}"
            n_fail += msg.startswith("FAILED")
            print(f"[{i}/{len(rows)}] {acc} {msg}", flush=True)
    print(f"done; {n_fail} failed")


if __name__ == "__main__":
    main()
