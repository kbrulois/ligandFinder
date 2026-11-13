# Every page against the independent scans: the popup payload's CNN score, from
# the surface scored live on the page, must equal the per-gene scan's all-50 mean.
import json,os,re,sys,glob,csv
import numpy as np
A2=os.path.expanduser("~/AF2_analysis"); O=f"{A2}/pepend_pages_allsec_arms2"
genes=[l.strip() for f in ("genes.txt","giant_genes.txt") for l in open(f"{O}/{f}") if l.strip()]
giants={l.strip() for l in open(f"{O}/giant_genes.txt") if l.strip()}
arms={"AAswap":"curated_AAswap","ESM-C only":"esmonly36"}
Z={(lb,T):np.load(f"{A2}/lf_winxgb_scan_{T}_cv5all50_nomask_flatten_{a}_allsec_arms2.npz",allow_pickle=True) for lb,a in arms.items() for T in "NC"}
bad=[]; worst=0.0; nser={}; checked=0; total_anchors=0; acc={}
for k,g in enumerate(genes):
    p=f"{O}/{g}.html"
    try:
        with open(p,"rb") as fh:
            pay=None
            for line in fh:
                if line.startswith(b"<script>window.LF_PEP = "):
                    pay=line[len(b"<script>window.LF_PEP = "):].rstrip()
                    pay=pay[:-1] if pay.endswith(b";") else pay
                    break
            fh.seek(-16,2); tail=fh.read()
        if pay is None: bad.append((g,"no popup payload")); continue
        if b"</html>" not in tail: bad.append((g,"truncated")); continue
        d=json.loads(pay)
    except Exception as e:
        bad.append((g,f"read error {e}")); continue
    S=d["series"]; nser[len(S)]=nser.get(len(S),0)+1; acc[g]=d["accession"]
    want={f"{lb} all 50, {T} terminus" for lb in arms for T in "NC"}
    if set(S)!=want: bad.append((g,"series missing: "+", ".join(sorted(want-set(S)))))
    for key,v in S.items():
        lb,T=key.split(" all 50, ")[0],key[-10]
        z=Z[(lb,T)]; an=z[f"{g}__anchor"]; cnn=z[f"{g}__cnn"]
        if str(z[f"{g}__accession"])!=d["accession"]: bad.append((g,f"accession {key}")); continue
        idx=np.searchsorted(an,np.asarray(v["anchor"])); 
        if (idx>=len(an)).any() or not np.array_equal(an[idx],np.asarray(v["anchor"])): bad.append((g,f"anchors {key}")); continue
        e=float(np.abs(cnn[idx]-np.asarray(v["g"],dtype="f8")).max()); worst=max(worst,e); checked+=1; total_anchors+=len(v["anchor"])
        if e>1.5e-4: bad.append((g,f"{key}: max|g-scan|={e:.4g}"))
        cap=500 if g in giants else 2000
        if g in ("ANO8","C3","C5","EGF","COPA"): cap=600
        if len(v["anchor"])!=min(cap,len(an)): bad.append((g,f"{key}: {len(v['anchor'])} anchors, expected {min(cap,len(an))}"))
    if k%1000==0: print(k,g,flush=True)
print("pages read",len(genes),"series per page",nser); print("surfaces compared",checked,"anchors",total_anchors,"worst |page-scan|",worst)
print("problems",len(bad)); [print(" ",b) for b in bad[:40]]
json.dump(acc,open(sys.argv[1],"w"))
