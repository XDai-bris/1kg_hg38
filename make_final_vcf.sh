#!/usr/bin/env bash
set -euo pipefail
export PS1="${PS1-}"

# ---- Load Python on BlueCrystal (safe) ----
set +u
if command -v module >/dev/null 2>&1; then
  source /etc/profile 2>/dev/null || true
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh || true
  module load languages/python/3.12.3
fi
set -u
command -v python3 >/dev/null 2>&1 || { echo "❌ python3 not found"; exit 127; }
echo "🐍 Using $(python3 --version) at $(which python3)"

# ---- Inputs ----
POS_FILE="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq/1kg_hg38_pos.txt"
ALL_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome_ALL/ALL_genome_freq.afreq"
AFR_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome/AFR_genome_freq.afreq"
AMR_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome/AMR_genome_freq.afreq"
EAS_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome/EAS_genome_freq.afreq"
EUR_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome/EUR_genome_freq.afreq"
SAS_AFREQ="/user/home/xd14188/repo/1kg_hg38/merged_genome/SAS_genome_freq.afreq"

# ---- Outputs ----
OUT_DIR="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq"
VCF_OUT="${OUT_DIR}/with_popAF.vcf"
LOG_OUT="${OUT_DIR}/with_popAF.log"
QC_OUT="${OUT_DIR}/with_popAF_qc.tsv"
mkdir -p "${OUT_DIR}"

# ---- Sanity checks ----
for f in "$POS_FILE" "$ALL_AFREQ" "$AFR_AFREQ" "$AMR_AFREQ" "$EAS_AFREQ" "$EUR_AFREQ" "$SAS_AFREQ"; do
  [[ -s "$f" ]] || { echo "❌ Missing or empty: $f"; exit 2; }
done

echo "🚀 Building VCF for all rows in ${POS_FILE} …"

python3 - << 'PYCODE'
import sys, gzip

POS_FILE="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq/1kg_hg38_pos.txt"
VCF_OUT ="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq/with_popAF.vcf"
LOG_OUT ="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq/with_popAF.log"
QC_OUT  ="/user/home/xd14188/repo/1kg_hg38/final_vcf_outputWithAfreq/with_popAF_qc.tsv"

AF_FILES={
 "ALL":"/user/home/xd14188/repo/1kg_hg38/merged_genome_ALL/ALL_genome_freq.afreq",
 "AFR":"/user/home/xd14188/repo/1kg_hg38/merged_genome/AFR_genome_freq.afreq",
 "AMR":"/user/home/xd14188/repo/1kg_hg38/merged_genome/AMR_genome_freq.afreq",
 "EAS":"/user/home/xd14188/repo/1kg_hg38/merged_genome/EAS_genome_freq.afreq",
 "EUR":"/user/home/xd14188/repo/1kg_hg38/merged_genome/EUR_genome_freq.afreq",
 "SAS":"/user/home/xd14188/repo/1kg_hg38/merged_genome/SAS_genome_freq.afreq",
}

def open_auto(p):
    return gzip.open(p,"rt") if p.endswith(".gz") else open(p,"rt")

def load_afreq(p):
    """
    Parse plink2 --freq output.
    Accept ALT_FREQS and similar names.
    Always take first ALT and first ALT_FREQS value.
    Store only REF/ALT/AF to save RAM.
    """
    d={}
    with open_auto(p) as f:
        # Header detection
        header = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                if line.startswith("##"):
                    continue
                header=line.lstrip("#").strip().split()
                break
            else:
                header=line.strip().split()
                break

        if header is None:
            raise RuntimeError(f"Could not find header in {p}")

        cols={h:i for i,h in enumerate(header)}
        # fallback lowercase
        for k in ("CHROM","ID","REF","ALT"):
            if k not in cols and k.lower() in cols:
                cols[k]=cols[k.lower()]

        # accept multiple AF column names
        af_key=None
        for cand in ("ALT_FREQS","ALT_FREQ","A1_FREQ","ALT_AF","ALTFRQ","FREQ","AF"):
            if cand in cols:
                af_key=cand
                break
        if af_key is None:
            raise RuntimeError(f"No AF-like column in {p} (header={header})")

        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts=line.split()
            if len(parts) < len(header):
                continue
            vid=parts[cols["ID"]]
            ref=parts[cols["REF"]]
            alt=parts[cols["ALT"]]
            afs=parts[cols[af_key]]

            alt_first = alt.split(",")[0]
            af_first  = afs.split(",")[0]

            d[vid] = {
                "REF": ref,
                "ALT": alt_first,
                "AF":  af_first,
            }
    return d

# Load AF dicts in memory (one per pop)
af = {pop: load_afreq(path) for pop,path in AF_FILES.items()}

def ffmt_zero(x):
    """
    Format AF values:
    - Missing or unparsable → 0
    - Otherwise → float formatting
    """
    if x is None:
        return "0"
    if isinstance(x,str):
        if x in ("", ".", "NA"):
            return "0"
        try:
            v=float(x)
        except Exception:
            return "0"
    else:
        try:
            v=float(x)
        except Exception:
            return "0"
    s=f"{v:.6f}"
    return s.rstrip("0").rstrip(".") if "." in s else s

vcf_header=[
 "##fileformat=VCFv4.2",
 "##source=make_final_vcf",
 '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency in ALL">',
 '##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="ALT allele frequency in AFR">',
 '##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="ALT allele frequency in AMR">',
 '##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="ALT allele frequency in EAS">',
 '##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="ALT allele frequency in EUR">',
 '##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="ALT allele frequency in SAS">',
 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
]

# Stats + QC tracking
n_variants=0
stats = {
    pop: {"present":0, "missing":0, "exact":0, "swapped":0, "mismatch":0}
    for pop in AF_FILES.keys()
}
swap_notes=[]
hard_mismatches=[]

with open(VCF_OUT,"w") as vfh, open(POS_FILE) as posfh, open(QC_OUT,"w") as qcfh:
    # write VCF header
    vfh.write("\n".join(vcf_header) + "\n")

    # write QC header
    qcfh.write("\t".join([
        "ID","POP","ISSUE_TYPE",
        "POS_REF","POS_ALT",
        "AFREF_REF","AFREF_ALT",
        "RAW_AF","USED_AF"
    ]) + "\n")

    # skip header of POS file
    pos_header = posfh.readline()

    for line in posfh:
        if not line.strip():
            continue
        chrom,pos,ref,alt,vid = line.split()
        n_variants += 1

        info_vals = {}

        for pop,key in [
            ("ALL","AF"),
            ("AFR","AFR_AF"),
            ("AMR","AMR_AF"),
            ("EAS","EAS_AF"),
            ("EUR","EUR_AF"),
            ("SAS","SAS_AF"),
        ]:
            rec = af[pop].get(vid)
            if rec is None:
                # not present in this afreq → treat as AF=0
                info_vals[key] = "0"
                stats[pop]["missing"] += 1
                continue

            stats[pop]["present"] += 1
            rec_ref = rec["REF"]
            rec_alt = rec["ALT"]
            raw_af  = rec["AF"]

            # Try to parse AF once
            try:
                af_val = float(raw_af)
            except Exception:
                af_val = None

            # Case 1: exact match REF/ALT → keep AF as is
            if rec_ref == ref and rec_alt == alt:
                stats[pop]["exact"] += 1
                used_af = ffmt_zero(af_val)
                info_vals[key] = used_af

            # Case 2: swapped alleles (same set, reversed order) → flip AF: 1 - AF
            elif rec_ref == alt and rec_alt == ref:
                stats[pop]["swapped"] += 1
                if af_val is not None:
                    flipped = 1.0 - af_val
                    used_af = ffmt_zero(flipped)
                else:
                    used_af = "0"
                info_vals[key] = used_af
                swap_notes.append(
                    f"{vid}: {pop} swapped REF/ALT ({rec_ref}/{rec_alt} vs {ref}/{alt}), AF flipped"
                )
                qcfh.write("\t".join([
                    vid, pop, "SWAPPED_AF_FLIPPED",
                    ref, alt,
                    rec_ref, rec_alt,
                    raw_af, used_af
                ]) + "\n")

            # Case 3: genuine mismatch (different allele set) → set 0 and warn
            else:
                stats[pop]["mismatch"] += 1
                used_af = "0"
                info_vals[key] = used_af
                msg = f"{vid}: {pop} allele mismatch: {rec_ref}/{rec_alt} vs {ref}/{alt}"
                hard_mismatches.append(msg)
                qcfh.write("\t".join([
                    vid, pop, "HARD_MISMATCH_AF_ZERO",
                    ref, alt,
                    rec_ref, rec_alt,
                    raw_af, used_af
                ]) + "\n")

        # fill any missing keys with "0" (should not happen, but just in case)
        for _, key in [("ALL","AF"),("AFR","AFR_AF"),("AMR","AMR_AF"),
                       ("EAS","EAS_AF"),("EUR","EUR_AF"),("SAS","SAS_AF")]:
            info_vals.setdefault(key, "0")

        info = (
            f"AF={info_vals['AF']};"
            f"AFR_AF={info_vals['AFR_AF']};"
            f"AMR_AF={info_vals['AMR_AF']};"
            f"EAS_AF={info_vals['EAS_AF']};"
            f"EUR_AF={info_vals['EUR_AF']};"
            f"SAS_AF={info_vals['SAS_AF']}"
        )

        vfh.write("\t".join([chrom,pos,vid,ref,alt,"100","PASS",info]) + "\n")

# Write log summary
with open(LOG_OUT,"w") as lfh:
    lfh.write(f"Wrote VCF for {n_variants} variants\n")
    lfh.write(f"QC report (swaps + hard mismatches): {QC_OUT}\n\n")
    lfh.write("Per-population stats (counts of variants by AF status):\n")
    lfh.write("POP\tPRESENT_IN_AFREQ\tMISSING_IN_AFREQ\tEXACT_MATCH\tSWAPPED_AF_FLIPPED\tHARD_MISMATCH_AF_ZERO\n")
    for pop in AF_FILES.keys():
        s = stats[pop]
        lfh.write(
            f"{pop}\t{s['present']}\t{s['missing']}\t"
            f"{s['exact']}\t{s['swapped']}\t{s['mismatch']}\n"
        )
    lfh.write("\n")

    if swap_notes:
        lfh.write("Swapped allele cases (AF flipped):\n")
        for w in swap_notes:
            lfh.write("  - " + w + "\n")
        lfh.write("\n")

    if hard_mismatches:
        lfh.write("Hard mismatches (allele sets differ, AF set to 0):\n")
        for w in hard_mismatches:
            lfh.write("  - " + w + "\n")
    else:
        lfh.write("No hard mismatches detected.\n")

print(f"✅ Output VCF: {VCF_OUT}")
print(f"📝 Log:        {LOG_OUT}")
print(f"🧾 QC report:  {QC_OUT}")
PYCODE

echo "✅ Done. Files:"
ls -l "${VCF_OUT}" "${LOG_OUT}" || true


# ---- Compress and index final VCF to BCF ----
echo "📦 Compressing and indexing VCF to BCF..."
module load bcftools/1.19-openblas-2333
module load  htslib/1.19.1-v2zi
# Write header (all lines starting with #)
grep '^#' with_popAF.vcf > with_popAF.sorted.vcf

# Append sorted body (non-header)
grep -v '^#' with_popAF.vcf \
  | LC_ALL=C sort -T /tmp -k1,1V -k2,2n \
  >> with_popAF.sorted.vcf
head with_popAF.sorted.vcf
tail with_popAF.sorted.vcf
wc -l with_popAF.vcf with_popAF.sorted.vcf
bgzip with_popAF.sorted.vcf
file with_popAF.sorted.vcf.gz
zcat with_popAF.sorted.vcf.gz | head
tabix -p vcf with_popAF.sorted.vcf.gz

bcftools view -Ob -o with_popAF.bcf with_popAF.sorted.vcf.gz
bcftools index -f with_popAF.bcf

ls -lh with_popAF.bcf
bcftools view -h with_popAF.bcf | head
bcftools view with_popAF.bcf | head