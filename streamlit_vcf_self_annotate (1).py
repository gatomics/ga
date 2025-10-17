import os, gzip, subprocess, tempfile
import streamlit as st
import pandas as pd
from datetime import datetime

# ------------------ Auth (optional) ------------------
APP_PASSWORD = os.getenv("APP_PASSWORD")
if APP_PASSWORD:
    st.set_page_config(page_title="Genmatix/GATomics — VCF → Report", layout="wide")
    st.title("🔐 Login")
    pw = st.text_input("Password", type="password")
    if pw != APP_PASSWORD:
        st.stop()

st.set_page_config(page_title="Genmatix/GATomics — VCF → Report (Self-Annotating)", layout="wide")
st.markdown("<h1>🧬 Genmatix/GATomics — VCF → Report</h1><p><i>Self‑annotating (VEP) | Single-sample VCF</i></p>", unsafe_allow_html=True)

# ------------------ Inputs ------------------
col1, col2 = st.columns([2,1], gap="large")
with col1:
    vcf_file = st.file_uploader("Upload VCF (.vcf or .vcf.gz)", type=["vcf","gz"])
    sample_name = st.text_input("Sample name (optional)", "")
    hpo_input = st.text_area("HPO terms (optional)", placeholder="HP:0003774 HP:0000083")
with col2:
    genome_build = st.selectbox("Genome build (for VEP)", ["GRCh38","GRCh37"], index=0)
    impact_keep = st.multiselect("Impacts to keep", ["HIGH","MODERATE","LOW","MODIFIER"], default=["HIGH","MODERATE"])
    af_thresh = st.number_input("Max AF (keep ≤)", min_value=0.0, max_value=1.0, value=0.01, step=0.001)
    zyg_keep = st.multiselect("Zygosity filter", ["HET","HOM"], default=["HET","HOM"])
    panel_file = st.file_uploader("Gene panel boost (txt/csv, optional)", type=["txt","csv"])
run = st.button("⚡ Annotate & Generate Report")

IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

# ------------------ Helpers ------------------
def parse_info(info_str):
    out = {}
    for item in info_str.split(";"):
        if not item: 
            continue
        if "=" in item:
            k,v = item.split("=",1)
            out[k]=v
        else:
            out[item]=True
    return out

def detect_has_ann(vcf_lines):
    for ln in vcf_lines:
        if isinstance(ln, bytes):
            ln = ln.decode("utf-8","ignore")
        if ln.startswith("##INFO=<ID=ANN") or "ANN=" in ln:
            return True
    return False

def run_vep_annotate(input_path, output_path, build):
    cmd = ["/bin/bash", "/app/annotate_vep.sh", input_path, output_path, build]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return res

def parse_vcf_to_df(vcf_text_lines):
    header_cols = []
    records = []
    for line in vcf_text_lines:
        if isinstance(line, bytes):
            line = line.decode("utf-8","ignore")
        line = line.rstrip("\n")
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_cols = line.lstrip("#").split("\t"); 
            continue
        parts = line.split("\t")
        if len(parts) < 8: 
            continue
        row = dict(zip(header_cols[:8], parts[:8]))
        info = parse_info(parts[7])
        row.update({f"INFO_{k}": v for k,v in info.items()})
        # ANN/CSQ — pick first transcript if exists
        ann = info.get("ANN")
        if ann:
            rec = ann.split(",")[0].split("|")
            def g(i): return rec[i] if i < len(rec) else ""
            row.update({
                "Gene": g(3), "Consequence": g(1), "Impact": g(2),
                "HGVSc": g(9), "HGVSp": g(10)
            })
        # ClinSig + AF heuristics
        row["CLNSIG"] = info.get("CLNSIG","")
        row["AF"] = None
        for k in ["gnomad_af","gnomad4_af","AF","popmax_af","gnomad_exomes_af","gnomad_genomes_af"] +                    [k for k in info.keys() if k.lower().endswith("af") or k=="AF"]:
            val = info.get(k)
            if val:
                try:
                    row["AF"] = float(str(val).split(",")[0]); break
                except: pass
        # FORMAT/SAMPLE for zygosity
        fmt = parts[8] if len(parts)>8 else None
        smp = parts[9] if len(parts)>9 else None
        if fmt and smp:
            fks = fmt.split(":"); fvs = smp.split(":")
            m = {k:v for k,v in zip(fks,fvs)}
            gt = m.get("GT")
            row["GT"] = gt
        records.append(row)
    return pd.DataFrame(records)

def clnsig_rank(x):
    x = (x or "").lower()
    if "pathogenic" in x and "benign" not in x: return 0
    if "likely_pathogenic" in x or "likely pathogenic" in x: return 1
    if "vus" in x or "uncertain" in x: return 2
    if "likely_benign" in x or "benign" in x: return 4
    return 3

def zyg_from_gt(gt):
    if not gt: return None
    gt = gt.replace("|","/")
    if gt in ("0/1","1/0","0/2","2/0","1/2","2/1"): return "HET"
    if gt in ("1/1","2/2","3/3"): return "HOM"
    return None

def parse_panel(file) -> set:
    if not file: return set()
    name = file.name.lower()
    try:
        if name.endswith(".csv"):
            df = pd.read_csv(file)
            cols = [c for c in df.columns if c.lower() in ("gene","symbol","hgnc","hgnc_symbol")]
            if cols:
                return set(df[cols[0]].astype(str).str.strip().str.upper())
            return set(df.iloc[:,0].astype(str).str.strip().str.upper())
        else:
            genes = []
            for line in file.getvalue().decode("utf-8","ignore").splitlines():
                g = line.strip()
                if g: genes.append(g)
            return set([g.upper() for g in genes])
    except Exception as e:
        st.warning(f"Could not parse panel file: {e}")
        return set()

def prioritize(df, impacts_keep, af_max, zyg_keep, panel_genes):
    if df.empty: return df
    df["ImpactLevel"] = df.get("Impact","MODIFIER").map(lambda x: IMPACT_ORDER.get(x, 3))
    df["CLNSIG_rank"] = df["CLNSIG"].map(clnsig_rank)
    df["AF_pass"] = df["AF"].map(lambda x: True if pd.isna(x) else (x <= af_max))
    df["ZYG"] = df["GT"].map(zyg_from_gt)
    if impacts_keep:
        df = df[df["Impact"].isin(impacts_keep)]
    if zyg_keep:
        df = df[df["ZYG"].isin(zyg_keep) | df["ZYG"].isna()]
    if panel_genes:
        df["BOOST"] = df["Gene"].astype(str).str.upper().isin(panel_genes).astype(int)
    else:
        df["BOOST"] = 0
    # Rank: BOOST desc, AF_pass desc, ImpactLevel asc, ClinSig asc, QUAL desc
    df = df.sort_values(["BOOST","AF_pass","ImpactLevel","CLNSIG_rank","QUAL"],
                        ascending=[False, False, True, True, False], na_position="last")
    return df

def make_html_report(sample_label, hpo_str, df_all, shortlist):
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    cols = [c for c in ["CHROM","POS","REF","ALT","QUAL","FILTER","Gene","Consequence","Impact","HGVSc","HGVSp","CLNSIG","AF","GT"] if c in shortlist.columns]
    html = f"""
    <html><head><meta charset='utf-8'><title>Variant Report</title>
    <style>
      body {{ font-family: Arial, Helvetica, sans-serif; margin: 24px; }}
      h1,h2 {{ margin-bottom: 0; }}
      table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
      th, td {{ border: 1px solid #ccc; padding: 6px 8px; }}
      th {{ background: #f6f6f6; text-align: left; }}
      .meta {{ color: #555; }}
    </style></head><body>
      <h1>Genmatix/GATomics — Variant Report</h1>
      <p class="meta">Generated: {now}</p>
      <h2>Sample</h2><p><b>{sample_label}</b></p>
      <h2>Phenotype (HPO)</h2><p>{hpo_str or '—'}</p>
      <h2>Summary</h2><p>Total parsed variants: <b>{len(df_all)}</b> · Shortlist: <b>{len(shortlist)}</b></p>
      <h2>Top Candidates (up to 100)</h2>
      {shortlist[cols].to_html(index=False)}
      <p class="meta">Notes: Self-annotated with Ensembl VEP ({'GRCh38' if True else 'GRCh37'}). Ranking favors panel genes, AF pass, HIGH/MOD impact, pathogenic ClinSig, high QUAL.</p>
    </body></html>
    """
    return html

def html_to_pdf_bytes(html_str):
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.pdfgen import canvas
        from reportlab.lib.units import cm
        from reportlab.lib.utils import simpleSplit
        import io, re
        buffer = io.BytesIO()
        c = canvas.Canvas(buffer, pagesize=A4)
        width, height = A4
        x_margin, y_margin = 1.8*cm, 1.8*cm
        textobject = c.beginText(x_margin, height - y_margin)
        textobject.setFont("Helvetica", 10)
        text = re.sub("<[^<]+?>", "", html_str)
        lines = text.splitlines()
        for line in lines:
            wrapped = simpleSplit(line, "Helvetica", 10, width - 2*x_margin)
            for w in wrapped:
                textobject.textLine(w)
                if textobject.getY() < y_margin:
                    c.drawText(textobject); c.showPage()
                    textobject = c.beginText(x_margin, height - y_margin)
                    textobject.setFont("Helvetica", 10)
        c.drawText(textobject); c.save()
        pdf_bytes = buffer.getvalue(); buffer.close()
        return pdf_bytes
    except Exception:
        return None

# ------------------ Run ------------------
if run and vcf_file is not None:
    # Save upload to temp
    suffix = ".vcf.gz" if vcf_file.name.endswith(".gz") else ".vcf"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp_in:
        tmp_in.write(vcf_file.read())
        tmp_in_path = tmp_in.name

    # Peek to see if ANN exists
    if tmp_in_path.endswith(".gz"):
        with gzip.open(tmp_in_path, "rt", errors="ignore") as f:
            lines = f.readlines()[:1000]
    else:
        with open(tmp_in_path, "rt", errors="ignore") as f:
            lines = f.readlines()[:1000]

    has_ann = detect_has_ann(lines)
    annotated_path = tmp_in_path

    if not has_ann:
        st.info("No ANN detected — running VEP annotation (offline cache)…")
        out_path = tempfile.NamedTemporaryFile(delete=False, suffix=".annot.vcf").name
        res = run_vep_annotate(tmp_in_path, out_path, genome_build)
        if res.returncode != 0:
            st.error("VEP annotation failed. Check Docker image/paths.\n\nSTDOUT:\n" + res.stdout + "\n\nSTDERR:\n" + res.stderr)
            st.stop()
        annotated_path = out_path

    # Read annotated VCF text
    with open(annotated_path, "rt", errors="ignore") as f:
        vcf_lines = f.readlines()

    # Parse & prioritize
    df = parse_vcf_to_df(vcf_lines)
    if df.empty:
        st.error("Parsed 0 variants after annotation. Verify input VCF and genome build.")
        st.stop()

    panel_genes = parse_panel(panel_file)
    df_rank = prioritize(df, impact_keep, af_thresh, zyg_keep, panel_genes)
    shortlist = df_rank[df_rank["AF_pass"]].head(100) if "AF_pass" in df_rank.columns else df_rank.head(100)

    st.subheader("Shortlist (top ≤ 100)")
    cols = [c for c in ["CHROM","POS","REF","ALT","QUAL","FILTER","Gene","Consequence","Impact","HGVSc","HGVSp","CLNSIG","AF","GT"] if c in df_rank.columns]
    st.dataframe(shortlist[cols], use_container_width=True)

    # Downloads
    st.download_button("⬇️ Download shortlist CSV", data=shortlist[cols].to_csv(index=False), file_name="shortlist.csv", mime="text/csv")

    # Reports
    sample_label = sample_name or vcf_file.name
    hpo_str = (hpo_input or "").strip().replace("\n"," ").replace(","," ")
    html = make_html_report(sample_label, hpo_str, df, shortlist)
    st.download_button("⬇️ Download HTML report", data=html, file_name="vcf_report.html", mime="text/html")

    pdf_bytes = html_to_pdf_bytes(html)
    if pdf_bytes:
        st.download_button("⬇️ Download PDF report", data=pdf_bytes, file_name="vcf_report.pdf", mime="application/pdf")
    else:
        st.info("PDF generation requires reportlab. (It’s installed in the Docker image.)")

    st.success("Done ✅")
else:
    st.caption("Select genome build, upload VCF, and click “Annotate & Generate Report”.")
