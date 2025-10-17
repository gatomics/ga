
import streamlit as st
import pandas as pd
from io import StringIO
from datetime import datetime
import base64

st.set_page_config(page_title="Genmatix/GATomics — VCF → Report", layout="wide")
st.title("🧬 Genmatix/GATomics — VCF → Variant Report")
st.caption("MVP+ | Single-sample VCF with ANN annotations (SnpEff/VEP). Optional SpliceAI/REVEL/CADD parsed if present.")

# ------------------------------ Inputs ------------------------------
col1, col2 = st.columns([2,1])
with col1:
    vcf_file = st.file_uploader("Upload VCF (.vcf or .vcf.gz)", type=["vcf","gz"])
    sample_name = st.text_input("Sample name (optional)", "")
    hpo_input = st.text_area("HPO terms (optional, comma- or space-separated)", placeholder="e.g., HP:0003774 HP:0000083")
with col2:
    panel_file = st.file_uploader("Optional: Gene panel for boosting (txt/csv, one symbol per line or 'gene' column)", type=["txt","csv"])
    zygosity = st.multiselect("Zygosity filter", ["HET","HOM"], default=["HET","HOM"])
    impact_keep = st.multiselect("Impacts to keep", ["HIGH","MODERATE","LOW","MODIFIER"], default=["HIGH","MODERATE"])
    af_thresh = st.number_input("Max AF (keep ≤)", min_value=0.0, max_value=1.0, value=0.01, step=0.001)

run = st.button("⚡ Generate Report")

IMPACT_ORDER = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

def parse_info(info_str):
    out = {}
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = True
    return out

def parse_ann_field(ann_val):
    recs = []
    for rec in ann_val.split(","):
        parts = rec.split("|")
        def g(i):
            return parts[i] if i < len(parts) else ""
        recs.append({
            "ANN_Allele": g(0),
            "Consequence": g(1),
            "Impact": g(2),
            "Gene": g(3),
            "Gene_ID": g(4),
            "Feature_Type": g(5),
            "Feature": g(6),
            "Transcript_BIOTYPE": g(7),
            "Rank": g(8),
            "HGVSc": g(9),
            "HGVSp": g(10),
        })
    return recs

def best_ann(anns):
    if not anns:
        return {}
    anns_sorted = sorted(anns, key=lambda a: IMPACT_ORDER.get(a.get("Impact","MODIFIER"), 3))
    return anns_sorted[0]

def extract_float(val):
    if val is None:
        return None
    try:
        return float(str(val).split(",")[0])
    except:
        return None

def parse_vcf_lines(lines):
    header_cols = []
    records = []
    sv_records = []
    for line in lines:
        if isinstance(line, bytes):
            line = line.decode("utf-8", errors="ignore")
        line = line.rstrip("\n")
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header_cols = line.lstrip("#").split("\t")
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        base = dict(zip(header_cols[:8], parts[:8]))
        if len(parts) > 8:
            base["FORMAT"] = parts[8]
        if len(parts) > 9:
            base["SAMPLE"] = parts[9]
        info = parse_info(base.get("INFO",""))
        base.update({f"INFO_{k}": v for k,v in info.items()})
        svtype = info.get("SVTYPE")
        if svtype:
            ann_val = info.get("ANN")
            gene = ""
            cons = ""
            if ann_val:
                ba = best_ann(parse_ann_field(ann_val))
                gene = ba.get("Gene","")
                cons = ba.get("Consequence","")
            rec = {
                "CHROM": base.get("CHROM"),
                "POS": base.get("POS"),
                "END": info.get("END"),
                "SVTYPE": svtype,
                "Gene": gene,
                "Consequence": cons,
                "QUAL": base.get("QUAL"),
                "FILTER": base.get("FILTER"),
            }
            sv_records.append(rec)
            continue

        ann_val = info.get("ANN")
        if ann_val:
            ba = best_ann(parse_ann_field(ann_val))
            base.update(ba)

        clnsig = info.get("CLNSIG","") or info.get("CLNSIGCONF","")
        base["CLNSIG"] = clnsig

        af = None
        for k in ["gnomad_af","gnomad4_af","AF","popmax_af","gnomad_exomes_af","gnomad_genomes_af"] +                 [k for k in info.keys() if k.lower().endswith("af") or k.upper()=="AF"]:
            v = info.get(k)
            af = extract_float(v)
            if af is not None:
                break
        base["AF"] = af

        spliceai = info.get("SpliceAI") or info.get("SPLICEAI")
        if spliceai:
            try:
                fields = spliceai.split("|")
                ds = []
                for i in range(2, 6):
                    if i < len(fields):
                        ds.append(extract_float(fields[i]))
                ds = [x for x in ds if x is not None]
                base["SpliceAI_DSmax"] = max(ds) if ds else None
            except:
                base["SpliceAI_DSmax"] = None
        else:
            base["SpliceAI_DSmax"] = None

        revel = info.get("REVEL") or info.get("REVEL_SCORE")
        base["REVEL"] = extract_float(revel)

        cadd = info.get("CADD") or info.get("CADD_PHRED") or info.get("CADD13_PHRED")
        base["CADD_PHRED"] = extract_float(cadd)

        if base.get("FORMAT") and base.get("SAMPLE"):
            fmt = base["FORMAT"].split(":")
            smp = base["SAMPLE"].split(":")
            fm = {k:v for k,v in zip(fmt, smp)}
            base["GT"] = fm.get("GT")
            for key in ["DP","GQ","AD","ADF","ADR","VAF"]:
                if key in fm:
                    base[f"SAMP_{key}"] = fm[key]

        records.append(base)
    return records, sv_records

def clnsig_rank(x):
    x = (x or "").lower()
    if "pathogenic" in x and "benign" not in x:
        return 0
    if "likely_pathogenic" in x or "likely pathogenic" in x:
        return 1
    if "vus" in x or "uncertain" in x:
        return 2
    if "likely_benign" in x or "benign" in x:
        return 4
    return 3

def zygosity_from_gt(gt):
    if not gt:
        return None
    gt = gt.replace("|","/")
    if gt in ("0/1","1/0","0/2","2/0","1/2","2/1"):
        return "HET"
    if gt in ("1/1","2/2","3/3"):
        return "HOM"
    return None

def prioritize(df, impacts_keep, af_max, panel_genes, zyg_keep):
    if df.empty:
        return df
    df["ImpactLevel"] = df.get("Impact","MODIFIER").map(lambda x: IMPACT_ORDER.get(x, 3))
    df["CLNSIG_rank"] = df["CLNSIG"].map(clnsig_rank)
    df["AF_pass"] = df["AF"].map(lambda x: True if x is None else (x <= af_max))
    df["ZYG"] = df["GT"].map(zygosity_from_gt)
    if zyg_keep:
        df = df[df["ZYG"].isin(zyg_keep) | df["ZYG"].isna()]
    df = df[df["Impact"].isin(impacts_keep)]
    if panel_genes:
        df["BOOST"] = df["Gene"].astype(str).str.upper().isin(panel_genes).astype(int)
    else:
        df["BOOST"] = 0
    sort_cols = ["BOOST","AF_pass","ImpactLevel","CLNSIG_rank","SpliceAI_DSmax","REVEL","CADD_PHRED","QUAL"]
    ascending = [False, False, True, True, False, False, False, False]
    df = df.sort_values(by=sort_cols, ascending=ascending, na_position="last")
    return df

def parse_panel(file) -> set:
    if file is None:
        return set()
    name = file.name.lower()
    try:
        if name.endswith(".csv"):
            df = pd.read_csv(file)
            cols = [c for c in df.columns if c.lower() in ("gene","symbol","hgnc","hgnc_symbol")]
            if cols:
                return set(df[cols[0]].astype(str).str.strip().str.upper().tolist())
            return set(df.iloc[:,0].astype(str).str.strip().str.upper().tolist())
        else:
            genes = []
            for line in file.getvalue().decode("utf-8", errors="ignore").splitlines():
                g = line.strip()
                if g:
                    genes.append(g)
            return set([g.upper() for g in genes])
    except Exception as e:
        st.warning(f"Could not parse panel file: {e}")
        return set()

def make_html_report(sample_label, hpo_str, df_all, shortlist, comphets, sv_df):
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    sel_cols = ["CHROM","POS","REF","ALT","QUAL","FILTER","Gene","Consequence","Impact","HGVSc","HGVSp",
                "CLNSIG","AF","GT","SpliceAI_DSmax","REVEL","CADD_PHRED"]
    sel_cols = [c for c in sel_cols if c in shortlist.columns]
    topn = shortlist.head(20)[sel_cols]
    comp_table = None
    if comphets is not None and not comphets.empty:
        comp_table = comphets.copy()
    sv_html = ""
    if sv_df is not None and not sv_df.empty:
        sv_html = sv_df.to_html(index=False)
    html = f"""
    <html><head><meta charset='utf-8'><title>Variant Report</title>
    <style>
      body {{ font-family: Arial, Helvetica, sans-serif; margin: 24px; }}
      h1,h2 {{ margin-bottom: 0; }}
      table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
      th, td {{ border: 1px solid #ccc; padding: 6px 8px; }}
      th {{ background: #f6f6f6; text-align: left; }}
      .meta {{ color: #555; }}
      .tag {{ background:#eef; border:1px solid #dde; padding:2px 6px; border-radius:10px; }}
    </style>
    </head><body>
      <h1>Genmatix/GATomics — Variant Report</h1>
      <p class="meta">Generated: {now}</p>
      <h2>Sample</h2>
      <p><b>{sample_label}</b></p>
      <h2>Phenotype (HPO)</h2>
      <p>{hpo_str or '—'}</p>
      <h2>Summary</h2>
      <p>Total parsed variants: <b>{len(df_all)}</b><br>
         Shortlist shown: <b>{len(shortlist)}</b></p>
      <h2>Top Candidates (up to 20)</h2>
      {topn.to_html(index=False)}
      <h2>Compound-het candidates (≥2 het variants in same gene)</h2>
      {comp_table.to_html(index=False) if comp_table is not None else '<p>—</p>'}
      <h2>Structural variants (SV)</h2>
      {sv_html or '<p>—</p>'}
      <p class="meta">Notes: Ranking = Panel boost → AF pass → Impact (HIGH/MOD) → ClinSig → SpliceAI → REVEL → CADD → QUAL. VCF must include ANN.</p>
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
        x_margin, y_margin = 2*cm, 2*cm
        textobject = c.beginText(x_margin, height - y_margin)
        textobject.setFont("Helvetica", 10)
        text = re.sub("<[^<]+?>", "", html_str)
        lines = text.splitlines()
        for line in lines:
            wrapped = simpleSplit(line, "Helvetica", 10, width - 2*x_margin)
            for w in wrapped:
                textobject.textLine(w)
                if textobject.getY() < y_margin:
                    c.drawText(textobject)
                    c.showPage()
                    textobject = c.beginText(x_margin, height - y_margin)
                    textobject.setFont("Helvetica", 10)
        c.drawText(textobject)
        c.save()
        pdf_bytes = buffer.getvalue()
        buffer.close()
        return pdf_bytes
    except Exception:
        return None

panel_genes = parse_panel(panel_file)
if run and vcf_file is not None:
    st.info("Parsing VCF…")
    content = vcf_file.read()
    try:
        import gzip
        if vcf_file.name.endswith(".gz"):
            content = gzip.decompress(content)
    except Exception:
        pass
    lines = content.decode("utf-8", errors="ignore").splitlines()
    recs, sv_recs = parse_vcf_lines(lines)
    if not recs and not sv_recs:
        st.error("No variants parsed. Ensure the VCF is valid and not multi-sample.")
        st.stop()
    df = pd.DataFrame(recs) if recs else pd.DataFrame(columns=["CHROM","POS"])
    sv_df = pd.DataFrame(sv_recs) if sv_recs else pd.DataFrame(columns=["CHROM","POS","END","SVTYPE"])

    df_filt = df.copy()
    df_filt = prioritize(df_filt, impact_keep, af_thresh, panel_genes, zygosity)

    shortlist = df_filt[df_filt["AF_pass"]] if "AF_pass" in df_filt.columns else df_filt
    shortlist = shortlist.head(100)

    st.subheader("Shortlist (top ≤ 100 variants)")
    show_cols = ["CHROM","POS","REF","ALT","QUAL","FILTER","Gene","Consequence","Impact","HGVSc","HGVSp",
                 "CLNSIG","AF","GT","SpliceAI_DSmax","REVEL","CADD_PHRED"]
    show_cols = [c for c in show_cols if c in shortlist.columns]
    st.dataframe(shortlist[show_cols], use_container_width=True)

    comp = None
    if not df.empty and "Gene" in df.columns:
        tmp = df.copy()
        tmp["ZYG"] = tmp["GT"].map(zygosity_from_gt)
        comp = (tmp[tmp["ZYG"]=="HET"]
                .groupby("Gene")
                .size().reset_index(name="Het_variant_count"))
        comp = comp[comp["Het_variant_count"]>=2].sort_values("Het_variant_count", ascending=False)
        if not comp.empty:
            st.subheader("Compound-het candidates")
            st.dataframe(comp, use_container_width=True)

    if not sv_df.empty:
        st.subheader("Structural variants (SV)")
        st.dataframe(sv_df, use_container_width=True)

    csv = shortlist[show_cols].to_csv(index=False) if not shortlist.empty else "".encode()
    st.download_button("⬇️ Download Shortlist CSV", data=csv, file_name="vcf_shortlist.csv", mime="text/csv")

    sample_label = sample_name or vcf_file.name
    hpo_str = hpo_input.strip().replace("\n"," ").replace(","," ")
    html = make_html_report(sample_label, hpo_str, df, shortlist, comp, sv_df)
    st.download_button("⬇️ Download HTML report", data=html, file_name="vcf_report.html", mime="text/html")

    pdf_bytes = html_to_pdf_bytes(html)
    if pdf_bytes:
        st.download_button("⬇️ Download PDF report", data=pdf_bytes, file_name="vcf_report.pdf", mime="application/pdf")
    else:
        st.info("PDF generation requires reportlab. Use HTML report if PDF is unavailable.")

    st.success("Done ✅")
else:
    st.caption("Tip: Include ANN annotations and (optionally) SpliceAI/REVEL/CADD in the VCF INFO for richer ranking.")
