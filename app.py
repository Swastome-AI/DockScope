"""
DockScope - CB-Dock2 Blind-Docking Pocket Analytics Suite
=========================================================
Upload CB-Dock2 `curpocket.txt` files (one by one or in bulk with a metadata
sheet), validate them, and get binding-affinity, cavity, residue and
prediction analyses with interactive Plotly figures, CSV tables and a
one-click ZIP report.

Run locally:   streamlit run app.py
"""
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from itertools import combinations

from datetime import datetime

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import stats

# ----------------------------------------------------------------------------
# Page + constants
# ----------------------------------------------------------------------------
st.set_page_config(page_title="DockScope | CB-Dock2 Analytics", page_icon="🧬", layout="wide")

APP_NAME = "DockScope"
APP_TAGLINE = "CB-Dock2 Blind-Docking Pocket Analytics Suite"

REQUIRED = ["ID", "Cavities_volume", "center_x", "center_y", "center_z",
            "size_x", "size_y", "size_z", "score"]
MANUAL = "✍️ Manual upload (one file at a time)"
BULK = "📦 Bulk upload (many files + metadata sheet)"
BEST_LOWEST = "Lowest Vina score (recommended)"
BEST_ID1 = "Pocket ID 1 (CB-Dock2 top-ranked cavity)"
R_KCAL = 0.0019872041  # kcal / (mol K)

BLACK = "black"
FOOTER_LINES = [
    "Dockscope is actively maintained and developed by Swastome-AI under Swastome Private Limited.",
    "© 2026 Swastome Private Limited. All rights reserved. Dockscope and the Swastome logo are the property of Swastome Private Limited. "
    "No part of this app may be copied, modified or distributed without written permission. "
    "If you are using Dockscop please acknowledge and refer ",
]

RES_TYPES = {
    "Hydrophobic": {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "GLY"},
    "Aromatic": {"PHE", "TYR", "TRP"},
    "Polar": {"SER", "THR", "ASN", "GLN", "CYS"},
    "Positive (+)": {"ARG", "LYS", "HIS"},
    "Negative (−)": {"ASP", "GLU"},
}
CLASS_ORDER = ["Very strong (≤ −9)", "Strong (−9 to −7)", "Moderate (−7 to −5)", "Weak (> −5)"]

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.4rem;}
    .ds-title {font-size: 2.3rem; font-weight: 800; margin-bottom: 0;}
    .ds-sub {font-size: 1.05rem; color: #5b6770; margin-top: 0;}
    .ds-footer {font-size: 0.78rem; line-height: 1.45; color: #6b7780; text-align: center; padding: 0.4rem 1rem 1.2rem;}
    div[data-testid="stMetric"] {background: rgba(120,140,160,0.08); border-radius: 10px; padding: 8px 12px;}
    </style>
    """,
    unsafe_allow_html=True,
)


# ----------------------------------------------------------------------------
# Session state
# ----------------------------------------------------------------------------
def init_state():
    defaults = dict(
        ligands=[], receptors=[],
        manual_records=[], manual_done=False,
        bulk_records=[], bulk_submitted=False, bulk_msgs=[], bulk_sig=None,
        ack_mixed=False, results=None, zip_bytes=None, zip_info="",
        flash=None, last_n=5, pdf_bytes=None, pdf_for=None,
    )
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


init_state()


# ----------------------------------------------------------------------------
# Parsing + validation
# ----------------------------------------------------------------------------
def parse_curpocket(raw: bytes):
    """Read a CB-Dock2 curpocket.txt into a DataFrame. Returns (df|None, errors)."""
    try:
        text = raw.decode("utf-8", errors="ignore")
        df = pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python")
    except Exception as e:  # noqa: BLE001
        return None, [f"File cannot be read as a table ({e})"]
    df.columns = [str(c).strip() for c in df.columns]
    errs = []
    missing = [c for c in REQUIRED if c not in df.columns]
    contact = next((c for c in df.columns if c.lower().startswith("contact")), None)
    if missing:
        errs.append("Missing column(s): " + ", ".join(missing))
    if contact is None:
        errs.append("Missing contact-residue column (contact_resi...)")
    if errs:
        return None, errs
    df = df.rename(columns={contact: "contact_residues"})
    for c in REQUIRED:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    if df[REQUIRED].isna().any().any():
        return None, ["Empty or non-numeric values found in numeric columns"]
    df["contact_residues"] = df["contact_residues"].astype(str)
    return df, []


def check_ids(df, declared):
    ids = df["ID"].astype(int).tolist()
    found = len(ids)
    problems = []
    if found != declared:
        problems.append(f"Pocket ID count mismatch: file has {found} pockets (IDs) but you declared {declared} cavities")
    if len(set(ids)) != found:
        problems.append("Duplicate pocket IDs inside the file")
    elif sorted(ids) != list(range(1, found + 1)):
        problems.append("Pocket IDs are not a continuous 1…N series")
    return problems


def make_record(fname, raw, receptor, ligand, declared, extra_issues=None):
    df, errs = parse_curpocket(raw)
    issues = list(extra_issues or []) + list(errs)
    found = None
    if df is not None:
        found = len(df)
        issues += check_ids(df, declared)
    return dict(file=fname, receptor=receptor, ligand=ligand, declared=int(declared), found=found,
                status="ok" if not issues else "flag", issues=issues, df=df if not issues else None)


def records_frame(recs):
    return pd.DataFrame([{
        "Status": "✅" if r["status"] == "ok" else "❌",
        "File": r["file"], "Receptor": r["receptor"], "Ligand": r["ligand"],
        "Declared cavities": r["declared"],
        "IDs found": str(r["found"]) if r["found"] is not None else "—",
        "Issues": "; ".join(r["issues"]) if r["issues"] else "None",
    } for r in recs])


# ----------------------------------------------------------------------------
# Bulk metadata
# ----------------------------------------------------------------------------
def read_meta(f):
    raw = f.getvalue()
    name = f.name.lower()
    try:
        if name.endswith((".xlsx", ".xls")):
            m = pd.read_excel(io.BytesIO(raw), dtype=str)
        else:
            m = pd.read_csv(io.BytesIO(raw), sep=None, engine="python", dtype=str)
    except Exception as e:  # noqa: BLE001
        return None, [f"Metadata file cannot be read ({e})"]
    mapping = {}
    for c in m.columns:
        n = re.sub(r"[^a-z]", "", str(c).lower())
        if "ligand" in n:
            mapping["Ligand"] = c
        elif "receptor" in n or "protein" in n:
            mapping["Receptor"] = c
        elif "cavit" in n or "pocket" in n:
            mapping["Cavities"] = c
        elif "file" in n:
            mapping["File_name"] = c
    need = ["Ligand", "Receptor", "Cavities", "File_name"]
    miss = [k for k in need if k not in mapping]
    if miss:
        return None, [f"Metadata is missing column(s): {', '.join(miss)}. Required: Ligand, Receptor, Number_of_cavities, File_name"]
    m = m[[mapping[k] for k in need]].copy()
    m.columns = need
    for c in need:
        m[c] = m[c].astype(str).str.strip().replace({"nan": np.nan, "": np.nan})
    return m, []


def validate_bulk(uploads, meta_file):
    msgs, recs = [], []
    m, errs = read_meta(meta_file)
    if errs:
        return [], errs
    up, dup_up = {}, []
    for f in uploads:
        if f.name in up:
            dup_up.append(f.name)
        up[f.name] = f.getvalue()
    if dup_up:
        msgs.append("Same file uploaded more than once: " + ", ".join(sorted(set(dup_up))))
    for col in ["Ligand", "Receptor"]:
        groups = defaultdict(set)
        for v in m[col].dropna():
            groups[re.sub(r"\s+", "", v.lower())].add(v)
        for variants in groups.values():
            if len(variants) > 1:
                msgs.append(f"Inconsistent spelling in {col} column: {sorted(variants)} – use one spelling everywhere")
    seen_files, seen_combo = set(), set()
    for idx, row in m.iterrows():
        line = idx + 2
        fn, lig, rec, cav = row["File_name"], row["Ligand"], row["Receptor"], row["Cavities"]
        issues = []
        if pd.isna(fn) or pd.isna(lig) or pd.isna(rec) or pd.isna(cav):
            recs.append(dict(file=fn if not pd.isna(fn) else f"(row {line})", receptor=rec if not pd.isna(rec) else "?",
                             ligand=lig if not pd.isna(lig) else "?", declared=0, found=None, status="flag",
                             issues=[f"Metadata row {line} has an empty field"], df=None))
            continue
        try:
            n = int(float(cav))
            if n < 1:
                raise ValueError
        except ValueError:
            recs.append(dict(file=fn, receptor=rec, ligand=lig, declared=0, found=None, status="flag",
                             issues=[f"Number of cavities '{cav}' is not a positive whole number"], df=None))
            continue
        if fn in seen_files:
            issues.append("File name listed more than once in metadata")
        if (rec, lig) in seen_combo:
            issues.append("Receptor + Ligand combination listed more than once")
        seen_files.add(fn)
        seen_combo.add((rec, lig))
        if fn not in up:
            issues.append("File not uploaded (name must match exactly, including .txt)")
            recs.append(dict(file=fn, receptor=rec, ligand=lig, declared=n, found=None, status="flag", issues=issues, df=None))
        else:
            recs.append(make_record(fn, up[fn], rec, lig, n, extra_issues=issues))
    for fn in up:
        if fn not in seen_files:
            recs.append(dict(file=fn, receptor="?", ligand="?", declared=0, found=None, status="flag",
                             issues=["Uploaded file is not listed in the metadata sheet"], df=None))
    return recs, msgs


def bulk_widget_sig():
    files = st.session_state.get("bulk_files") or []
    meta = st.session_state.get("bulk_meta")
    return (tuple(sorted((f.name, f.size) for f in files)), (meta.name, meta.size) if meta else None)


# ----------------------------------------------------------------------------
# Status / gating
# ----------------------------------------------------------------------------
def current_state():
    mode = st.session_state.get("mode", MANUAL)
    if mode == MANUAL:
        recs, finished, extra = st.session_state.manual_records, st.session_state.manual_done, []
    else:
        recs = st.session_state.bulk_records
        finished = st.session_state.bulk_submitted and st.session_state.bulk_sig == bulk_widget_sig()
        extra = st.session_state.bulk_msgs if finished else []
    bad = [r for r in recs if r["status"] != "ok"]
    ok = [r for r in recs if r["status"] == "ok"]
    declared = sorted({r["declared"] for r in recs if r["declared"]})
    mixed = len(declared) > 1
    ready = bool(finished and recs and not bad and not extra and (not mixed or st.session_state.ack_mixed))
    if ready:
        label = "✅ Validated"
    elif finished and (bad or extra):
        label = "❌ Fix issues"
    else:
        label = "⏳ Pending"
    return dict(mode=mode, recs=recs, ok=ok, bad=bad, finished=finished, extra=extra,
                declared=declared, mixed=mixed, ready=ready, label=label)


def make_sig(S):
    items = [(r["file"], r["receptor"], r["ligand"], r["declared"],
              int(pd.util.hash_pandas_object(r["df"], index=False).sum())) for r in S["ok"]]
    return json.dumps([S["mode"], items, st.session_state.get("best_mode"), st.session_state.get("temp_k")])


def render_banner(S):
    recs = S["recs"]
    st.markdown("##### 📋 File statistics")
    c = st.columns(6)
    c[0].metric("Files", len(recs))
    c[1].metric("Receptors", len({r["receptor"] for r in recs if r["receptor"] != "?"}))
    c[2].metric("Ligands", len({r["ligand"] for r in recs if r["ligand"] != "?"}))
    c[3].metric("Cavities / file", ", ".join(map(str, S["declared"])) or "—")
    c[4].metric("Valid files", f"{len(S['ok'])}/{len(recs)}")
    c[5].metric("Status", S["label"])
    st.divider()


def gate(S):
    if not S["ready"]:
        st.warning("🔒 This tab is locked. Finish **Tab 1 – Upload & Validate** with all files green (✅) first.")
        return None
    res = st.session_state.results
    if res is None or res["sig"] != make_sig(S):
        st.warning("🔒 Validation passed. Now click **🚀 Run prediction analysis** in **Tab 1** to unlock this tab.")
        return None
    return res


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------
def fmt_ki(nm):
    if nm < 1000:
        return f"{nm:.1f} nM"
    if nm < 1e6:
        return f"{nm / 1e3:.1f} µM"
    return f"{nm / 1e6:.1f} mM"


def affinity_class(s):
    if s <= -9:
        return CLASS_ORDER[0]
    if s <= -7:
        return CLASS_ORDER[1]
    if s <= -5:
        return CLASS_ORDER[2]
    return CLASS_ORDER[3]


def mm(s):
    rng = s.max() - s.min()
    return (s - s.min()) / rng if rng > 0 else pd.Series(0.5, index=s.index)


def res_type(name):
    for t, members in RES_TYPES.items():
        if name in members:
            return t
    return "Other"


def split_res(s):
    out = []
    for tok in str(s).split(","):
        tok = tok.strip()
        if not tok or tok == "nan":
            continue
        p = tok.split(":")
        out.append(dict(Residue=tok, ResName=p[0], ResNum=p[1] if len(p) > 1 else "", Chain=p[2] if len(p) > 2 else ""))
    return out


def build_pockets(recs):
    frames = []
    for r in recs:
        d = r["df"].copy().rename(columns={"ID": "Pocket_ID", "Cavities_volume": "Volume_A3", "score": "Score"})
        d["Pocket_ID"] = d["Pocket_ID"].astype(int)
        d.insert(0, "Complex", f"{r['receptor']} + {r['ligand']}")
        d.insert(0, "Ligand", r["ligand"])
        d.insert(0, "Receptor", r["receptor"])
        d["File"] = r["file"]
        d["Score_Rank"] = d["Score"].rank(method="min").astype(int)
        d["Delta_vs_Best"] = d["Score"] - d["Score"].min()
        frames.append(d)
    P = pd.concat(frames, ignore_index=True)
    sd = P["Score"].std()
    P["Score_Z"] = (P["Score"] - P["Score"].mean()) / sd if sd and sd > 0 else 0.0
    return P


def build_best(P, best_mode, temp_k):
    rt = R_KCAL * temp_k
    rows = []
    for (rec, lig), g in P.groupby(["Receptor", "Ligand"], sort=False):
        gs = g.sort_values("Score")
        b = g[g["Pocket_ID"] == g["Pocket_ID"].min()].iloc[0] if best_mode == BEST_ID1 else gs.iloc[0]
        sc = gs["Score"].values
        rows.append(dict(
            Receptor=rec, Ligand=lig, Complex=f"{rec} + {lig}", Best_Pocket_ID=int(b["Pocket_ID"]),
            Best_Score=float(b["Score"]), Lowest_Score_In_File=float(sc[0]),
            Second_Best_Score=float(sc[1]) if len(sc) > 1 else np.nan,
            Gap_To_Second=float(sc[1] - sc[0]) if len(sc) > 1 else 0.0,
            N_Pockets=len(g), Mean_Score=float(g["Score"].mean()),
            Competing_Pockets_1kcal=int(((g["Score"] - sc[0]) <= 1.0).sum() - 1),
            Pocket_Breadth=float((g["Score"] <= -6.0).mean()),
            Volume_A3=float(b["Volume_A3"]),
            Center_x=float(b["center_x"]), Center_y=float(b["center_y"]), Center_z=float(b["center_z"]),
            Affinity_Class=affinity_class(float(b["Score"])),
            Pred_Ki_nM=float(np.exp(float(b["Score"]) / rt) * 1e9),
            File=b["File"],
        ))
    best = pd.DataFrame(rows)
    best["Pred_Ki_uM"] = best["Pred_Ki_nM"] / 1e3
    best["Priority_Score"] = 100 * (0.70 * mm(-best["Best_Score"]) + 0.15 * mm(best["Gap_To_Second"].fillna(0))
                                    + 0.15 * mm(best["Pocket_Breadth"]))
    best["Priority_Rank"] = best["Priority_Score"].rank(ascending=False, method="min").astype(int)
    return best


def corr_text(r, p):
    if np.isnan(r):
        return "Not enough data (need ≥ 3 pockets with varying values)"
    strength = "Weak" if abs(r) < 0.3 else "Moderate" if abs(r) < 0.6 else "Strong"
    if abs(r) < 0.3:
        direction = "no clear volume trend"
    elif r < 0:
        direction = "larger cavities tend to give stronger (more negative) scores"
    else:
        direction = "larger cavities tend to give weaker scores"
    sig = "statistically notable (p < 0.05)" if p < 0.05 else "not statistically significant"
    return f"{strength} correlation; {direction}; {sig}"


def build_corr(P):
    def row(level, label, d):
        n = len(d)
        base = dict(Level=level, Group=label, N_pockets=n)
        if n < 3 or d["Volume_A3"].nunique() < 2 or d["Score"].nunique() < 2:
            return {**base, "Pearson_r": np.nan, "Pearson_p": np.nan, "Spearman_rho": np.nan,
                    "Spearman_p": np.nan, "Slope_kcal_per_100A3": np.nan, "Interpretation": corr_text(np.nan, np.nan)}
        pr, pp = stats.pearsonr(d["Volume_A3"], d["Score"])
        sr, sp = stats.spearmanr(d["Volume_A3"], d["Score"])
        slope = np.polyfit(d["Volume_A3"], d["Score"], 1)[0] * 100
        return {**base, "Pearson_r": pr, "Pearson_p": pp, "Spearman_rho": sr, "Spearman_p": sp,
                "Slope_kcal_per_100A3": slope, "Interpretation": corr_text(pr, pp)}
    rows = [row("Overall", "All pockets", P)]
    rows += [row("Ligand", k, g) for k, g in P.groupby("Ligand")]
    rows += [row("Receptor", k, g) for k, g in P.groupby("Receptor")]
    return pd.DataFrame(rows)


def build_pairs(best, temp_k):
    rt = R_KCAL * temp_k
    rows = []
    for rec, g in best.groupby("Receptor"):
        for (_, a), (_, b) in combinations(list(g.iterrows()), 2):
            ddg = b["Best_Score"] - a["Best_Score"]
            pref = a["Ligand"] if ddg > 0 else b["Ligand"] if ddg < 0 else "Tie"
            verdict = "Comparable" if abs(ddg) < 0.5 else "Slight preference" if abs(ddg) < 1.0 else "Clear preference"
            rows.append(dict(Receptor=rec, Ligand_A=a["Ligand"], Ligand_B=b["Ligand"], Pair=f"{a['Ligand']} vs {b['Ligand']}",
                             Score_A=a["Best_Score"], Score_B=b["Best_Score"], DDG_B_minus_A=ddg,
                             Stronger_Binder=pref, Ki_Fold_Difference=float(np.exp(abs(ddg) / rt)), Verdict=verdict))
    return pd.DataFrame(rows)


def build_residues(P, best):
    rows = []
    for _, b in best.iterrows():
        g = P[P["Complex"] == b["Complex"]]
        br = g[g["Pocket_ID"] == b["Best_Pocket_ID"]].iloc[0]
        cnt = Counter()
        for s in g["contact_residues"]:
            cnt.update({x["Residue"] for x in split_res(s)})
        for x in split_res(br["contact_residues"]):
            rows.append(dict(Receptor=b["Receptor"], Ligand=b["Ligand"], Complex=b["Complex"],
                             Pocket_ID=b["Best_Pocket_ID"], **x, Type=res_type(x["ResName"]),
                             N_Pockets_Contacting=cnt[x["Residue"]]))
    return pd.DataFrame(rows)


def build_jaccard(best, residues):
    rows = []
    for rec, g in best.groupby("Receptor"):
        for (_, a), (_, b) in combinations(list(g.iterrows()), 2):
            ra = set(residues[residues["Complex"] == a["Complex"]]["Residue"])
            rb = set(residues[residues["Complex"] == b["Complex"]]["Residue"])
            union = ra | rb
            j = len(ra & rb) / len(union) if union else np.nan
            dist = float(np.linalg.norm(np.array([a["Center_x"], a["Center_y"], a["Center_z"]]) -
                                        np.array([b["Center_x"], b["Center_y"], b["Center_z"]])))
            verdict = ("Largely overlapping site (likely competitive)" if j >= 0.5 else
                       "Partially overlapping site" if j >= 0.2 else "Distinct sites")
            rows.append(dict(Receptor=rec, Pair=f"{a['Ligand']} vs {b['Ligand']}", Jaccard_Similarity=j,
                             N_Shared_Residues=len(ra & rb), Shared_Residues=", ".join(sorted(ra & rb)),
                             Pocket_Centre_Distance_A=dist, Verdict=verdict))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
PALETTE = px.colors.qualitative.Bold
# Light-to-mid blue/green scale (dark navy removed) so BLACK text on heatmap cells is always legible
HEAT_SCALE = px.colors.sample_colorscale("YlGnBu_r", list(np.linspace(0.3, 1.0, 9)))
CLASS_COLORS = dict(zip(CLASS_ORDER, ["#4daf4a", "#a6d96a", "#fdd05e", "#f4837d"]))


def bundle_fig(fig):
    """Pre-compute the HTML once so Streamlit does not rebuild it on every click."""
    return {"fig": fig, "html": fig.to_html(include_plotlyjs="cdn").encode()}


def style(fig, title, x=None, y=None, h=520):
    """Common layout. All text is black on a white background so it stays readable in light AND dark app themes."""
    fig.update_layout(template="plotly_white", title=dict(text=f"<b>{title}</b>", x=0.02, font=dict(color=BLACK, size=18)),
                      font=dict(size=14, color=BLACK), paper_bgcolor="white", plot_bgcolor="white",
                      height=h, margin=dict(l=70, r=30, t=80, b=70),
                      legend=dict(bordercolor="#bbbbbb", borderwidth=1, bgcolor="white",
                                  font=dict(color=BLACK), title=dict(font=dict(color=BLACK))))
    fig.update_xaxes(tickfont=dict(color=BLACK), title_font=dict(color=BLACK), linecolor=BLACK)
    fig.update_yaxes(tickfont=dict(color=BLACK), title_font=dict(color=BLACK), linecolor=BLACK)
    fig.update_annotations(font=dict(color=BLACK))
    fig.update_coloraxes(colorbar=dict(tickfont=dict(color=BLACK), title=dict(font=dict(color=BLACK))))
    if x:
        fig.update_xaxes(title_text=f"<b>{x}</b>")
    if y:
        fig.update_yaxes(title_text=f"<b>{y}</b>")
    return fig


def fig_3d(P, receptor):
    d = P[P["Receptor"] == receptor]
    fig = px.scatter_3d(d, x="center_x", y="center_y", z="center_z", color="Score", size="Volume_A3",
                        symbol="Ligand", hover_data=["Complex", "Pocket_ID", "Volume_A3", "Score"],
                        color_continuous_scale=HEAT_SCALE, size_max=28)
    ax3d = dict(title=dict(font=dict(color=BLACK)), tickfont=dict(color=BLACK), color=BLACK)
    fig.update_layout(template="plotly_white", title=dict(text=f"<b>3D pocket map – {receptor}</b>", x=0.02, font=dict(color=BLACK, size=18)),
                      font=dict(size=13, color=BLACK), paper_bgcolor="white", height=620, margin=dict(l=0, r=0, t=70, b=0),
                      scene=dict(xaxis=ax3d, yaxis=ax3d, zaxis=ax3d),
                      legend=dict(font=dict(color=BLACK), bgcolor="white"))
    fig.update_coloraxes(colorbar=dict(tickfont=dict(color=BLACK), title=dict(font=dict(color=BLACK))))
    return fig


def build_figures(res):
    P, best = res["P"], res["best"]
    corr, pairs, comp, jac = res["corr"], res["pairs"], res["comp"], res["jac"]
    F = {}
    nrec = P["Receptor"].nunique()

    f = px.bar(best, x="Receptor", y="Best_Score", color="Ligand", barmode="group", text="Best_Score",
               color_discrete_sequence=PALETTE, hover_data=["Best_Pocket_ID", "Affinity_Class"])
    f.update_traces(texttemplate="%{text:.1f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
    f.update_yaxes(autorange="reversed")
    F["01_binding_affinity_bar"] = style(f, "Best Binding Affinity per Receptor–Ligand Complex", "Target protein / receptor",
                                         "Top Vina score (kcal/mol)  ↑ stronger")

    piv = best.pivot(index="Receptor", columns="Ligand", values="Best_Score")
    f = px.imshow(piv, text_auto=".1f", color_continuous_scale=HEAT_SCALE, aspect="auto",
                  labels=dict(x="Ligand", y="Receptor", color="Score"))
    f.update_traces(textfont=dict(color=BLACK, size=13))
    F["02_binding_affinity_heatmap"] = style(f, "Binding Affinity Heatmap (dark = stronger)", h=max(380, 110 * len(piv) + 160))

    f = px.scatter(P, x="Volume_A3", y="Score", color="Ligand", symbol="Receptor", color_discrete_sequence=PALETTE,
                   hover_data=["Complex", "Pocket_ID"])
    f.update_traces(marker=dict(size=13, line=dict(width=1, color="black")))
    if len(P) >= 3 and P["Volume_A3"].nunique() > 1:
        co = np.polyfit(P["Volume_A3"], P["Score"], 1)
        xs = np.linspace(P["Volume_A3"].min(), P["Volume_A3"].max(), 50)
        f.add_trace(go.Scatter(x=xs, y=co[0] * xs + co[1], mode="lines", name="Overall trend",
                               line=dict(color="black", dash="dash")))
    f.update_yaxes(autorange="reversed")
    F["03_volume_vs_score_scatter"] = style(f, "Cavity Volume vs Binding Score (all pockets)", "Cavity volume (Å³)",
                                            "Vina score (kcal/mol)  ↑ stronger", h=560)

    f = px.scatter(P, x="Volume_A3", y="Score", color="Ligand", facet_col="Receptor", facet_col_wrap=min(3, nrec),
                   color_discrete_sequence=PALETTE, hover_data=["Pocket_ID"])
    f.update_traces(marker=dict(size=11, line=dict(width=1, color="black")))
    f.update_yaxes(autorange="reversed")
    f.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    F["04_volume_vs_score_by_receptor"] = style(f, "Cavity Volume vs Score – split by Receptor", "Volume (Å³)", "Score (kcal/mol)",
                                                h=380 * int(np.ceil(nrec / 3)) + 80)

    cc = corr.dropna(subset=["Spearman_rho"])
    if len(cc):
        f = px.bar(cc, x="Group", y="Spearman_rho", color="Level", text="Spearman_rho", color_discrete_sequence=PALETTE)
        f.update_traces(texttemplate="%{text:.2f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
        f.update_yaxes(range=[-1.1, 1.1])
        F["05_volume_score_correlation"] = style(f, "Spearman Correlation: Cavity Volume vs Score", "Group", "Spearman ρ")

    f = px.line(P.sort_values("Pocket_ID"), x="Pocket_ID", y="Score", color="Ligand", facet_col="Receptor",
                facet_col_wrap=min(3, nrec), markers=True, color_discrete_sequence=PALETTE)
    f.update_yaxes(autorange="reversed")
    f.update_xaxes(dtick=1)
    f.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    F["06_pocket_score_profile"] = style(f, "Score Profile across Docking Pockets", "Pocket ID", "Score (kcal/mol)",
                                         h=380 * int(np.ceil(nrec / 3)) + 80)

    pv = P.pivot_table(index="Complex", columns="Pocket_ID", values="Score")
    f = px.imshow(pv, text_auto=".1f", color_continuous_scale=HEAT_SCALE, aspect="auto",
                  labels=dict(x="Pocket ID", y="Complex", color="Score"))
    f.update_xaxes(type="category")
    f.update_traces(textfont=dict(color=BLACK, size=12))
    F["07_pocket_heatmap"] = style(f, "Pocket-by-Pocket Score Heatmap", h=max(400, 60 * len(pv) + 160))

    f = px.bar(best, x="Receptor", y="Gap_To_Second", color="Ligand", barmode="group", text="Gap_To_Second",
               color_discrete_sequence=PALETTE)
    f.update_traces(texttemplate="%{text:.2f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
    F["08_site_dominance_gap"] = style(f, "Binding-Site Dominance (gap between best and 2nd-best pocket)", "Receptor",
                                       "Score gap (kcal/mol)")

    f = px.box(P, x="Receptor", y="Score", color="Ligand", points="all", color_discrete_sequence=PALETTE)
    f.update_yaxes(autorange="reversed")
    F["09_score_distribution_box"] = style(f, "Distribution of Pocket Scores", "Receptor", "Score (kcal/mol)  ↑ stronger")

    if len(comp):
        f = px.bar(comp, x="Complex", y="Percent", color="Type", barmode="stack", text="Percent",
                   color_discrete_sequence=px.colors.qualitative.Pastel)
        f.update_traces(texttemplate="%{text:.0f}%", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
        F["10_residue_composition"] = style(f, "Chemical Character of Best-Pocket Contact Residues", "Complex",
                                            "Share of contact residues (%)", h=560)

    if len(jac):
        f = px.bar(jac, x="Receptor", y="Jaccard_Similarity", color="Pair", barmode="group", text="Jaccard_Similarity",
                   color_discrete_sequence=PALETTE)
        f.update_traces(texttemplate="%{text:.2f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
        f.update_yaxes(range=[0, 1.1])
        F["11_ligand_site_overlap"] = style(f, "Binding-Site Overlap between Ligands (Jaccard of contact residues)", "Receptor",
                                            "Jaccard similarity (0–1)")

    f = px.bar(best, x="Receptor", y="Pred_Ki_nM", color="Ligand", barmode="group", text="Pred_Ki_nM",
               color_discrete_sequence=PALETTE, log_y=True)
    f.update_traces(texttemplate="%{text:.3g}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
    F["12_predicted_Ki"] = style(f, "Estimated Inhibition Constant (Ki) from Vina Score", "Receptor", "Estimated Ki (nM, log scale)  ↓ stronger")

    bs = best.sort_values("Priority_Score")
    f = px.bar(bs, x="Priority_Score", y="Complex", orientation="h", color="Affinity_Class", text="Priority_Score",
               category_orders={"Affinity_Class": CLASS_ORDER},
               color_discrete_map=CLASS_COLORS)
    f.update_traces(texttemplate="%{text:.0f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
    F["13_priority_ranking"] = style(f, "Complex Priority Ranking (composite score, 0–100)", "Priority score", "Complex",
                                     h=max(450, 45 * len(bs) + 160))

    if len(pairs):
        f = px.bar(pairs, x="Receptor", y="DDG_B_minus_A", color="Pair", barmode="group", text="DDG_B_minus_A",
                   color_discrete_sequence=PALETTE)
        f.update_traces(texttemplate="%{text:.2f}", textposition="outside", textfont=dict(color=BLACK), marker_line_color="black", marker_line_width=1)
        f.add_hline(y=0, line_dash="dash", line_color="black")
        F["14_ligand_differential_DDG"] = style(f, "Ligand Selectivity: ΔΔG = Score(B) − Score(A)  (positive → A binds stronger)",
                                                "Receptor", "ΔΔG (kcal/mol)")
    F = {k: bundle_fig(v) for k, v in F.items()}
    F["3d_maps"] = {rec: bundle_fig(fig_3d(P, rec)) for rec in P["Receptor"].unique()}
    return F


def run_analysis(S, best_mode, temp_k, cb):
    cb(0.05, "Step 1/7 · Assembling pocket table from validated files…")
    P = build_pockets(S["ok"])
    cb(0.20, "Step 2/7 · Selecting best pockets & converting scores to estimated Ki…")
    best = build_best(P, best_mode, temp_k)
    cb(0.35, "Step 3/7 · Volume–score correlation statistics…")
    corr = build_corr(P)
    cb(0.50, "Step 4/7 · Ligand-vs-ligand selectivity (ΔΔG)…")
    pairs = build_pairs(best, temp_k)
    cb(0.62, "Step 5/7 · Contact-residue profiling & site-overlap analysis…")
    residues = build_residues(P, best)
    comp = pd.DataFrame()
    if len(residues):
        comp = residues.groupby(["Complex", "Type"]).size().reset_index(name="N")
        comp["Percent"] = 100 * comp["N"] / comp.groupby("Complex")["N"].transform("sum")
    jac = build_jaccard(best, residues)
    cb(0.75, "Step 6/7 · Outlier detection & summary tables…")
    outliers = P[P["Score_Z"].abs() >= 2].sort_values("Score")[["Receptor", "Ligand", "Pocket_ID", "Volume_A3", "Score", "Score_Z"]]
    agg_lig = best.groupby("Ligand").agg(N_Receptors=("Receptor", "nunique"), Mean_Best_Score=("Best_Score", "mean"),
                                         Strongest=("Best_Score", "min"), Weakest=("Best_Score", "max")).reset_index()
    agg_rec = best.groupby("Receptor").agg(N_Ligands=("Ligand", "nunique"), Mean_Best_Score=("Best_Score", "mean"),
                                           Strongest=("Best_Score", "min"), Weakest=("Best_Score", "max")).reset_index()
    res = dict(P=P, best=best, corr=corr, pairs=pairs, residues=residues, comp=comp, jac=jac, outliers=outliers,
               agg_lig=agg_lig, agg_rec=agg_rec, sig=make_sig(S), best_mode=best_mode, temp_k=temp_k,
               val=records_frame(S["recs"]))
    cb(0.88, "Step 7/7 · Rendering interactive figures (this takes a moment)…")
    res["figs"] = build_figures(res)
    cb(1.0, "Analysis complete ✔")
    return res


# ----------------------------------------------------------------------------
# Narrative
# ----------------------------------------------------------------------------
def narrative(res):
    best, P, corr = res["best"], res["P"], res["corr"]
    top = best.loc[best["Best_Score"].idxmin()]
    lines = [
        f"This analysis covers **{len(best)} receptor–ligand complexes** ({best['Receptor'].nunique()} receptor(s) × "
        f"{best['Ligand'].nunique()} ligand(s)) and **{len(P)} docking pockets** in total.",
        f"The strongest predicted interaction is **{top['Complex']}** with **{top['Best_Score']:.1f} kcal/mol** "
        f"(pocket {int(top['Best_Pocket_ID'])}, estimated Ki ≈ {fmt_ki(top['Pred_Ki_nM'])}).",
    ]
    if best["Ligand"].nunique() >= 2:
        m = best.groupby("Ligand")["Best_Score"].mean().sort_values()
        lines.append(f"On average **{m.index[0]}** binds more strongly (mean {m.iloc[0]:.2f} kcal/mol) than "
                     f"{', '.join(m.index[1:])} ({', '.join(f'{v:.2f}' for v in m.iloc[1:])} kcal/mol).")
    if best["Receptor"].nunique() >= 2:
        m = best.groupby("Receptor")["Best_Score"].mean().sort_values()
        lines.append(f"**{m.index[0]}** is the most 'ligandable' receptor in this set (mean best score {m.iloc[0]:.2f} kcal/mol); "
                     f"**{m.index[-1]}** is the weakest ({m.iloc[-1]:.2f}).")
    cls = best["Affinity_Class"].value_counts()
    lines.append("Affinity classes: " + ", ".join(f"{cls.get(c, 0)} × {c}" for c in CLASS_ORDER) + ".")
    o = corr[corr["Level"] == "Overall"].iloc[0]
    if not np.isnan(o["Pearson_r"]):
        lines.append(f"Across all pockets, volume vs score gave Pearson r = {o['Pearson_r']:.2f} (p = {o['Pearson_p']:.3g}): {o['Interpretation'].lower()}.")
    dom = best.loc[best["Gap_To_Second"].idxmax()]
    lines.append(f"**{dom['Complex']}** shows the most dominant single site (best pocket is {dom['Gap_To_Second']:.1f} kcal/mol better than the next).")
    if len(res["jac"]):
        n_over = int((res["jac"]["Jaccard_Similarity"] >= 0.5).sum())
        lines.append(f"{n_over} of {len(res['jac'])} ligand pairs bind largely the same site (residue Jaccard ≥ 0.5).")
    return lines


# ----------------------------------------------------------------------------
# UI helpers
# ----------------------------------------------------------------------------
def block(title, body, kind="info"):
    box = {"info": st.info, "warn": st.warning, "ok": st.success}[kind]
    box(f"**{title}**\n\n{body}")


def show_fig(bundle, name):
    cfg = {"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 4, "filename": name}}
    st.plotly_chart(bundle["fig"], config=cfg, key=f"pc_{name}", theme=None)
    st.download_button("⬇️ Interactive HTML", bundle["html"], file_name=f"{name}.html", mime="text/html", key=f"dlh_{name}")


def show_table(df, name, title=None, decimals=3):
    if title:
        st.markdown(f"**{title}**")
    st.dataframe(df.round(decimals), hide_index=True)
    st.download_button("⬇️ Download table (.csv)", df.to_csv(index=False).encode(), file_name=f"{name}.csv",
                       mime="text/csv", key=f"dlt_{name}")


def summary_block(S):
    recs = S["recs"]
    show_df = records_frame(recs)
    st.dataframe(show_df, hide_index=True)
    if S["bad"]:
        st.error(f"❌ {len(S['bad'])} file(s) failed validation. **Remove them and re-upload corrected files.** "
                 "The next tabs stay locked until every file is green.")
        for r in S["bad"]:
            st.markdown(f"- **{r['file']}** → " + "; ".join(r["issues"]))
    for m in S["extra"]:
        st.error("❌ " + m)
    if S["mixed"]:
        st.warning(f"⚠️ Different 'Number of cavities' values were declared across files ({', '.join(map(str, S['declared']))}). "
                   "Within a single analysis all files should use the same number of cavities – mixed counts give biased, "
                   "non-comparable results.")
        st.checkbox("I understand the risk and want to continue anyway", key="ack_mixed")
    okr = S["ok"]
    if okr:
        rcp = sorted({r["receptor"] for r in okr})
        lig = sorted({r["ligand"] for r in okr})
        mat = pd.DataFrame({l: [("✅" if any(r["receptor"] == rc and r["ligand"] == l for r in okr) else "—") for rc in rcp]
                            for l in lig}, index=rcp)
        mat.index.name = "Receptor ↓ / Ligand →"
        st.markdown("**Coverage matrix** (✅ = file present)")
        st.dataframe(mat)
        if (mat == "—").any().any():
            st.warning("⚠️ Not every receptor–ligand combination is present. Ligand-vs-ligand comparisons will only use complete pairs.")
    if S["ready"]:
        st.success(f"✅ 100% successful – all {len(recs)} file(s) validated. Scroll down and run the prediction analysis.")


# ----------------------------------------------------------------------------
# Tab 1 – Upload
# ----------------------------------------------------------------------------
def add_name(kind, name):
    name = (name or "").strip()
    lst = st.session_state[kind]
    if not name:
        return
    if name.lower().replace(" ", "") in [x.lower().replace(" ", "") for x in lst]:
        st.session_state.flash = ("warning", f"'{name}' (or a spelling variant) is already registered.")
    else:
        lst.append(name)
    st.rerun()


def manual_ui():
    ss = st.session_state
    if ss.flash:
        getattr(st, ss.flash[0])(ss.flash[1])
        ss.flash = None
    st.markdown("#### Step A · Register your names once")
    st.caption("Names are typed only once and then chosen from drop-downs, so spelling variants can never creep in.")
    c1, c2 = st.columns(2)
    with c1, st.form("form_rec", clear_on_submit=True):
        n = st.text_input("Receptor name (e.g. CotA, KatE, KstD)")
        if st.form_submit_button("➕ Add receptor"):
            add_name("receptors", n)
    with c2, st.form("form_lig", clear_on_submit=True):
        n = st.text_input("Ligand name (e.g. Genistein, Daidzein)")
        if st.form_submit_button("➕ Add ligand"):
            add_name("ligands", n)
    st.write("**Receptors:** " + (" ".join(f"`{x}`" for x in ss.receptors) or "_none yet_") +
             "   |   **Ligands:** " + (" ".join(f"`{x}`" for x in ss.ligands) or "_none yet_"))
    used = {r["receptor"] for r in ss.manual_records} | {r["ligand"] for r in ss.manual_records}
    removable = [x for x in ss.receptors + ss.ligands if x not in used]
    if removable:
        with st.expander("Made a typo? Remove an unused registered name"):
            rm = st.selectbox("Name to remove", removable, key="rm_name")
            if st.button("🗑️ Remove name"):
                for k in ("receptors", "ligands"):
                    if rm in ss[k]:
                        ss[k].remove(rm)
                st.rerun()
    if not ss.receptors or not ss.ligands:
        st.info("Register at least one receptor and one ligand to unlock the file uploader.")
        return
    st.markdown("#### Step B · Submit one `curpocket.txt` at a time")
    with st.form("form_file", clear_on_submit=True):
        c1, c2, c3 = st.columns(3)
        rec = c1.selectbox("Receptor", ss.receptors)
        lig = c2.selectbox("Ligand", ss.ligands)
        n = c3.number_input("Number of cavities for docking", min_value=1, max_value=100, value=int(ss.last_n), step=1)
        up = st.file_uploader("Upload curpocket.txt", type=["txt"], accept_multiple_files=False)
        if st.form_submit_button("✅ Submit file", type="primary"):
            if up is None:
                ss.flash = ("error", "Please choose a file before submitting.")
            else:
                extra = []
                if any(r["receptor"] == rec and r["ligand"] == lig and r["status"] == "ok" for r in ss.manual_records):
                    extra.append(f"{rec} + {lig} was already uploaded – remove the old entry first")
                if any(r["file"] == up.name for r in ss.manual_records):
                    extra.append("A file with this name was already submitted")
                r = make_record(up.name, up.getvalue(), rec, lig, n, extra_issues=extra)
                ss.manual_records.append(r)
                ss.last_n = int(n)
                ss.manual_done = False
                ss.flash = ("success", f"✅ {up.name} ({rec} + {lig}) passed all checks.") if r["status"] == "ok" else \
                           ("error", f"❌ {up.name} failed: " + "; ".join(r["issues"]))
            st.rerun()
    st.markdown("#### Submitted files")
    if not ss.manual_records:
        st.caption("No files submitted yet.")
    for i, r in enumerate(ss.manual_records):
        c = st.columns([0.5, 3, 1.6, 1.6, 1.3, 1.3, 1.1])
        c[0].markdown("✅" if r["status"] == "ok" else "❌")
        c[1].write(r["file"])
        c[2].write(r["receptor"])
        c[3].write(r["ligand"])
        c[4].write(f"declared {r['declared']}")
        c[5].write(f"found {r['found'] if r['found'] is not None else '—'}")
        if c[6].button("Remove", key=f"rm_{i}_{r['file']}"):
            ss.manual_records.pop(i)
            ss.manual_done = False
            st.rerun()
        if r["status"] != "ok":
            st.error("; ".join(r["issues"]) + " → remove this file and re-upload the corrected one.")
    st.markdown("#### Step C · Finish")
    if st.button("🏁 Done – no more files", type="primary", disabled=not ss.manual_records):
        if any(r["status"] != "ok" for r in ss.manual_records):
            ss.flash = ("error", "Cannot finish: remove the red (❌) files and re-upload them first.")
        else:
            ss.manual_done = True
        st.rerun()


def bulk_ui():
    ss = st.session_state
    st.markdown("#### Step A · Prepare the metadata sheet")
    st.caption("One row per file. Columns: **Ligand, Receptor, Number_of_cavities, File_name** (file name must match the .txt exactly). "
               "CSV, TSV or Excel accepted.")
    tmpl = "Ligand,Receptor,Number_of_cavities,File_name\nDaidzein,KatE,5,Daidzein - Catalase.txt\nGenistein,KatE,5,Genistein - Catalase.txt\n"
    st.download_button("⬇️ Download metadata template (.csv)", tmpl.encode(), file_name="metadata_template.csv", mime="text/csv")
    st.markdown("#### Step B · Upload all files + metadata")
    c1, c2 = st.columns([2, 1])
    files = c1.file_uploader("All curpocket .txt files", type=["txt"], accept_multiple_files=True, key="bulk_files")
    meta = c2.file_uploader("Metadata file", type=["csv", "tsv", "txt", "xlsx", "xls"], key="bulk_meta")
    st.write(f"Files selected: **{len(files or [])}** · Metadata: **{'yes' if meta else 'no'}**")
    if st.button("🔍 Submit & check file stats", type="primary", disabled=not (files and meta)):
        with st.spinner("Validating files against metadata…"):
            recs, msgs = validate_bulk(files, meta)
        ss.bulk_records, ss.bulk_msgs, ss.bulk_submitted, ss.bulk_sig = recs, msgs, True, bulk_widget_sig()
        st.rerun()
    if ss.bulk_submitted and ss.bulk_sig != bulk_widget_sig():
        st.warning("Uploaded files or metadata changed since the last check – click **Submit & check file stats** again.")


def tab_upload(S):
    st.markdown("### 1️⃣ Upload, validate & run the analysis")
    block("What happens here", "You provide your CB-Dock2 `curpocket.txt` files. DockScope checks that every file has the expected columns, "
          "that the number of pocket IDs equals the number of cavities you declared, and that names are consistent. Only when everything is "
          "green can you run the prediction analysis that unlocks the other tabs.")
    c1, c2 = st.columns(2)
    with c1:
        st.warning("**⚠️ Warnings**\n\n"
                   "- Docking scores are *in-silico predictions*, not measured affinities.\n"
                   "- **Use the same 'Number of cavities for docking' for every file in one analysis.** One file with 5 and another with 10 "
                   "cavities gives non-comparable, biased results.\n"
                   "- File names in a bulk metadata sheet must match the uploaded `.txt` names exactly.\n"
                   "- Remove and re-upload any red (❌) file – you cannot proceed until all are green.")
    with c2:
        st.info("**📝 Notes**\n\n"
                "- CB-Dock2 numbers pockets by cavity size, so **ID 1 is the largest cavity, not automatically the best score**. "
                "Default 'best pocket' here = lowest Vina score (switch in the sidebar).\n"
                "- Estimated Ki values assume Ki ≈ exp(ΔG/RT) – useful for ranking, not for absolute numbers.\n"
                "- With only 3–10 pockets per complex, correlations are indicative only.\n"
                "- Use the camera icon in each plot toolbar for a 4× high-resolution PNG.")
    st.info("**📥 Where do I get the input file?** Run blind docking on the CB-Dock2 server "
            "(http://cao.labshare.cn:11080/cb-dock2/php/blinddock.php), download the **curpocket.txt** file from the results, "
            "and upload it here without editing. Files from other tools will not pass validation.")
    st.markdown("#### Choose how to upload")
    st.radio("Upload mode", [MANUAL, BULK], key="mode", horizontal=True, label_visibility="collapsed")
    mode = st.session_state.mode
    if mode == MANUAL:
        manual_ui()
    else:
        bulk_ui()
    S = current_state()
    st.divider()
    st.markdown("#### File-check summary")
    if not S["finished"]:
        st.info("⏳ Waiting for you to finish uploading" + (" and click **Done**." if mode == MANUAL else " and submit."))
    else:
        summary_block(S)
    if S["ready"]:
        st.markdown("#### Run prediction analysis")
        if st.button("🚀 Run prediction analysis", type="primary"):
            bar = st.progress(0, text="Starting…")
            res = run_analysis(S, st.session_state.best_mode, st.session_state.temp_k, lambda f, t: bar.progress(f, text=t))
            st.session_state.results = res
            st.session_state.zip_bytes = None
            st.session_state.pdf_bytes = None
            st.success("Analysis finished – open the next tabs (2️⃣ → 7️⃣) to explore results.")
        elif st.session_state.results is not None and st.session_state.results["sig"] == make_sig(S):
            st.success("Analysis is up to date – explore tabs 2️⃣ → 7️⃣.")


# ----------------------------------------------------------------------------
# Analysis tabs
# ----------------------------------------------------------------------------
def tab_affinity(res):
    best, F = res["best"], res["figs"]
    st.markdown("### 2️⃣ Binding affinity")
    block("What is happening", f"For each receptor–ligand file, the **best pocket** ({res['best_mode'].split(' (')[0].lower()}) is selected and its AutoDock-Vina "
          "score (kcal/mol) is reported. More negative = stronger predicted binding. Scores are also converted into an estimated Ki and "
          "an affinity class.")
    top, low = best.loc[best["Best_Score"].idxmin()], best.loc[best["Best_Score"].idxmax()]
    block("Key findings", f"- Strongest: **{top['Complex']}** at {top['Best_Score']:.1f} kcal/mol (≈ {fmt_ki(top['Pred_Ki_nM'])}).\n"
          f"- Weakest: **{low['Complex']}** at {low['Best_Score']:.1f} kcal/mol.\n- Mean best score overall: {best['Best_Score'].mean():.2f} kcal/mol.", "ok")
    show_fig(F["01_binding_affinity_bar"], "01_binding_affinity_bar")
    show_fig(F["02_binding_affinity_heatmap"], "02_binding_affinity_heatmap")
    show_table(best[["Priority_Rank", "Receptor", "Ligand", "Best_Pocket_ID", "Best_Score", "Affinity_Class", "Pred_Ki_nM", "Volume_A3", "File"]]
               .sort_values("Priority_Rank"), "best_pocket_table", "Best-pocket table")
    c1, c2 = st.columns(2)
    with c1:
        show_table(res["agg_lig"], "summary_by_ligand", "Summary by ligand")
    with c2:
        show_table(res["agg_rec"], "summary_by_receptor", "Summary by receptor")
    block("How to interpret", "- **≤ −9 kcal/mol**: very strong predicted binding; **−9 to −7**: strong; **−7 to −5**: moderate; **> −5**: weak / likely non-specific.\n"
          "- Differences < ~0.5 kcal/mol between ligands are within docking noise – treat them as *comparable*.\n"
          "- A bar that is clearly higher (more negative) for one ligand across several receptors suggests a consistently better binder.\n"
          "- Always follow up top hits with MD simulation or experiment before drawing conclusions.")


def tab_volume(res):
    corr, F = res["corr"], res["figs"]
    st.markdown("### 3️⃣ Cavity volume vs binding score")
    block("What is happening", "Every detected cavity (pocket) is a point. We test whether larger cavities systematically give better scores using "
          "Pearson (linear) and Spearman (rank-based) correlation, overall and per ligand / receptor.")
    o = corr[corr["Level"] == "Overall"].iloc[0]
    block("Key finding", o["Interpretation"] + (f" (Pearson r = {o['Pearson_r']:.2f}, n = {int(o['N_pockets'])})" if not np.isnan(o["Pearson_r"]) else ""), "ok")
    show_fig(F["03_volume_vs_score_scatter"], "03_volume_vs_score_scatter")
    show_fig(F["04_volume_vs_score_by_receptor"], "04_volume_vs_score_by_receptor")
    if "05_volume_score_correlation" in F:
        show_fig(F["05_volume_score_correlation"], "05_volume_score_correlation")
    show_table(corr, "volume_score_correlation_table", "Correlation statistics")
    show_table(res["P"][["Receptor", "Ligand", "Pocket_ID", "Volume_A3", "Score", "Score_Rank"]], "all_pockets_table", "All pockets")
    block("How to interpret", "- **Negative r**: bigger cavities → more negative (stronger) scores. Vina often rewards size because large pockets allow more contacts.\n"
          "- **r near 0**: score is driven by pocket chemistry, not size – a good sign for specific binding.\n"
          "- A pocket that is *small yet scores very well* (bottom-right of the stronger-up plot) is a particularly specific site worth inspecting.\n"
          "- Small n (≤ 10 points) → rely on the scatter plot more than on p-values.")


def tab_landscape(res):
    best, F, P = res["best"], res["figs"], res["P"]
    st.markdown("### 4️⃣ Pocket landscape")
    block("What is happening", "Within each complex we compare *all* pockets: how scores spread, how much better the best pocket is than the runner-up "
          "(**dominance gap**), and how many alternative pockets are within 1 kcal/mol (**competing sites**).")
    dom = best.loc[best["Gap_To_Second"].idxmax()]
    pro = best.loc[best["Competing_Pockets_1kcal"].idxmax()]
    block("Key findings", f"- Most dominant single site: **{dom['Complex']}** (gap {dom['Gap_To_Second']:.2f} kcal/mol).\n"
          f"- Most promiscuous: **{pro['Complex']}** with {int(pro['Competing_Pockets_1kcal'])} other pocket(s) within 1 kcal/mol of the best.\n"
          f"- Outlier pockets (|z| ≥ 2): {len(res['outliers'])}.", "ok")
    show_fig(F["06_pocket_score_profile"], "06_pocket_score_profile")
    show_fig(F["07_pocket_heatmap"], "07_pocket_heatmap")
    show_fig(F["08_site_dominance_gap"], "08_site_dominance_gap")
    show_fig(F["09_score_distribution_box"], "09_score_distribution_box")
    show_table(best[["Receptor", "Ligand", "Best_Score", "Second_Best_Score", "Gap_To_Second", "Competing_Pockets_1kcal", "Pocket_Breadth", "N_Pockets"]],
               "pocket_landscape_table", "Dominance & competition metrics")
    if len(res["outliers"]):
        show_table(res["outliers"], "outlier_pockets", "Outlier pockets (|z-score| ≥ 2)")
    block("How to interpret", "- **Large gap (> 1.5 kcal/mol)**: one clearly preferred site → ligand likely binds specifically.\n"
          "- **Small gap + many competing pockets**: ligand can sit in several places → blind docking is ambiguous; consider targeted docking per site.\n"
          "- **Pocket breadth** = fraction of pockets scoring ≤ −6 kcal/mol; high breadth means the ligand is 'sticky' across the protein surface.\n"
          "- Outliers are pockets whose score is unusually good/bad relative to the entire data set.")


def tab_residues(res):
    P, best, F, residues, jac = res["P"], res["best"], res["figs"], res["residues"], res["jac"]
    st.markdown("### 5️⃣ Binding site & residues")
    block("What is happening", "The contact residues of each complex's best pocket are profiled by chemical character (hydrophobic, aromatic, polar, charged). "
          "If you have more than one ligand per receptor, we also check whether ligands compete for the same site (Jaccard overlap of residues "
          "and pocket-centre distance). A 3-D map shows where every pocket sits on the receptor.")
    if len(jac):
        n_ov = int((jac["Jaccard_Similarity"] >= 0.5).sum())
        block("Key finding", f"{n_ov} of {len(jac)} ligand pair(s) share largely the same binding site.", "ok")
    if "10_residue_composition" in F:
        show_fig(F["10_residue_composition"], "10_residue_composition")
    if "11_ligand_site_overlap" in F:
        show_fig(F["11_ligand_site_overlap"], "11_ligand_site_overlap")
        show_table(jac, "ligand_site_overlap_table", "Ligand site-overlap table")
    st.markdown("**3-D pocket explorer**")
    rec = st.selectbox("Receptor", sorted(P["Receptor"].unique()), key="sel_3d")
    show_fig(F["3d_maps"][rec], f"3D_pocket_map_{re.sub(r'[^A-Za-z0-9_-]+', '_', rec)}")
    st.markdown("**Contact residues of the best pocket**")
    cplx = st.selectbox("Complex", sorted(residues["Complex"].unique()), key="sel_res")
    show_table(residues[residues["Complex"] == cplx].drop(columns=["Receptor", "Ligand"]), "residues_" + re.sub(r"[^A-Za-z0-9_-]+", "_", cplx))
    show_table(residues, "contact_residues_all_complexes", "All complexes")
    block("How to interpret", "- Many **aromatic (Phe/Tyr/Trp)** and hydrophobic residues → π-stacking / hydrophobic packing, typical for flavonoid-like ligands such as genistein and daidzein.\n"
          "- **Polar / charged residues** (Ser, Thr, Asn, Arg, Lys, Asp, Glu, His) suggest hydrogen bonds or salt bridges with the hydroxyl groups.\n"
          "- **Jaccard ≥ 0.5** → ligands compete for the same pocket; **< 0.2** → they likely use different sites.\n"
          "- Residues present in *many* pockets (N_Pockets_Contacting) sit where cavities overlap – they are often the structural core of the site.")


def tab_predictions(res):
    best, pairs, F = res["best"], res["pairs"], res["figs"]
    st.markdown("### 6️⃣ Predictions & ranking")
    block("What is happening", "Vina scores are converted into **estimated Ki** (Ki ≈ exp(ΔG / RT) at "
          f"{res['temp_k']:.1f} K), ligands are compared head-to-head (**ΔΔG**, fold-difference in Ki), and every complex receives a composite "
          "**priority score** = 70% binding strength + 15% site dominance + 15% pocket breadth (each min–max scaled).")
    top = best.sort_values("Priority_Rank").iloc[0]
    block("Key finding", f"Top-priority complex: **{top['Complex']}** (priority {top['Priority_Score']:.0f}/100, {top['Best_Score']:.1f} kcal/mol, ≈ {fmt_ki(top['Pred_Ki_nM'])}).", "ok")
    show_fig(F["13_priority_ranking"], "13_priority_ranking")
    show_fig(F["12_predicted_Ki"], "12_predicted_Ki")
    show_table(best[["Priority_Rank", "Complex", "Priority_Score", "Best_Score", "Pred_Ki_nM", "Pred_Ki_uM", "Affinity_Class"]]
               .sort_values("Priority_Rank"), "priority_and_Ki_table", "Ranking table")
    if len(pairs):
        show_fig(F["14_ligand_differential_DDG"], "14_ligand_differential_DDG")
        show_table(pairs, "ligand_selectivity_table", "Ligand-vs-ligand comparison")
        if len(pairs) >= 3:
            t, p = stats.ttest_1samp(pairs["DDG_B_minus_A"], 0.0)
            st.caption(f"One-sample t-test of ΔΔG vs 0 across {len(pairs)} comparisons: t = {t:.2f}, p = {p:.3g} "
                       "(indicative only; receptors are not independent in blind docking).")
    else:
        st.info("Ligand-vs-ligand analysis needs ≥ 2 ligands docked to the same receptor.")
    block("How to interpret", "- **Estimated Ki**: lower = tighter. A 1.4 kcal/mol difference ≈ 10-fold change in Ki at 298 K.\n"
          "- **Priority score** is a *screening heuristic* for deciding which complexes to validate first – it is not a probability of activity.\n"
          "- **ΔΔG > 0** means ligand A binds stronger than B; |ΔΔG| < 0.5 = comparable, 0.5–1.0 = slight, > 1.0 = clear preference.\n"
          "- Vina has ~2–3 kcal/mol typical error; rankings are more reliable than absolute values.", "warn")


def build_summary_md(res):
    lines = [f"# {APP_NAME} – analysis report", "", f"Best-pocket definition: {res['best_mode']}  ", f"Temperature for Ki: {res['temp_k']} K", "",
             "## Executive summary", ""] + [f"- {x}" for x in narrative(res)]
    lines += ["", "## Notes", "- Docking scores are predictions, not measurements.",
              "- All files should use the same number of docking cavities.",
              "- Estimated Ki and priority scores are screening heuristics."]
    return "\n".join(lines)


def build_zip(res, prog):
    buf = io.BytesIO()
    tables = {"file_validation": res["val"], "all_pockets": res["P"].drop(columns=["contact_residues"]), "best_pocket_table": res["best"],
              "volume_score_correlation": res["corr"], "ligand_selectivity": res["pairs"], "contact_residues": res["residues"],
              "ligand_site_overlap": res["jac"], "outlier_pockets": res["outliers"], "summary_by_ligand": res["agg_lig"],
              "summary_by_receptor": res["agg_rec"]}
    figs = res["figs"]
    total = len(tables) + len(figs) + 2
    done = 0
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, df in tables.items():
            if len(df):
                z.writestr(f"tables/{name}.csv", df.to_csv(index=False))
            done += 1
            prog(done / total, f"Writing table {name}…")
        for name, bundle in figs.items():
            if name == "3d_maps":
                continue
            z.writestr(f"plots_interactive_html/{name}.html", bundle["html"])
            done += 1
            prog(done / total, f"Zipping figure {name}…")
        for rec, bundle in figs.get("3d_maps", {}).items():
            z.writestr(f"plots_interactive_html/3D_pocket_map_{re.sub(r'[^A-Za-z0-9_-]+', '_', rec)}.html", bundle["html"])
        z.writestr("SUMMARY.md", build_summary_md(res))
        done += 1
        readme = ("DockScope export\n================\ntables/ - CSV tables\nplots_interactive_html/ - interactive Plotly figures (open in a browser; "
                  "use the camera icon in the toolbar to save a high-resolution PNG)\nSUMMARY.md - executive summary\n")
        z.writestr("README.txt", readme)
        done += 1
        prog(1.0, "ZIP ready ✔")
    n_plots = len(figs) - 1 + len(figs.get("3d_maps", {}))
    return buf.getvalue(), n_plots


# ----------------------------------------------------------------------------
# PDF report (ReportLab + Matplotlib – no browser / kaleido needed)
# ----------------------------------------------------------------------------
def _pdf_fonts():
    """Register DejaVu Sans (ships with Matplotlib) so Å, Δ, ≈, ≤, −, ρ … render correctly in the PDF."""
    import os
    import matplotlib
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    d = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
    files = {"DSans": "DejaVuSans.ttf", "DSans-Bold": "DejaVuSans-Bold.ttf",
             "DSans-Italic": "DejaVuSans-Oblique.ttf", "DSans-BoldItalic": "DejaVuSans-BoldOblique.ttf"}
    for name, fn in files.items():
        if name not in pdfmetrics.getRegisteredFontNames():
            pdfmetrics.registerFont(TTFont(name, os.path.join(d, fn)))
    pdfmetrics.registerFontFamily("DSans", normal="DSans", bold="DSans-Bold", italic="DSans-Italic", boldItalic="DSans-BoldItalic")


def _pdf_figures(res):
    """Static, black-on-white Matplotlib versions of the 14 analysis figures. Returns {key: (caption, png_bytes)}."""
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.ticker import MaxNLocator
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    P, best, corr, pairs, comp, jac = res["P"], res["best"], res["corr"], res["pairs"], res["comp"], res["jac"]
    cmap = LinearSegmentedColormap.from_list("ds", matplotlib.colormaps["YlGnBu_r"](np.linspace(0.3, 1.0, 256)))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    out = {}

    def colors_for(names):
        names = list(dict.fromkeys(names))
        return {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(names)}

    def polish(ax, title=None, xl=None, yl=None, grid="y"):
        ax.set_facecolor("white")
        ax.grid(axis=grid, color="#e3e3e3", lw=0.7)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        if title:
            ax.set_title(title, loc="left", fontweight="bold", fontsize=10.5)
        if xl:
            ax.set_xlabel(xl, fontweight="bold")
        if yl:
            ax.set_ylabel(yl, fontweight="bold")

    def save(key, caption, fig):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=200, bbox_inches="tight", facecolor="white")
        out[key] = (caption, buf.getvalue())

    def grouped(key, caption, df, x, y, hue, fmt, xl, yl, log=False, invert=False, ylim=None, hline=None, legend_title=None):
        fig = Figure(figsize=(7.2, 4.3))
        ax = fig.subplots()
        cats, hues = list(dict.fromkeys(df[x])), list(dict.fromkeys(df[hue]))
        col, n = colors_for(hues), len(hues)
        width = 0.8 / n
        for j, h in enumerate(hues):
            vals = df[df[hue] == h].drop_duplicates(x).set_index(x).reindex(cats)[y].values.astype(float)
            xs = np.arange(len(cats)) + (j - (n - 1) / 2) * width
            ax.bar(xs, vals, width, color=col[h], edgecolor="black", linewidth=0.6, label=h, log=log)
            for xv, v in zip(xs, vals):
                if np.isnan(v):
                    continue
                up = (v >= 0) != invert
                ax.annotate(fmt.format(v), (xv, v), xytext=(0, 3 if up else -3), textcoords="offset points",
                            ha="center", va="bottom" if up else "top", fontsize=7, color="black")
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats, rotation=0 if len(cats) <= 4 else 30, ha="center" if len(cats) <= 4 else "right")
        if invert:
            ax.invert_yaxis()
        if ylim:
            ax.set_ylim(*ylim)
        if hline is not None:
            ax.axhline(hline, color="black", ls="--", lw=0.9)
        polish(ax, None, xl, yl)
        ax.margins(y=0.12)
        ax.legend(title=legend_title or hue, loc="upper center", bbox_to_anchor=(0.5, -0.2 if len(cats) <= 4 else -0.3),
                  ncol=min(n, 3), frameon=True, edgecolor="#bbbbbb", fontsize=8, title_fontsize=8)
        save(key, caption, fig)

    def heat(key, caption, piv, xl, yl, cbar, h):
        fig = Figure(figsize=(7.2, h))
        ax = fig.subplots()
        im = ax.imshow(piv.values.astype(float), cmap=cmap, aspect="auto")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels([str(c) for c in piv.columns], rotation=0 if piv.shape[1] <= 8 else 45, ha="center" if piv.shape[1] <= 8 else "right")
        ax.set_yticks(range(piv.shape[0]))
        ax.set_yticklabels([str(i) for i in piv.index])
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if pd.notna(v):
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=8, color="black")
        ax.set_xlabel(xl, fontweight="bold")
        ax.set_ylabel(yl, fontweight="bold")
        for sp in ax.spines.values():
            sp.set_visible(False)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(cbar, color="black")
        cb.ax.tick_params(colors="black")
        save(key, caption, fig)

    nrec = P["Receptor"].nunique()
    ncol = min(3, nrec)
    nrow = int(np.ceil(nrec / ncol))
    lig_col = colors_for(P["Ligand"])

    def facets(key, caption, kind):
        fig = Figure(figsize=(7.2, 2.9 * nrow + 0.6))
        axes = fig.subplots(nrow, ncol, sharey=True, squeeze=False)
        for k, (rec, g) in enumerate(P.groupby("Receptor", sort=False)):
            ax = axes[k // ncol][k % ncol]
            for lig, gl in g.groupby("Ligand", sort=False):
                if kind == "scatter":
                    ax.scatter(gl["Volume_A3"], gl["Score"], s=34, color=lig_col[lig], edgecolor="black", lw=0.6, zorder=3)
                else:
                    gl = gl.sort_values("Pocket_ID")
                    ax.plot(gl["Pocket_ID"], gl["Score"], marker="o", ms=4, color=lig_col[lig], lw=1.4)
            polish(ax, rec, "Cavity volume (Å³)" if kind == "scatter" else "Pocket ID",
                   "Score (kcal/mol)" if k % ncol == 0 else None, grid="both")
            ax.title.set_fontsize(9)
            ax.tick_params(labelsize=7)
            if kind == "line":
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        for k in range(nrec, nrow * ncol):
            axes[k // ncol][k % ncol].set_visible(False)
        axes[0][0].invert_yaxis()
        fig.legend(handles=[Line2D([], [], marker="o", ls="", color=c, mec="black", label=l) for l, c in lig_col.items()],
                   title="Ligand", loc="lower center", ncol=min(len(lig_col), 4), frameon=True, edgecolor="#bbbbbb", fontsize=8,
                   title_fontsize=8, bbox_to_anchor=(0.5, -0.04))
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        save(key, caption, fig)

    # 01 / 02 – binding affinity
    grouped("01", "Best binding affinity per receptor–ligand complex (higher bar = stronger binding)", best, "Receptor", "Best_Score",
            "Ligand", "{:.1f}", "Target protein / receptor", "Top Vina score (kcal/mol)  ↑ stronger", invert=True)
    heat("02", "Binding-affinity heat-map (kcal/mol; darker = stronger)", best.pivot(index="Receptor", columns="Ligand", values="Best_Score"),
         "Ligand", "Receptor", "Score (kcal/mol)", max(3.0, 0.6 * best["Receptor"].nunique() + 1.8))

    # 03 – scatter (all pockets)
    fig = Figure(figsize=(7.2, 4.6))
    ax = fig.subplots()
    recs = list(dict.fromkeys(P["Receptor"]))
    for (rec, lig), g in P.groupby(["Receptor", "Ligand"], sort=False):
        ax.scatter(g["Volume_A3"], g["Score"], s=48, marker=markers[recs.index(rec) % len(markers)], color=lig_col[lig],
                   edgecolor="black", lw=0.7, zorder=3)
    if len(P) >= 3 and P["Volume_A3"].nunique() > 1:
        co = np.polyfit(P["Volume_A3"], P["Score"], 1)
        xs = np.linspace(P["Volume_A3"].min(), P["Volume_A3"].max(), 50)
        ax.plot(xs, co[0] * xs + co[1], "k--", lw=1.1)
    ax.invert_yaxis()
    polish(ax, None, "Cavity volume (Å³)", "Vina score (kcal/mol)  ↑ stronger", grid="both")
    handles = [Line2D([], [], marker="o", ls="", color=c, mec="black", label=l) for l, c in lig_col.items()]
    handles += [Line2D([], [], marker=markers[i % len(markers)], ls="", color="white", mec="black", label=r) for i, r in enumerate(recs)]
    handles += [Line2D([], [], ls="--", color="black", label="Overall trend")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=4, frameon=True, edgecolor="#bbbbbb", fontsize=8)
    save("03", "Cavity volume vs Vina score for every pocket (colour = ligand, symbol = receptor)", fig)

    facets("04", "Cavity volume vs score, split by receptor", "scatter")

    # 05 – correlation
    cc = corr.dropna(subset=["Spearman_rho"])
    if len(cc):
        fig = Figure(figsize=(7.2, 4.0))
        ax = fig.subplots()
        lv_col = colors_for(cc["Level"])
        ax.bar(range(len(cc)), cc["Spearman_rho"], color=[lv_col[l] for l in cc["Level"]], edgecolor="black", lw=0.6)
        for i, v in enumerate(cc["Spearman_rho"]):
            ax.annotate(f"{v:.2f}", (i, v), xytext=(0, 3 if v >= 0 else -3), textcoords="offset points", ha="center",
                        va="bottom" if v >= 0 else "top", fontsize=7.5)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylim(-1.1, 1.1)
        ax.set_xticks(range(len(cc)))
        ax.set_xticklabels(cc["Group"], rotation=30, ha="right")
        polish(ax, None, "Group", "Spearman ρ")
        ax.legend(handles=[Patch(facecolor=c, edgecolor="black", label=l) for l, c in lv_col.items()], title="Level", fontsize=8,
                  title_fontsize=8, frameon=True, edgecolor="#bbbbbb", loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3)
        save("05", "Spearman correlation between cavity volume and score", fig)

    facets("06", "Score profile across docking pockets", "line")

    pv = P.pivot_table(index="Complex", columns="Pocket_ID", values="Score")
    heat("07", "Pocket-by-pocket score heat-map (kcal/mol)", pv, "Pocket ID", "Complex", "Score (kcal/mol)", max(3.2, 0.42 * len(pv) + 1.8))

    grouped("08", "Binding-site dominance: gap between the best and 2nd-best pocket", best, "Receptor", "Gap_To_Second", "Ligand",
            "{:.2f}", "Receptor", "Score gap (kcal/mol)")

    # 09 – box plot
    fig = Figure(figsize=(7.2, 4.3))
    ax = fig.subplots()
    cats, hues = list(dict.fromkeys(P["Receptor"])), list(dict.fromkeys(P["Ligand"]))
    width = 0.8 / len(hues)
    for j, h in enumerate(hues):
        data, pos = [], []
        for i, c in enumerate(cats):
            v = P[(P["Receptor"] == c) & (P["Ligand"] == h)]["Score"].values
            if len(v):
                data.append(v)
                pos.append(i + (j - (len(hues) - 1) / 2) * width)
        if data:
            bp = ax.boxplot(data, positions=pos, widths=width * 0.85, patch_artist=True, showfliers=False,
                            medianprops=dict(color="black"), boxprops=dict(edgecolor="black"),
                            whiskerprops=dict(color="black"), capprops=dict(color="black"))
            for b in bp["boxes"]:
                b.set_facecolor(lig_col[h])
            for v, p_ in zip(data, pos):
                ax.scatter(np.full(len(v), p_), v, s=10, color="black", zorder=4)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats)
    ax.invert_yaxis()
    polish(ax, None, "Receptor", "Score (kcal/mol)  ↑ stronger")
    ax.legend(handles=[Patch(facecolor=c, edgecolor="black", label=l) for l, c in lig_col.items()], title="Ligand", fontsize=8,
              title_fontsize=8, frameon=True, edgecolor="#bbbbbb", loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=min(len(lig_col), 4))
    save("09", "Distribution of pocket scores per receptor and ligand", fig)

    # 10 – residue composition
    if len(comp):
        pc = comp.pivot_table(index="Complex", columns="Type", values="Percent", fill_value=0)
        tcol = {t: px.colors.qualitative.Pastel[i % 11] for i, t in enumerate(pc.columns)}
        fig = Figure(figsize=(7.2, 4.6))
        ax = fig.subplots()
        bottom = np.zeros(len(pc))
        for t in pc.columns:
            ax.bar(range(len(pc)), pc[t].values, bottom=bottom, color=tcol[t], edgecolor="black", lw=0.6, label=t)
            for i, v in enumerate(pc[t].values):
                if v >= 6:
                    ax.text(i, bottom[i] + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=7, color="black")
            bottom += pc[t].values
        ax.set_xticks(range(len(pc)))
        ax.set_xticklabels(pc.index, rotation=30, ha="right", fontsize=7.5)
        polish(ax, None, "Complex", "Share of contact residues (%)")
        ax.legend(title="Type", fontsize=8, title_fontsize=8, frameon=True, edgecolor="#bbbbbb", loc="upper center",
                  bbox_to_anchor=(0.5, -0.38), ncol=3)
        save("10", "Chemical character of best-pocket contact residues", fig)

    if len(jac):
        grouped("11", "Binding-site overlap between ligands (Jaccard similarity of contact residues)", jac, "Receptor", "Jaccard_Similarity",
                "Pair", "{:.2f}", "Receptor", "Jaccard similarity (0–1)", ylim=(0, 1.15))

    grouped("12", "Estimated Ki from Vina score (log scale; lower = tighter)", best, "Receptor", "Pred_Ki_nM", "Ligand", "{:.3g}",
            "Receptor", "Estimated Ki (nM, log scale)  ↓ stronger", log=True)

    # 13 – priority ranking
    bs = best.sort_values("Priority_Score")
    fig = Figure(figsize=(7.2, max(3.2, 0.34 * len(bs) + 1.4)))
    ax = fig.subplots()
    ax.barh(range(len(bs)), bs["Priority_Score"], color=[CLASS_COLORS[c] for c in bs["Affinity_Class"]], edgecolor="black", lw=0.6)
    for i, v in enumerate(bs["Priority_Score"]):
        ax.text(v + 1, i, f"{v:.0f}", va="center", fontsize=7.5)
    ax.set_yticks(range(len(bs)))
    ax.set_yticklabels(bs["Complex"], fontsize=7.5)
    ax.set_xlim(0, 108)
    polish(ax, None, "Priority score (0–100)", None, grid="x")
    ax.legend(handles=[Patch(facecolor=c, edgecolor="black", label=l) for l, c in CLASS_COLORS.items()], title="Affinity class",
              fontsize=7.5, title_fontsize=8, frameon=True, edgecolor="#bbbbbb", loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    save("13", "Composite priority ranking of complexes (70% affinity, 15% dominance, 15% breadth)", fig)

    if len(pairs):
        grouped("14", "Ligand selectivity: ΔΔG = Score(B) − Score(A); positive means A binds stronger", pairs, "Receptor", "DDG_B_minus_A",
                "Pair", "{:.2f}", "Receptor", "ΔΔG (kcal/mol)", hline=0)
    return out


def _md_to_rl(t):
    """Tiny markdown (**bold**, *italic*) → ReportLab paragraph markup."""
    from xml.sax.saxutils import escape
    t = escape(str(t))
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    return re.sub(r"(?<!\*)\*(.+?)\*(?!\*)", r"<i>\1</i>", t)


def build_pdf(res, author):
    """Create the PDF report (A4) and return its bytes. Every page carries the Swastome footer."""
    _pdf_fonts()
    import matplotlib
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib.utils import ImageReader
    from reportlab.platypus import Image, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
    from xml.sax.saxutils import escape

    with matplotlib.rc_context({"font.family": "DejaVu Sans", "font.size": 8.5, "text.color": "black", "axes.labelcolor": "black",
                                "axes.edgecolor": "black", "xtick.color": "black", "ytick.color": "black", "axes.titlecolor": "black",
                                "figure.facecolor": "white", "axes.facecolor": "white", "legend.labelcolor": "black"}):
        figs = _pdf_figures(res)

    best, P, corr = res["best"], res["P"], res["corr"]
    W, _H = A4
    LM = RM = 1.8 * cm
    avail = W - LM - RM
    ink, teal, grey = colors.HexColor("#111111"), colors.HexColor("#0b5d6b"), colors.HexColor("#555555")
    base = dict(fontName="DSans", textColor=ink)
    S = {
        "title": ParagraphStyle("t", fontName="DSans-Bold", fontSize=26, leading=30, textColor=teal),
        "sub": ParagraphStyle("s", fontSize=11, leading=15, textColor=grey, fontName="DSans"),
        "h1": ParagraphStyle("h1", fontName="DSans-Bold", fontSize=14, leading=18, textColor=teal, spaceBefore=6, spaceAfter=6),
        "h2": ParagraphStyle("h2", fontName="DSans-Bold", fontSize=10.5, leading=14, textColor=ink, spaceBefore=8, spaceAfter=4),
        "p": ParagraphStyle("p", fontSize=9, leading=13, spaceAfter=4, **base),
        "li": ParagraphStyle("li", fontSize=9, leading=13, leftIndent=12, bulletIndent=2, spaceAfter=2, **base),
        "cap": ParagraphStyle("cap", fontName="DSans-Italic", fontSize=7.8, leading=10, textColor=grey, alignment=TA_CENTER, spaceAfter=8),
        "th": ParagraphStyle("th", fontName="DSans-Bold", fontSize=7, leading=9, textColor=colors.white),
        "td": ParagraphStyle("td", fontName="DSans", fontSize=7, leading=9, textColor=ink),
        "small": ParagraphStyle("sm", fontName="DSans", fontSize=7.5, leading=10, textColor=grey, spaceAfter=4),
    }

    def para(t, st="p"):
        return Paragraph(_md_to_rl(t), S[st])

    def bullets(items):
        return [Paragraph(_md_to_rl(x), S["li"], bulletText="•") for x in items]

    def cellfmt(v):
        if isinstance(v, (float, np.floating)):
            if pd.isna(v):
                return "—"
            return f"{v:.3g}" if abs(v) >= 1e4 or (v != 0 and abs(v) < 1e-3) else f"{v:.2f}"
        return str(v)

    def table(df, max_rows=40):
        d = df.head(max_rows)
        head = [Paragraph(escape(str(c).replace("_", " ")), S["th"]) for c in d.columns]
        rows = [[Paragraph(escape(cellfmt(v)), S["td"]) for v in r] for r in d.itertuples(index=False)]
        lens = [min(28, max(6, len(str(c)), *(len(cellfmt(v)) for v in d[c]))) for c in d.columns]
        widths = [avail * l / sum(lens) for l in lens]
        t = Table([head] + rows, colWidths=widths, repeatRows=1)
        t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), teal), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                               ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#eef4f5")]),
                               ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c8d3d6")),
                               ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5)]))
        out = [t]
        if len(df) > max_rows:
            out.append(Paragraph(f"Showing the first {max_rows} of {len(df)} rows – the full table is in the ZIP package.", S["small"]))
        return out

    def table_t(title, df):
        """Heading + table kept together so a heading is never stranded at the bottom of a page."""
        return [KeepTogether([Paragraph(title, S["h2"])] + table(df))]

    def figure(key):
        if key not in figs:
            return []
        cap, png = figs[key]
        iw, ih = ImageReader(io.BytesIO(png)).getSize()
        w = avail
        h = w * ih / iw
        if h > 12.5 * cm:
            h = 12.5 * cm
            w = h * iw / ih
        return [KeepTogether([Image(io.BytesIO(png), width=w, height=h), Paragraph(escape(f"Figure {key}. {cap}"), S["cap"])])]

    top = best.loc[best["Best_Score"].idxmin()]
    story = [Paragraph(APP_NAME, S["title"]), Paragraph(APP_TAGLINE + " – analysis report", S["sub"]), Spacer(1, 10)]
    meta = [["Prepared by", author], ["Date", datetime.now().strftime("%d %B %Y, %H:%M")],
            ["Best-pocket definition", res["best_mode"]], ["Temperature for Ki", f"{res['temp_k']:.2f} K"],
            ["Complexes / pockets analysed", f"{len(best)} / {len(P)}"],
            ["Strongest complex", f"{top['Complex']} ({top['Best_Score']:.1f} kcal/mol, Ki ≈ {fmt_ki(top['Pred_Ki_nM'])})"]]
    mt = Table([[Paragraph(f"<b>{escape(a)}</b>", S["td"]), Paragraph(escape(b), S["td"])] for a, b in meta], colWidths=[avail * 0.3, avail * 0.7])
    mt.setStyle(TableStyle([("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#c8d3d6")), ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#eef4f5")),
                            ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
    story += [mt, Spacer(1, 10), Paragraph("Executive summary", S["h1"])] + bullets(narrative(res))

    cls = best["Affinity_Class"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int).reset_index()
    cls.columns = ["Affinity class", "Complexes"]
    story += table_t("Complexes per affinity class", cls)
    story += table_t("Top-5 complexes by priority", 
        best.sort_values("Priority_Rank").head(5)[["Priority_Rank", "Complex", "Best_Score", "Affinity_Class"]])

    story += [Spacer(1, 8), Paragraph("1. Binding affinity", S["h1"]),
              para(f"For each receptor–ligand file the best pocket ({res['best_mode'].split(' (')[0].lower()}) is selected and its AutoDock-Vina score "
                   "(kcal/mol) is reported. More negative = stronger predicted binding.")]
    story += figure("01") + figure("02")
    story += table_t("Best-pocket table", 
        best[["Priority_Rank", "Receptor", "Ligand", "Best_Pocket_ID", "Best_Score", "Affinity_Class", "Pred_Ki_nM"]].sort_values("Priority_Rank"))
    story += bullets(["**≤ −9 kcal/mol** very strong; **−9 to −7** strong; **−7 to −5** moderate; **> −5** weak / likely non-specific.",
                      "Differences below ~0.5 kcal/mol between ligands are within docking noise."])

    story += [Spacer(1, 8), Paragraph("2. Cavity volume vs binding score", S["h1"]),
              para("Each detected cavity is one point. Pearson (linear) and Spearman (rank) correlations test whether larger cavities systematically score better.")]
    o = corr[corr["Level"] == "Overall"].iloc[0]
    story += [para("**Key finding:** " + o["Interpretation"] + (f" (Pearson r = {o['Pearson_r']:.2f}, n = {int(o['N_pockets'])})." if not np.isnan(o["Pearson_r"]) else "."))]
    story += figure("03") + figure("04") + figure("05")
    story += table_t("Correlation statistics", 
        corr[["Level", "Group", "N_pockets", "Pearson_r", "Spearman_rho", "Slope_kcal_per_100A3", "Interpretation"]])

    story += [Spacer(1, 8), Paragraph("3. Pocket landscape", S["h1"]),
              para("Within each complex all pockets are compared: score spread, the dominance gap between the best and runner-up pocket, "
                   "and the number of competing pockets within 1 kcal/mol.")]
    story += figure("06") + figure("07") + figure("08") + figure("09")
    story += table_t("Dominance & competition metrics", 
        best[["Receptor", "Ligand", "Best_Score", "Second_Best_Score", "Gap_To_Second", "Competing_Pockets_1kcal", "Pocket_Breadth", "N_Pockets"]])
    story += bullets(["**Gap > 1.5 kcal/mol**: one clearly preferred site → likely specific binding.",
                      "**Small gap + many competing pockets**: the ligand can sit in several places; consider targeted docking per site."])

    story += [Spacer(1, 8), Paragraph("4. Binding site & residues", S["h1"]),
              para("Contact residues of each best pocket are profiled by chemical character. For receptors with ≥ 2 ligands, site overlap is quantified "
                   "by the Jaccard similarity of contact residues.")]
    story += figure("10") + figure("11")
    if len(res["jac"]):
        story += table_t("Ligand site-overlap table", 
            res["jac"][["Receptor", "Pair", "Jaccard_Similarity", "N_Shared_Residues", "Pocket_Centre_Distance_A", "Verdict"]])

    story += [Spacer(1, 8), Paragraph("5. Predictions & ranking", S["h1"]),
              para(f"Vina scores are converted into estimated Ki (Ki ≈ exp(ΔG / RT) at {res['temp_k']:.1f} K); ligands are compared head-to-head (ΔΔG); "
                   "every complex receives a composite priority score = 70% binding strength + 15% site dominance + 15% pocket breadth.")]
    story += figure("12") + figure("13")
    story += table_t("Ranking table", 
        best[["Priority_Rank", "Complex", "Priority_Score", "Best_Score", "Pred_Ki_nM", "Pred_Ki_uM", "Affinity_Class"]].sort_values("Priority_Rank"))
    if len(res["pairs"]):
        story += figure("14")
        story += table_t("Ligand-vs-ligand comparison", 
            res["pairs"][["Receptor", "Pair", "Score_A", "Score_B", "DDG_B_minus_A", "Stronger_Binder", "Ki_Fold_Difference", "Verdict"]])

    story += [Spacer(1, 8), Paragraph("Methods text (suggested)", S["h2"]),
              para(f"Blind docking was performed with CB-Dock2 (cavity detection + AutoDock Vina). {len(P)} cavities across {len(best)} receptor–ligand complexes "
                   f"were analysed. The best pocket per complex was defined as the {res['best_mode'].split(' (')[0].lower()}. Estimated Ki was derived as "
                   f"exp(ΔG/RT) at {res['temp_k']:.1f} K. Contact residues, pocket dominance and cavity-volume correlations were evaluated with {APP_NAME}."),
              Paragraph("Before you report", S["h2"])]
    story += bullets(["Docking scores are predictions; validate top hits with MD simulation or experiments.",
                      "Keep the same number of cavities for all files in one analysis.",
                      "Priority score and estimated Ki are screening heuristics; Vina typically carries ~2–3 kcal/mol error."])
    story += [Spacer(1, 6), para("Please cite the underlying tools: Liu et al., CB-Dock2, Nucleic Acids Res. 2022; Trott & Olson, AutoDock Vina, "
                                 "J Comput Chem 2010; Eberhardt et al., AutoDock Vina 1.2.0, J Chem Inf Model 2021.", "small")]

    fstyle = ParagraphStyle("f", fontName="DSans", fontSize=6.4, leading=8, textColor=grey, alignment=TA_CENTER)

    def decorate(canvas, doc):
        canvas.saveState()
        pw, ph = A4
        canvas.setFont("DSans", 7.5)
        canvas.setFillColor(grey)
        canvas.drawString(LM, ph - 1.05 * cm, f"{APP_NAME} analysis report")
        canvas.drawRightString(pw - RM, ph - 1.05 * cm, f"Prepared by: {author}"[:90])
        canvas.setStrokeColor(colors.HexColor("#c8d3d6"))
        canvas.setLineWidth(0.5)
        canvas.line(LM, ph - 1.25 * cm, pw - RM, ph - 1.25 * cm)
        y = 1.15 * cm
        for line in reversed(FOOTER_LINES):
            p = Paragraph(escape(line), fstyle)
            _w, h = p.wrap(avail, 100)
            p.drawOn(canvas, LM, y)
            y += h + 2
        canvas.line(LM, y + 1, pw - RM, y + 1)
        canvas.setFont("DSans", 7.5)
        canvas.drawCentredString(pw / 2, 0.55 * cm, f"Page {canvas.getPageNumber()}")
        canvas.restoreState()

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, leftMargin=LM, rightMargin=RM, topMargin=1.8 * cm, bottomMargin=3.1 * cm,
                            title=f"{APP_NAME} analysis report", author=author, creator=APP_NAME)
    doc.build(story, onFirstPage=decorate, onLaterPages=decorate)
    return buf.getvalue()


def tab_report(res):
    best, P = res["best"], res["P"]
    st.markdown("### 7️⃣ Overall results & download")
    block("What is happening", "This tab condenses the whole analysis into headline numbers and a plain-language summary, and packages every table "
          "and figure into one ZIP file for your records, supplementary material or sharing with collaborators.")
    top = best.loc[best["Best_Score"].idxmin()]
    c = st.columns(5)
    c[0].metric("Complexes analysed", len(best))
    c[1].metric("Pockets analysed", len(P))
    c[2].metric("Best score", f"{top['Best_Score']:.1f} kcal/mol", top["Complex"], delta_color="off")
    c[3].metric("Mean best score", f"{best['Best_Score'].mean():.2f}")
    c[4].metric("Strong / very strong", int(best["Affinity_Class"].isin(CLASS_ORDER[:2]).sum()))
    st.markdown("#### 🧾 Executive summary")
    for x in narrative(res):
        st.markdown(f"- {x}")
    cls = best["Affinity_Class"].value_counts().reindex(CLASS_ORDER).fillna(0).reset_index()
    cls.columns = ["Affinity_Class", "N_complexes"]
    c1, c2 = st.columns([1, 1])
    with c1:
        show_table(cls, "affinity_class_counts", "Complexes per affinity class")
    with c2:
        show_table(best.sort_values("Priority_Rank").head(5)[["Priority_Rank", "Complex", "Best_Score", "Affinity_Class"]], "top5_complexes", "Top-5 complexes")
    st.markdown("#### ✍️ Suggested methods text")
    st.code(f"Blind docking was performed with CB-Dock2 (cavity detection + AutoDock Vina). {len(P)} cavities across {len(best)} receptor–ligand "
            f"complexes were analysed. The best pocket per complex was defined as the {res['best_mode'].split(' (')[0].lower()}. Estimated Ki was derived "
            f"as exp(ΔG/RT) at {res['temp_k']:.1f} K. Contact residues, pocket dominance and cavity-volume correlations were evaluated with DockScope.", language="text")
    st.caption("Please cite the underlying tools: Liu et al., CB-Dock2, Nucleic Acids Res. 2022; Trott & Olson, AutoDock Vina, J Comput Chem 2010; Eberhardt et al., AutoDock Vina 1.2.0, J Chem Inf Model 2021.")
    block("Before you report", "- Docking scores are predictions; validate top hits with MD simulation or experiments.\n"
          "- Keep the same number of cavities for all files in one analysis.\n- Priority score and Ki are screening heuristics.", "warn")
    st.markdown("#### 📦 Download everything")
    if st.button("🗜️ Build ZIP package (tables + interactive plots + summary)", type="primary"):
        bar = st.progress(0, text="Preparing…")
        data, nfig = build_zip(res, lambda f, t: bar.progress(min(f, 1.0), text=t))
        st.session_state.zip_bytes = data
        st.session_state.zip_info = f"{nfig} interactive HTML plots, CSV tables, summary"
    if st.session_state.zip_bytes:
        st.success("Package ready: " + st.session_state.zip_info)
        st.download_button("⬇️ Download DockScope_results.zip", st.session_state.zip_bytes, file_name="DockScope_results.zip", mime="application/zip")

    st.markdown("#### 📄 Download PDF report")
    st.caption("A formatted A4 report with the executive summary, all 14 figures (black text on white), key tables, methods text and the "
               "Swastome footer on every page. The name you enter is printed on the report header of every page.")
    author = st.text_input("Your name (printed on the PDF)", key="pdf_author", placeholder="e.g. Dr. Jane Doe").strip()
    if st.button("📄 Generate PDF report", disabled=not author, help="Enter your name first" if not author else None):
        try:
            with st.spinner("Building PDF – rendering figures, this takes a few seconds…"):
                st.session_state.pdf_bytes = build_pdf(res, author)
                st.session_state.pdf_for = author
        except ImportError as e:
            st.session_state.pdf_bytes = None
            st.error(f"PDF libraries are missing ({e}). Add `reportlab` and `matplotlib` to requirements.txt and restart the app.")
        except Exception as e:  # noqa: BLE001
            st.session_state.pdf_bytes = None
            st.error(f"Could not build the PDF: {e}")
    if st.session_state.pdf_bytes:
        if st.session_state.get("pdf_for") == author:
            fname = re.sub(r"[^A-Za-z0-9_-]+", "_", author).strip("_") or "report"
            st.success("PDF ready.")
            st.download_button("⬇️ Download DockScope_report.pdf", st.session_state.pdf_bytes,
                               file_name=f"DockScope_report_{fname}.pdf", mime="application/pdf", type="primary")
        else:
            st.info("The name was changed – click **Generate PDF report** again to print the new name.")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
with st.sidebar:
    st.markdown(f"## 🧬 {APP_NAME}")
    st.caption(APP_TAGLINE)
    st.selectbox("Best-pocket definition", [BEST_LOWEST, BEST_ID1], key="best_mode")
    st.number_input("Temperature for Ki estimate (K)", min_value=250.0, max_value=400.0, value=298.15, step=0.5, key="temp_k")
    st.divider()
    if st.button("♻️ Reset everything"):
        for k in list(st.session_state.keys()):
            del st.session_state[k]
        st.rerun()
    st.caption("Tip: changing these settings after an analysis requires clicking **Run prediction analysis** again.")

st.markdown(f'<p class="ds-title">🧬 {APP_NAME}</p><p class="ds-sub">{APP_TAGLINE}</p>', unsafe_allow_html=True)
st.markdown(
    "**DockScope turns raw CB-Dock2 `curpocket.txt` outputs into a publication-ready docking report.** Upload files one-by-one or in bulk, "
    "let the tool verify every file, then explore binding affinity, cavity–score relationships, pocket landscapes, contact-residue chemistry, "
    "estimated Ki and ranking predictions – all with interactive Plotly figures, downloadable tables and a one-click ZIP package."
)
st.info(
    "🔗 **Input source:** DockScope works on the output of **CB-Dock2** "
    "([cao.labshare.cn/cb-dock2](http://cao.labshare.cn:11080/cb-dock2/php/blinddock.php)). Run your blind docking there, "
    "then download the **`curpocket.txt`** file from the results page and upload it here."
)

S = current_state()
tabs = st.tabs(["1️⃣ Upload & Validate", "2️⃣ Binding Affinity", "3️⃣ Volume vs Score", "4️⃣ Pocket Landscape",
                "5️⃣ Site & Residues", "6️⃣ Predictions", "7️⃣ Report & Download"])

with tabs[0]:
    render_banner(S)
    tab_upload(S)
    S = current_state()  # refresh after possible changes in this run

for i, fn in enumerate([tab_affinity, tab_volume, tab_landscape, tab_residues, tab_predictions, tab_report], start=1):
    with tabs[i]:
        render_banner(S)
        res = gate(S)
        if res is not None:
            fn(res)

# ----------------------------------------------------------------------------
# App footer
# ----------------------------------------------------------------------------
st.divider()
st.markdown('<div class="ds-footer">' + "<br>".join(FOOTER_LINES) + "</div>", unsafe_allow_html=True)
