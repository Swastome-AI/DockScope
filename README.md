# 🧬 DockScope – CB-Dock2 Blind-Docking Pocket Analytics Suite

DockScope v2.0 is live and available - https://dockscope-v2.streamlit.app/

Turn raw CB-Dock2 `curpocket.txt` files into a full docking report: validated inputs, interactive Plotly figures,
CSV tables, estimated Ki, ranking and a one-click ZIP package.

**Input source:** DockScope works on output from CB-Dock2 – http://cao.labshare.cn:11080/cb-dock2/php/blinddock.php . Run blind docking there and download the `curpocket.txt` file from the results page.

## Features
- **Manual upload**: register receptor/ligand names once, pick them from drop-downs, declare the number of cavities, upload one file at a time. Each file is checked (columns, ID count = declared cavities, duplicates). Failed files are flagged red and must be removed and re-uploaded.
- **Bulk upload**: many `.txt` files + one metadata sheet (`Ligand, Receptor, Number_of_cavities, File_name`). Exact file-name matching, spelling-variant detection, ID-count check, "100% successful" gate.
- **Locked workflow**: tabs 2–7 stay locked until every file is green and the analysis has been run (progress bar).
- **Analyses**: best-pocket affinity, heatmap, volume-vs-score correlation, pocket landscape (dominance gap, competing sites, outliers), residue chemistry, ligand site overlap, 3-D pocket map, estimated Ki, ΔΔG ligand selectivity, composite priority ranking.
- **Downloads**: CSV for every table, HTML/4x PNG for every plot (toolbar camera), and a final ZIP (tables + HTML + PNG + summary).

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy via GitHub + Streamlit Community Cloud
1. Create a GitHub repo and push `app.py`, `requirements.txt`, `.streamlit/config.toml`, `README.md`.
2. Go to https://share.streamlit.io → **New app** → pick the repo, branch `main`, main file `app.py` → **Deploy**.
3. If PNG export inside the ZIP fails on the server, the ZIP still contains HTML plots and tables (kaleido 0.2.1 is pinned to avoid needing Chrome).

## Metadata template
```
Ligand,Receptor,Number_of_cavities,File_name
Daidzein,KatE,5,Daidzein - Catalase.txt
Genistein,KatE,5,Genistein - Catalase.txt
```

## Notes
- Use the **same number of cavities** for every file in one analysis.
- CB-Dock2 orders pockets by cavity size, so ID 1 is not necessarily the best score; default best pocket = lowest Vina score (switchable in the sidebar).
- Docking scores are predictions; Ki and priority score are screening heuristics.
