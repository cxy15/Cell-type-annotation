#!/usr/bin/env python3
"""Step 4: Apply LLM annotations, UMAP embedding plot; overwrite working h5ad."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import omicverse as ov
import pandas as pd
import scanpy as sc

from pipeline_common import exchange_dir, figures_dir, get_out_dir, merge_progress, setup_logging, working_adata_path

ov.plot_set()


def load_annotation(path: Path) -> tuple[dict, dict, dict, str]:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return (
        data["celltype_to_plot"],
        data["markers_to_plot"],
        data["color_dict"],
        data.get("tissue", ""),
    )


def main() -> None:
    log = setup_logging("04_plot")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--adata",
        default=None,
        help="Working adata from step 1 (default: .workflow_adata.h5ad under output dir)",
    )
    args = ap.parse_args()

    if args.out:
        os.environ["CELLTYPE_PIPELINE_OUT"] = str(Path(args.out).resolve())

    out = get_out_dir()
    ex = exchange_dir()
    figd = figures_dir()

    ann_path = ex / "llm_annotation.json"
    if not ann_path.exists():
        raise FileNotFoundError(f"Run step 3 first: {ann_path}")
    adata_path = Path(args.adata) if args.adata else working_adata_path()
    if not adata_path.exists():
        raise FileNotFoundError(adata_path)

    celltype_to_plot, _markers_to_plot, color_dict, _tissue = load_annotation(ann_path)
    celltype_to_plot = {str(k): v for k, v in celltype_to_plot.items()}

    adata = sc.read_h5ad(adata_path)

    leiden_str = adata.obs["leiden"].astype(str)
    adata.obs["celltype_to_plot"] = leiden_str.map(celltype_to_plot)

    cats = list(color_dict.keys())
    adata.obs["celltype_to_plot"] = pd.Categorical(adata.obs["celltype_to_plot"], categories=cats, ordered=True)

    fig_umap, ax = plt.subplots(figsize=(5, 5))
    ov.pl.embedding(
        adata,
        basis="X_umap",
        color=["celltype_to_plot"],
        palette=color_dict,
        frameon="small",
        show=False,
        ax=ax,
    )
    umap_png = figd / "umap_celltype.png"
    fig_umap.savefig(umap_png, dpi=300, bbox_inches="tight")
    plt.close(fig_umap)
    log.info("Saved %s", umap_png)

    work_path = working_adata_path()
    adata.write_h5ad(work_path)
    log.info("Updated working AnnData %s", work_path)

    merge_progress(
        {
            "step": 4,
            "status": "done",
            "working_adata": str(work_path),
            "figures": {"umap": str(umap_png)},
        }
    )
    log.info("Step 4 complete.")


if __name__ == "__main__":
    main()
