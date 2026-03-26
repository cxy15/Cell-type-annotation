#!/usr/bin/env python3
"""Step 1: HVG subset, markers -> Exchange/all_markers.csv + .workflow_adata.h5ad.

Assumes input AnnData already contains X_umap (and typically neighbors/leiden) computed on the
same highly-variable gene space; this step does not compute UMAP.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from pipeline_common import exchange_dir, get_out_dir, merge_progress, setup_logging, working_adata_path

import omicverse as ov
import pandas as pd
import scanpy as sc

ov.settings.cpu_gpu_mixed_init()


def main() -> None:
    log = setup_logging("01_get_markers")
    ap = argparse.ArgumentParser(description="Extract cluster markers for downstream enrichment / LLM.")
    ap.add_argument(
        "--input",
        default=os.environ.get("INPUT_H5AD", "SC.h5ad"),
        help="Input h5ad path (default: INPUT_H5AD env or SC.h5ad)",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output root directory (default: CELLTYPE_PIPELINE_OUT or <repo>/outputs)",
    )
    ap.add_argument("--n-top-hvg", type=int, default=2000)
    ap.add_argument("--leiden-resolution", type=float, default=0.5)
    ap.add_argument("--foldchange", type=float, default=2.0)
    ap.add_argument("--top-genes", type=int, default=50, help="Markers per cluster from rank_genes_groups")
    args = ap.parse_args()

    if args.out:
        os.environ["CELLTYPE_PIPELINE_OUT"] = str(Path(args.out).resolve())
    out = get_out_dir()
    ex = exchange_dir()
    log.info("Output directory: %s", out)

    adata = sc.read_h5ad(args.input)

    if "leiden" not in adata.obs.columns:
        log.info("Running Leiden clustering (resolution=%s)", args.leiden_resolution)
        if "neighbors" not in adata.uns:
            sc.pp.neighbors(adata)
        sc.tl.leiden(adata, resolution=args.leiden_resolution)

    if "highly_variable" not in adata.var.columns:
        log.info("Computing highly variable genes (n_top_genes=%s)", args.n_top_hvg)
        ov.pp.highly_variable_genes(adata, n_top_genes=args.n_top_hvg)

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"]]

    log.info("Finding marker genes per Leiden cluster")
    marker_dict = ov.single.get_celltype_marker(
        adata,
        clustertype="leiden",
        rank=True,
        key="rank_genes_groups",
        foldchange=args.foldchange,
        topgenenumber=args.top_genes,
    )

    records = []
    for celltype, genes in marker_dict.items():
        for rank, gene in enumerate(genes, 1):
            records.append({"celltype": celltype, "gene": gene, "rank": rank})
    markers_path = ex / "all_markers.csv"
    pd.DataFrame(records).to_csv(markers_path, index=False)
    log.info("Wrote %s", markers_path)

    work_path = working_adata_path()
    adata.write_h5ad(work_path)
    log.info("Wrote working AnnData %s (removed after pipeline completes)", work_path)

    merge_progress(
        {
            "step": 1,
            "status": "done",
            "markers_csv": str(markers_path),
            "working_adata": str(work_path),
        }
    )
    log.info("Step 1 complete.")


if __name__ == "__main__":
    main()
