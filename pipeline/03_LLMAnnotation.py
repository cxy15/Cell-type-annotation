#!/usr/bin/env python3
"""Step 3: LLM annotation from markers + enrichment -> Exchange/llm_annotation.json."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from pipeline_common import exchange_dir, get_llm_config, merge_progress, save_json, setup_logging

import pandas as pd
from openai import OpenAI


def load_marker_dict(markers_csv: Path) -> dict[str, list[str]]:
    df = pd.read_csv(markers_csv)
    df = df.sort_values(["celltype", "rank"])
    out: dict[str, list[str]] = {}
    for ct, sub in df.groupby("celltype"):
        out[str(ct)] = sub["gene"].astype(str).tolist()
    return out


def load_pathway_dict(enrichment_csv: Path) -> dict[str, list[str]]:
    df = pd.read_csv(enrichment_csv)
    if "Cluster" not in df.columns:
        raise ValueError(f"Missing Cluster column in {enrichment_csv}")
    pathway_dict: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        cid = str(row["Cluster"])
        paths = []
        for j in range(1, 6):
            c = f"Top{j}_Pathway"
            if c in df.columns and pd.notna(row.get(c, None)):
                paths.append(str(row[c]))
        pathway_dict[cid] = paths
    return pathway_dict


def annotate_cells(
    tissue: str,
    marker_dict: dict[str, list[str]],
    pathway_dict: dict[str, list[str]],
    model: str,
    base_url: str,
    api_key: str,
    temperature: float = 0.2,
    timeout: int = 120,
):
    client = OpenAI(api_key=api_key, base_url=base_url, timeout=timeout)

    prompt = f"""
我正在进行一个单细胞数据注释，请你帮我完成。下面是详细信息：

组织：{tissue}

以下是Leiden识别的marker(存在顺序关系，靠前的基因更显著）：
{json.dumps(marker_dict, ensure_ascii=False, indent=2)}

以下是Leiden识别的通路(存在顺序关系，靠前的通路更显著）：
{json.dumps(pathway_dict, ensure_ascii=False, indent=2)}

注释要求：
1、选择最可信的细胞类型，如果细胞为低质量细胞、或者难以分类的细胞，统一注释为‘Others’，我会删除这些细胞
2、我们进行大群注释，请仅注释粗略的细胞类型，如仅注释出T细胞，不要细化为CD8T细胞等
3、请你根据如下格式给出三个字典给我使用，分别是：
   - leiden->细胞类型
   - 细胞类型->4个最具特征marker
   - 细胞类型->绘图颜色

只允许输出JSON：

{{
"celltype_to_plot": {{}},
"markers_to_plot": {{}},
"color_dict": {{}}
}}

不要输出解释，不要Markdown。
"""

    response = client.chat.completions.create(
        model=model,
        temperature=temperature,
        messages=[{"role": "user", "content": prompt}],
    )

    text = response.choices[0].message.content.strip()
    text = text.replace("```json", "").replace("```", "").strip()
    result = json.loads(text)
    return result["celltype_to_plot"], result["markers_to_plot"], result["color_dict"]


def main() -> None:
    log = setup_logging("03_llm")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None, help="CELLTYPE_PIPELINE_OUT override")
    ap.add_argument(
        "--tissue",
        default=None,
        help="Tissue description (required for each run; or set TISSUE env from run_workflow.sh)",
    )
    args = ap.parse_args()

    tissue = (args.tissue or os.environ.get("TISSUE", "")).strip()
    if not tissue:
        log.error("Missing tissue: pass --tissue or set TISSUE (run_workflow.sh sets it).")
        sys.exit(1)

    if args.out:
        os.environ["CELLTYPE_PIPELINE_OUT"] = str(Path(args.out).resolve())

    ex = exchange_dir()
    cfg = get_llm_config()

    markers_path = ex / "all_markers.csv"
    enrich_path = ex / "enrichment_summary.csv"
    if not markers_path.exists():
        raise FileNotFoundError(markers_path)
    if not enrich_path.exists():
        raise FileNotFoundError(enrich_path)

    marker_dict = load_marker_dict(markers_path)
    pathway_dict = load_pathway_dict(enrich_path)

    log.info("Calling LLM model=%s", cfg["model"])
    celltype_to_plot, markers_to_plot, color_dict = annotate_cells(
        tissue=tissue,
        marker_dict=marker_dict,
        pathway_dict=pathway_dict,
        model=cfg["model"],
        base_url=cfg["base_url"],
        api_key=cfg["api_key"],
    )

    payload = {
        "tissue": tissue,
        "celltype_to_plot": {str(k): v for k, v in celltype_to_plot.items()},
        "markers_to_plot": {str(k): v for k, v in markers_to_plot.items()},
        "color_dict": {str(k): v for k, v in color_dict.items()},
    }
    out_path = ex / "llm_annotation.json"
    save_json(out_path, payload)
    log.info("Wrote %s", out_path)

    merge_progress(
        {
            "step": 3,
            "status": "done",
            "llm_annotation": str(out_path),
        }
    )
    log.info("Step 3 complete.")


if __name__ == "__main__":
    main()
