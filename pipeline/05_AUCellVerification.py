#!/usr/bin/env python3
"""Step 5: LLM gold-standard markers + Omicverse geneset_aucell; write single final h5ad."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from pipeline_common import (
    exchange_dir,
    final_adata_path,
    get_llm_config,
    get_out_dir,
    get_species,
    merge_progress,
    save_json,
    setup_logging,
    slug_geneset_name,
    working_adata_path,
)

import numpy as np
import pandas as pd
import scanpy as sc
from openai import OpenAI

import omicverse as ov

ov.settings.cpu_gpu_mixed_init()


def species_display() -> str:
    return "小鼠" if get_species() == "mouse" else "人类"


def fetch_gold_standard_genes(
    cell_type: str,
    client: OpenAI,
    model: str,
    temperature: float,
    timeout: int,
) -> list[str]:
    sp = species_display()
    prompt = f"""
你是单细胞与流式细胞术领域的专家。请给出{sp}组织中「{cell_type}」这一细胞类型在文献与数据库中常用的权威（金标准）marker 基因，用于鉴定该类型。

只输出 JSON，格式严格如下：
{{"genes": ["GENE1", "GENE2", "..."]}}

要求：
- 8–20 个基因符号（SYMBOL），大写人类基因、小鼠基因按常规命名；
- 优先经典、共识度高的 marker；
- 不要解释，不要 Markdown，不要任何其他键。
"""
    r = client.chat.completions.create(
        model=model,
        temperature=temperature,
        timeout=timeout,
        messages=[{"role": "user", "content": prompt}],
    )
    text = r.choices[0].message.content.strip()
    text = text.replace("```json", "").replace("```", "").strip()
    data = json.loads(text)
    genes = data.get("genes", [])
    out = []
    for g in genes:
        if not isinstance(g, str):
            continue
        g = g.strip()
        if g:
            out.append(g)
    return out


def filter_genes_to_var(genes: list[str], var_names: pd.Index) -> list[str]:
    vn = set(var_names.astype(str))
    vn_lower = {x.lower(): x for x in vn}
    kept: list[str] = []
    for g in genes:
        if g in vn:
            kept.append(g)
        elif g.lower() in vn_lower:
            kept.append(vn_lower[g.lower()])
    seen = set()
    uniq = []
    for g in kept:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq


def _leiden_sort_key(val: object) -> tuple:
    try:
        return (0, int(float(str(val))))
    except (ValueError, TypeError):
        return (1, str(val))


def print_leiden_gold_standard_report(obs: pd.DataFrame) -> None:
    """Print per-Leiden pass/fail for AUCell verification (stdout, visible in bash logs)."""
    if "leiden" not in obs.columns:
        print("\n=== Leiden 金标准验证: obs 中无 leiden 列 ===\n")
        return
    if "aucell_match_argmax" not in obs.columns:
        print("\n=== Leiden 金标准验证: 无 AUCell 结果 ===\n")
        return

    leiden_vals = list(obs["leiden"].unique())
    leiden_vals.sort(key=_leiden_sort_key)

    ma = obs["aucell_match_argmax"].fillna(False).astype(bool)
    has_peak = "aucell_peak_for_assigned_gs" in obs.columns
    pk = obs["aucell_peak_for_assigned_gs"].fillna(False).astype(bool) if has_peak else None

    print("")
    print("========== Leiden 金标准验证（按 cluster）==========")
    print("  aucell_match_argmax: 各金标准 AUCell 列 argmax 是否等于当前注释")
    if has_peak:
        print("  aucell_peak_for_assigned_gs: 所属类型金标准列上是否为全局最高分")
    print("")

    all_pass_ma = True
    all_pass_pk = True

    for lv in leiden_vals:
        m = obs["leiden"] == lv
        n = int(m.sum())
        if n == 0:
            continue
        n_ma = int(ma[m].sum())
        pass_ma = n_ma == n
        if not pass_ma:
            all_pass_ma = False
        st_ma = "PASS" if pass_ma else "FAIL"

        if pk is not None:
            n_pk = int(pk[m].sum())
            pass_pk = n_pk == n
            if not pass_pk:
                all_pass_pk = False
            st_pk = "PASS" if pass_pk else "FAIL"
            print(f"  leiden {lv}:  argmax {st_ma} ({n_ma}/{n})  |  peak {st_pk} ({n_pk}/{n})")
        else:
            print(f"  leiden {lv}:  argmax {st_ma} ({n_ma}/{n})")

    print("")
    print(f"  汇总 — 全部 Leiden 在 argmax 规则下均通过: {'是' if all_pass_ma else '否'}")
    if pk is not None:
        print(f"  汇总 — 全部 Leiden 在 peak 规则下均通过: {'是' if all_pass_pk else '否'}")
    print("====================================================")
    print("")


def main() -> None:
    log = setup_logging("05_aucell")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--adata",
        default=None,
        help="Working annotated adata (default: .workflow_adata.h5ad)",
    )
    ap.add_argument("--temperature", type=float, default=0.1)
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument(
        "--reuse-gold-json",
        action="store_true",
        help="Load Exchange/gold_standard_markers.json instead of calling LLM",
    )
    args = ap.parse_args()

    if args.out:
        os.environ["CELLTYPE_PIPELINE_OUT"] = str(Path(args.out).resolve())

    out = get_out_dir()
    ex = exchange_dir()

    adata_path = Path(args.adata) if args.adata else working_adata_path()
    if not adata_path.exists():
        raise FileNotFoundError(adata_path)

    adata = sc.read_h5ad(adata_path)
    if "celltype_to_plot" not in adata.obs.columns:
        raise KeyError("adata.obs must contain celltype_to_plot (run step 4 first)")

    assigned = adata.obs["celltype_to_plot"].astype(str)
    unique_types = [x for x in pd.unique(assigned) if x and str(x) not in ("nan", "Others", "Other")]
    if not unique_types:
        log.debug("No non-Others cell types; skip AUCell.")
        final_path = final_adata_path()
        adata.write_h5ad(final_path)
        log.info("Wrote final (no AUCell) %s", final_path)
        merge_progress({"step": 5, "status": "skipped", "reason": "only Others/empty", "final_h5ad": str(final_path)})
        wp = working_adata_path()
        if wp.exists() and wp.resolve() != final_path.resolve():
            wp.unlink()
            log.info("Removed working file %s", wp)
        print("")
        print("========== Leiden 金标准验证 ==========")
        print("  (已跳过: 无非 Others 细胞类型，未计算 AUCell)")
        print("========================================")
        print("")
        return

    gold_path = ex / "gold_standard_markers.json"
    cell_type_to_slug: dict[str, str] = {}
    used_slugs: set[str] = set()

    if args.reuse_gold_json and gold_path.exists():
        log.info("Loading gold standard genes from %s", gold_path)
        with open(gold_path, encoding="utf-8") as f:
            gold_payload = json.load(f)
        per_ct = gold_payload.get("per_cell_type", gold_payload)
    elif args.reuse_gold_json:
        raise FileNotFoundError(f"--reuse-gold-json set but missing {gold_path}")
    else:
        cfg = get_llm_config()
        client = OpenAI(api_key=cfg["api_key"], base_url=cfg["base_url"], timeout=args.timeout)
        per_ct: dict[str, dict] = {}
        for ct in unique_types:
            log.info("LLM: gold-standard markers for %s", ct)
            genes = fetch_gold_standard_genes(
                ct, client, cfg["model"], args.temperature, args.timeout
            )
            slug = slug_geneset_name(ct, used_slugs)
            cell_type_to_slug[ct] = slug
            filtered = filter_genes_to_var(genes, adata.var_names)
            per_ct[ct] = {
                "genes_requested": genes,
                "genes_in_matrix": filtered,
                "geneset_name_slug": slug,
            }
            if len(filtered) < 5:
                log.debug("Few genes in matrix for %s: %s", ct, len(filtered))
        save_json(
            gold_path,
            {"species": get_species(), "per_cell_type": per_ct},
        )
        log.info("Wrote %s", gold_path)

    if args.reuse_gold_json and gold_path.exists():
        for ct, info in per_ct.items():
            slug = info.get("geneset_name_slug") or slug_geneset_name(ct, used_slugs)
            cell_type_to_slug[ct] = slug

    for ct in unique_types:
        info = per_ct.get(ct, {})
        slug = info.get("geneset_name_slug") or slug_geneset_name(ct, used_slugs)
        cell_type_to_slug[ct] = slug
        genes = info.get("genes_in_matrix") or filter_genes_to_var(
            info.get("genes_requested", []), adata.var_names
        )
        if not genes:
            log.debug("Skipping AUCell for %s: no genes in matrix", ct)
            continue
        log.info("geneset_aucell: %s -> %s (%s genes)", ct, slug, len(genes))
        ov.single.geneset_aucell(adata, geneset_name=slug, geneset=genes)

    slug_list = [cell_type_to_slug[c] for c in unique_types if c in cell_type_to_slug]
    obs = adata.obs
    score_cols = {sl: f"{sl}_aucell" for sl in slug_list if f"{sl}_aucell" in obs.columns}

    if not score_cols:
        raise RuntimeError("No *_aucell columns in obs — check geneset_aucell step.")

    mat = np.column_stack([obs[c].values.astype(float) for c in score_cols.values()])
    col_names = list(score_cols.keys())

    best_idx = np.argmax(mat, axis=1)
    best_slug = np.array([col_names[i] for i in best_idx], dtype=object)

    assigned_slug = np.array(
        [cell_type_to_slug.get(str(obs["celltype_to_plot"].iloc[i]), "") for i in range(adata.n_obs)],
        dtype=object,
    )

    match_argmax = (best_slug == assigned_slug) & (assigned_slug != "")

    peak_for_assigned = np.zeros(adata.n_obs, dtype=bool)
    score_on_assigned = np.full(adata.n_obs, np.nan)
    for i in range(adata.n_obs):
        sl = assigned_slug[i]
        if sl == "" or sl not in score_cols:
            continue
        col = score_cols[sl]
        vals = obs[col].values.astype(float)
        score_on_assigned[i] = vals[i]
        peak_for_assigned[i] = vals[i] >= np.max(vals) - 1e-12

    obs["aucell_best_slug"] = best_slug
    obs["aucell_match_argmax"] = match_argmax
    obs["aucell_peak_for_assigned_gs"] = peak_for_assigned
    obs["aucell_score_on_assigned_gs"] = score_on_assigned

    bc = pd.Index(adata.obs_names.astype(str))
    summary = pd.DataFrame(
        {
            "cell": bc,
            "celltype_to_plot": obs["celltype_to_plot"].astype(str).values,
            "assigned_gs_slug": assigned_slug,
            "best_gs_slug": best_slug,
            "aucell_match_argmax": match_argmax,
            "aucell_peak_for_assigned_gs": peak_for_assigned,
            "aucell_score_on_assigned_gs": score_on_assigned,
        }
    )
    for sl, col in score_cols.items():
        summary[f"aucell_{sl}"] = obs[col].values

    ver_csv = ex / "aucell_verification.csv"
    summary.to_csv(ver_csv, index=False)
    log.info("Wrote %s", ver_csv)

    print_leiden_gold_standard_report(adata.obs)

    final_path = final_adata_path()
    adata.write_h5ad(final_path)
    log.info("Wrote final annotated h5ad %s", final_path)

    wp = working_adata_path()
    if wp.exists() and wp.resolve() != final_path.resolve():
        wp.unlink()
        log.info("Removed working file %s", wp)

    merge_progress(
        {
            "step": 5,
            "status": "done",
            "gold_standard_markers": str(gold_path),
            "aucell_verification_csv": str(ver_csv),
            "final_h5ad": str(final_path),
        }
    )
    log.info("Step 5 complete.")


if __name__ == "__main__":
    main()
