#!/usr/bin/env bash
# Full pipeline: 01 markers -> 02 enrichment (R) -> 03 LLM -> 04 plots -> 05 AUCell
# Requires: each run: --species and --tissue. LLM vars from .env or interactive prompts.
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$PIPELINE_DIR/.." && pwd)"
cd "$PIPELINE_DIR"

ENV_FILE="${ENV_FILE:-$REPO_ROOT/.env}"
if [[ -f "$ENV_FILE" ]]; then
  set -a
  # shellcheck source=/dev/null
  source "$ENV_FILE"
  set +a
fi

if [[ -z "${LLM_BASE_URL:-}" ]]; then
  read -r -p "LLM_BASE_URL: " LLM_BASE_URL
  export LLM_BASE_URL
fi
if [[ -z "${LLM_MODEL:-}" ]]; then
  read -r -p "LLM_MODEL: " LLM_MODEL
  export LLM_MODEL
fi
if [[ ! -f "$ENV_FILE" ]] || ! grep -qE '^[[:space:]]*LLM_API_KEY=' "$ENV_FILE" 2>/dev/null; then
  read -r -p "LLM_API_KEY (optional): " LLM_API_KEY
  export LLM_API_KEY
else
  export LLM_API_KEY="${LLM_API_KEY:-}"
fi

if [[ -z "${LLM_BASE_URL:-}" ]] || [[ -z "${LLM_MODEL:-}" ]]; then
  echo "[run] ERROR: LLM_BASE_URL 与 LLM_MODEL 不能为空。" >&2
  exit 1
fi

SPECIES=""
TISSUE=""
INPUT_H5AD_ARG=""
OUT_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species)
      SPECIES="$2"
      shift 2
      ;;
    --tissue)
      TISSUE="$2"
      shift 2
      ;;
    --input)
      INPUT_H5AD_ARG="$2"
      shift 2
      ;;
    --out)
      OUT_OVERRIDE="$2"
      shift 2
      ;;
    -h|--help)
      echo "Usage: $0 --species human|mouse --tissue \"组织描述\" [--input path/to.h5ad] [--out outputs_dir]"
      exit 0
      ;;
    *)
      echo "Unknown option: $1 (use --help)" >&2
      exit 1
      ;;
  esac
done

if [[ -z "$SPECIES" ]]; then
  read -r -p "物种请选择 human 或 mouse: " SPECIES
fi
if [[ -z "$TISSUE" ]]; then
  read -r -p "组织/样本描述（自然语言，将传给 LLM）: " TISSUE
fi

case "$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')" in
  human|h)
    export SPECIES=human
    ;;
  mouse|m|mus|mm)
    export SPECIES=mouse
    ;;
  *)
    echo "[run] ERROR: --species 必须是 human 或 mouse（当前: $SPECIES）" >&2
    exit 1
    ;;
esac

if [[ -z "${TISSUE// /}" ]]; then
  echo "[run] ERROR: 组织描述不能为空。" >&2
  exit 1
fi
export TISSUE

if [[ -n "$OUT_OVERRIDE" ]]; then
  export CELLTYPE_PIPELINE_OUT="$OUT_OVERRIDE"
fi
export CELLTYPE_PIPELINE_OUT="${CELLTYPE_PIPELINE_OUT:-$REPO_ROOT/outputs}"
case "$CELLTYPE_PIPELINE_OUT" in
  /*) ;;
  *) export CELLTYPE_PIPELINE_OUT="$REPO_ROOT/$CELLTYPE_PIPELINE_OUT" ;;
esac

if [[ -n "$INPUT_H5AD_ARG" ]]; then
  INPUT_H5AD="$INPUT_H5AD_ARG"
elif [[ -n "${INPUT_H5AD:-}" ]]; then
  INPUT_H5AD="${INPUT_H5AD}"
else
  INPUT_H5AD="SC.h5ad"
fi
case "$INPUT_H5AD" in
  /*) ;;
  *) INPUT_H5AD="$REPO_ROOT/$INPUT_H5AD" ;;
esac
export INPUT_H5AD

if [[ ! -f "$INPUT_H5AD" ]]; then
  echo "[run] ERROR: INPUT_H5AD not found: $INPUT_H5AD" >&2
  exit 1
fi

mkdir -p "$CELLTYPE_PIPELINE_OUT/Exchange" "$CELLTYPE_PIPELINE_OUT/figures"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "SPECIES=$SPECIES TISSUE=$TISSUE"
log "INPUT_H5AD=$INPUT_H5AD"
log "CELLTYPE_PIPELINE_OUT=$CELLTYPE_PIPELINE_OUT"

log "Step 1/5: cluster markers + working h5ad (.workflow_adata.h5ad)"
python3 "$PIPELINE_DIR/01_GetCelltypeMarker.py" --input "$INPUT_H5AD" --out "$CELLTYPE_PIPELINE_OUT"

log "Step 2/5: GO enrichment (clusterProfiler)"
Rscript "$PIPELINE_DIR/02_MarkersEnrichment.R"

log "Step 3/5: LLM coarse cell-type annotation"
python3 "$PIPELINE_DIR/03_LLMAnnotation.py" --out "$CELLTYPE_PIPELINE_OUT" --tissue "$TISSUE"

log "Step 4/5: UMAP by cell type (update working h5ad)"
python3 "$PIPELINE_DIR/04_Plot.py" --out "$CELLTYPE_PIPELINE_OUT"

log "Step 5/5: gold-standard AUCell + 最终 h5ad — 各 Leiden 金标准验证结果见下方"
python3 "$PIPELINE_DIR/05_AUCellVerification.py" --out "$CELLTYPE_PIPELINE_OUT"

log "Finished (Leiden 验证表已在上文输出)."
log "  Final h5ad: ${CELLTYPE_PIPELINE_OUT}/${OUTPUT_H5AD_NAME:-SC_annotated.h5ad}"
log "  Tables/figures: $CELLTYPE_PIPELINE_OUT/Exchange/ , $CELLTYPE_PIPELINE_OUT/figures/"
