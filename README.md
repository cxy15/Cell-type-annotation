# 自动细胞注释

单细胞 **Leiden 聚类 → 富集分析 → LLM 粗注释 → UMAP 可视化 → 金标准 AUCell 校验** 的本地流水线，脚本位于 `pipeline/`，通过根目录 `run_workflow.sh` 一键执行。

## 环境要求

因本项目依赖复杂，故不使用虚拟环境，请激活对应的Conda环境后，运行脚本
主要包：`scanpy`、`omicverse`、`openai`、`pandas`、`matplotlib` 等（见 `requirements.txt`）。

## 快速开始

### 准备数据

项目会自动搜索 **AnnData**（`.h5ad`）文件，该文件可以被放在项目内任意路径，但项目同时只能有一个单细胞数据文件。流水线假设：

- 已在**同一套高变基因**空间上计算好 **`X_umap`**（及通常的 `leiden`）；
- 若缺少 `leiden`，步骤 1 会在必要时仅构建邻接图并做 Leiden，**不会重算 UMAP**。


### 运行完整流程

在项目根目录执行（需可执行权限：`chmod +x run_workflow.sh`）：

```bash
./run_workflow.sh --species human --tissue "人肺腺癌组织"
```

- `--species`：`human` 或 `mouse`（亦接受 `h` / `m` 等简写，脚本会规范化）
- `--tissue`：组织或样本的自然语言描述，会传入 LLM 用于注释

如果没有配置好LLM API，脚本会对LLM设置进行交互式询问：

   - `LLM_BASE_URL` — API 根地址，如 `https://api.openai.com/v1`
   - `LLM_MODEL` — 模型名
   - `LLM_API_KEY` — 密钥（本地无密钥可留空）

查看帮助：

```bash
./run_workflow.sh --help
# 或直接
./pipeline/run_workflow.sh --help
```

### 输出说明

默认在 `outputs/`（或 `--out` 指定目录）下生成：

| 路径 | 说明 |
|------|------|
| `SC_annotated.h5ad` | 最终注释结果（含 AUCell 等列）；若设置 `OUTPUT_H5AD_NAME` 则改名 |
| `Exchange/` | 中间表与 JSON：`all_markers.csv`、`enrichment_summary.csv`、`llm_annotation.json`、`gold_standard_markers.json`、`aucell_verification.csv`、`progress.json` |
| `figures/umap_celltype.png` | 按注释细胞类型着色的 UMAP |

运行过程中会短暂存在 `.workflow_adata.h5ad`（工作副本），成功结束后由步骤 5 删除，仅保留上述最终 h5ad。

## 项目结构

```
Celltype-annotation/
├── README.md                 # 本说明
├── requirements.txt          # Python 依赖
├── .env.example              # 环境变量模板 → 复制为 .env
├── run_workflow.sh           # 入口：转调 pipeline/run_workflow.sh
├── pipeline/
│   ├── run_workflow.sh       # 主流程编排（加载 .env、调用各步）
│   ├── pipeline_common.py    # 路径、日志、LLM 配置、工作/最终 h5ad 路径
│   ├── 01_GetCelltypeMarker.py   # 标记基因 + 工作 h5ad
│   ├── 02_MarkersEnrichment.R   # GO 富集（人/鼠）
│   ├── 03_LLMAnnotation.py      # LLM 粗注释
│   ├── 04_Plot.py               # UMAP 图
│   └── 05_AUCellVerification.py # 金标准 AUCell + 最终 h5ad
└── outputs/                  # 默认输出目录（运行后生成，可加入 .gitignore）
    ├── Exchange/
    ├── figures/
    └── SC_annotated.h5ad
```

（`SC.h5ad` 等数据文件由用户自行放置，不随仓库提供。）

## 流程概览

1. **01**：按 HVG 子集导出 cluster marker，写入工作 `.workflow_adata.h5ad`  
2. **02**：按物种做 GO 富集，得到通路摘要  
3. **03**：将 marker + 通路发给 LLM，得到 `llm_annotation.json`  
4. **04**：映射注释并绘制 UMAP  
5. **05**：按细胞类型请求金标准基因集，做 AUCell 最终验证注释准确性，最终将写出**唯一**成品 `SC_annotated.h5ad` 并删除工作用临时 h5ad  

程序并未设计中断机制；如果程序意外中断，可直接运行脚本单独重跑，但需确保前置文件已存在于 `Exchange/` 或对应路径中。

## 引用
寻找Marker、计算AUCell均基于Omicverse：
> OmicVerse: a framework for bridging and deepening insights across bulk and single-cell sequencing
> Zeng, Z., Ma, Y., Hu, L. et al.
> Nature Communication 2024 Jul 16. doi: 10.1038/s41467-024-50194-3.


