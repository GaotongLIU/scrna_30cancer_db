# scrna_30cancer_db

Subsampled (≤5000 cells/cell-type) single-cell reference across 30 cancer types, derived from the pan-cancer atlas:

> Systematic dissection of tumor-normal single-cell ecosystems across a thousand tumors of 30 cancer types
> Nature Communications (2024), DOI: [10.1038/s41467-024-48310-4](https://www.nature.com/articles/s41467-024-48310-4)

## What's in it

- 30 cancer types, **666,789** cells total (subsampled).
- Per cancer: `<STEM>_sub_trunc.h5ad` (AnnData) + `<STEM>_sub_trunc.h5ad.xz` (compressed twin, GitHub Release asset) + `<STEM>_inf_aver.csv` (cell2location reference) + `<STEM>_meta.csv` (cell-type × n_cells, dataset-collapsed).
- `cell_counts_all.csv`: cancer-level total cells overview.
- `manifest.csv`: full index (cancer → release URL, sizes, local paths).

## Quick use

```bash
# same host as the archive (absolute-path direct read)
python loader.py net --load

# anywhere else (auto-downloads .xz from GitHub Releases)
python loader.py SARC-EWING --load

python loader.py --list   # list all 30
```

```python
from loader import load_cancer
adata, inf_aver, meta = load_cancer('net')
```

## Cancer types

| code | organ | tissue | cells | datasets | cell types |
|---|---|---|---|---|---|
| NB | Adrenal Gland | Tumor | 11,600 | 1 | 7 |
| CHOL | Bile Duct | Tumor | 17,409 | 4 | 8 |
| BLCA | Bladder | Tumor | 10,227 | 2 | 8 |
| ALL | Bone Marrow | Tumor | 10,547 | 1 | 5 |
| CLL | Bone Marrow | Tumor | 2,600 | 1 | 4 |
| MM | Bone Marrow | Tumor | 23,307 | 1 | 5 |
| SARC-EWING | Bone | Tumor | 4,376 | 1 | 5 |
| SARC-OST | Bone | Tumor | 22,175 | 2 | 7 |
| GBM | Brain | Tumor | 25,053 | 2 | 8 |
| LGG | Brain | Tumor | 23,115 | 2 | 8 |
| BRCA | Breast | Tumor | 38,950 | 5 | 8 |
| CRC | Colon | Tumor | 38,129 | 6 | 8 |
| HNSC | Head and Neck | Tumor | 40,000 | 7 | 8 |
| RCC | Kidney | Tumor | 30,261 | 4 | 8 |
| WILM | Kidney | Tumor | 14,526 | 1 | 8 |
| HCC | Liver | Tumor | 33,434 | 5 | 8 |
| LC | Lung | Tumor | 35,648 | 6 | 8 |
| NHL | Lymph Node | Tumor | 25,335 | 4 | 8 |
| SARC-RHAB | Muscle | Tumor | 6,759 | 2 | 5 |
| NET | Neuroendocrine | Tumor | 17,193 | 2 | 8 |
| OV | Ovary | Tumor | 36,480 | 7 | 8 |
| PAAD | Pancreas | Tumor | 30,228 | 6 | 8 |
| PRAD | Prostate | Tumor | 29,109 | 4 | 8 |
| MEL | Skin | Tumor | 16,404 | 1 | 8 |
| SSCC | Skin | Tumor | 12,358 | 1 | 8 |
| STAD | Stomach | Tumor | 28,518 | 3 | 8 |
| SARC-SYN | Synovium | Tumor | 6,702 | 1 | 6 |
| THCA | Thyroid | Tumor | 36,191 | 3 | 8 |
| UCEC | Uterus | Tumor | 20,414 | 1 | 8 |
| UVM | Uvea | Tumor | 19,741 | 3 | 8 |

Compressed `.xz` total ≈ 2.55 GiB (Release assets); raw `.h5ad` kept locally only.

## Citation

```
Kang J, Lee J-H et al. Systematic dissection of tumor-normal single-cell ecosystems across a thousand tumors of 30 cancer types. Nature Communications 2024;15:4067. doi:10.1038/s41467-024-48310-4
```

## License

MIT (this package). The underlying data follows the original study's terms.

Repo: https://github.com/GaotongLIU/scrna_30cancer_db