#!/usr/bin/env python
"""loader.py — scRNA 30-cancer db 快速调用入口.

load_cancer(name) -> (adata, inf_aver_df, meta_df)
  - 同服务器（本地归档可达）: 直读绝对路径 <LOCAL_BASE>/<STEM>/<STEM>_sub_trunc.h5ad
  - 否则: 从 GitHub Release 下载 <STEM>_sub_trunc.h5ad.xz, 解压到缓存, 再读
  - inf_aver.csv / meta.csv: 同服务器读本地; 否则从 repo raw URL 拉取

CLI:
  python loader.py <cancer>            # 打印 shape + cell counts, 不落 adata
  python loader.py <cancer> --load     # 真正 load 并返回 adata
  python loader.py --list               # 列出 30 个癌种
"""
import argparse
import lzma
import os
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd

DEFAULT_LOCAL_BASE = (
    "/storeData/USER/data/02.Bioinformatics_for_STOmics/01.user/"
    "liugaotong/03.prj/11.public_h5ad/03.final"
)
REPO_OWNER = "GaotongLIU"
REPO_NAME = "scrna_30cancer_db"
DEFAULT_CACHE = "~/.cache/scrna_30cancer_db"


def _local_base():
    return Path(os.environ.get("SCRNA_DB_LOCAL_BASE", DEFAULT_LOCAL_BASE))


def _cache_dir():
    return Path(os.environ.get("SCRNA_DB_CACHE", DEFAULT_CACHE)).expanduser()


def _manifest_path():
    """manifest.csv 与本脚本同目录（03.final 或 clone 的 repo 根）。"""
    here = Path(__file__).resolve().parent
    return here / "manifest.csv"


def load_manifest():
    p = _manifest_path()
    if not p.is_file():
        raise FileNotFoundError(
            f"找不到 manifest.csv: {p}（设 SCRNA_DB_MANIFEST 或在脚本同目录放 manifest.csv）"
        )
    return pd.read_csv(p)


def resolve_cancer(name: str, manifest: pd.DataFrame):
    """输入 net/NET/sarc-ewing/SARC_EWING -> (code, stem, tag, row)。
    code 用连字符(SARC-EWING); stem 用下划线(SARC_EWING) 作文件名; tag 小写(sarc-ewing)。"""
    key = name.strip().upper().replace("_", "-")
    row = manifest[manifest["cancer_code"].str.upper() == key]
    if row.empty:
        raise ValueError(
            f"未知癌种 {name!r}; 可用: python loader.py --list"
        )
    row = row.iloc[0]
    code = row["cancer_code"]
    stem = code.replace("-", "_")
    tag = code.lower()
    return code, stem, tag, row


def _download(url: str, dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[fetch] {url} -> {dest}", file=sys.stderr)
    with urllib.request.urlopen(url) as r, open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh, 16 * 1024 * 1024)


def _resolve_h5ad_path(stem: str, tag: str):
    """本地可达 -> 直读 .h5ad; 否则下载 .xz 解压到缓存再读。返回 .h5ad 路径。"""
    base = _local_base()
    local_h5ad = base / stem / f"{stem}_sub_trunc.h5ad"
    if local_h5ad.is_file():
        return local_h5ad, "local"

    cache = _cache_dir()
    cached_h5ad = cache / f"{stem}_sub_trunc.h5ad"
    if cached_h5ad.is_file():
        return cached_h5ad, "cache"

    url = (
        f"https://github.com/{REPO_OWNER}/{REPO_NAME}/releases/download/"
        f"{tag}/{stem}_sub_trunc.h5ad.xz"
    )
    xz_tmp = cache / f"{stem}_sub_trunc.h5ad.xz"
    _download(url, xz_tmp)
    cached_h5ad.parent.mkdir(parents=True, exist_ok=True)
    with lzma.open(xz_tmp, "rb") as src, open(cached_h5ad, "wb") as dst:
        shutil.copyfileobj(src, dst, 16 * 1024 * 1024)
    xz_tmp.unlink()
    return cached_h5ad, "download"


def _resolve_csv(stem: str, tag: str, which: str):
    """which = 'inf_aver' 或 'meta'。本地优先, 否则 raw URL。返回 DataFrame。"""
    base = _local_base()
    local = base / stem / f"{stem}_{which}.csv"
    if local.is_file():
        return pd.read_csv(local), "local"
    url = (
        f"https://raw.githubusercontent.com/{REPO_OWNER}/{REPO_NAME}/main/"
        f"{stem}/{stem}_{which}.csv"
    )
    tmp = _cache_dir() / f"{stem}_{which}.csv"
    _download(url, tmp)
    return pd.read_csv(tmp), "download"


def load_cancer(name: str):
    """返回 (adata, inf_aver_df, meta_df)。adata 惰性加载 scanpy。"""
    manifest = load_manifest()
    code, stem, tag, _ = resolve_cancer(name, manifest)
    import scanpy as sc

    h5ad_path, src = _resolve_h5ad_path(stem, tag)
    print(f"[h5ad ] {h5ad_path} ({src})", file=sys.stderr)
    adata = sc.read_h5ad(h5ad_path)
    inf_df, s1 = _resolve_csv(stem, tag, "inf_aver")
    meta_df, s2 = _resolve_csv(stem, tag, "meta")
    print(f"[infer] inf_aver ({s1}) {inf_df.shape}", file=sys.stderr)
    print(f"[meta ] meta ({s2}) {meta_df.shape}", file=sys.stderr)
    return adata, inf_df, meta_df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cancer", nargs="?", help="癌种 code, 如 net / SARC-EWING")
    ap.add_argument("--load", action="store_true", help="真正 load adata (默认只打印摘要)")
    ap.add_argument("--list", action="store_true", help="列出全部癌种")
    args = ap.parse_args()

    if args.list:
        m = load_manifest()
        cols = ["cancer_code", "organ_origin", "tissue", "cancer_type",
                "n_cells_total", "n_datasets", "n_celltypes"]
        print(m[cols].to_string(index=False))
        return

    if not args.cancer:
        ap.error("需要癌种 code, 或用 --list")

    manifest = load_manifest()
    code, stem, tag, row = resolve_cancer(args.cancer, manifest)
    print(f"[cancer] {code}  organ={row['organ_origin']}  "
          f"tissue={row['tissue']}  tag={tag}")

    if not args.load:
        print(f"[meta ] n_cells_total={row['n_cells_total']}  "
              f"n_datasets={row['n_datasets']}  n_celltypes={row['n_celltypes']}")
        print("(加 --load 才真正读 h5ad)")
        return

    adata, inf_df, meta_df = load_cancer(args.cancer)
    print(f"[adata] shape={adata.shape}  obs_cols={list(adata.obs.columns)}")
    print(f"[infer] shape={inf_df.shape}")
    print(meta_df.to_string(index=False))


if __name__ == "__main__":
    main()
