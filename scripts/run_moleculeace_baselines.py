from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from xgboost import XGBRegressor

from pccfrl.fingerprints import MorganConfig, featurize_smiles
from pccfrl.metrics import cliff_pair_metrics, regression_metrics


DATASET_NAMES = [
    "chembl204_ki",
    "chembl214_ki",
    "chembl218_ec50",
    "chembl219_ki",
    "chembl228_ki",
    "chembl231_ki",
    "chembl233_ki",
    "chembl234_ki",
    "chembl235_ec50",
    "chembl236_ki",
    "chembl237_ec50",
    "chembl237_ki",
    "chembl238_ki",
    "chembl239_ec50",
    "chembl244_ki",
    "chembl262_ki",
    "chembl264_ki",
    "chembl287_ki",
    "chembl1862_ki",
    "chembl1871_ki",
    "chembl2034_ki",
    "chembl2047_ec50",
    "chembl2147_ki",
    "chembl2835_ki",
    "chembl2971_ki",
    "chembl3979_ec50",
    "chembl4005_ki",
    "chembl4203_ki",
    "chembl4616_ec50",
    "chembl4792_ki",
]


def parse_csv_arg(value: str) -> list[str]:
    if value.lower() in {"all", "*"}:
        return DATASET_NAMES
    return [v.strip() for v in value.split(",") if v.strip()]


def make_model(name: str, seed: int, threads: int, gpu_id: int | None = None):
    if name == "ridge":
        return Ridge(alpha=1.0, random_state=seed)
    if name == "extra_trees":
        return ExtraTreesRegressor(
            n_estimators=400,
            random_state=seed,
            n_jobs=threads,
            max_features="sqrt",
            min_samples_leaf=1,
        )
    if name == "rf":
        return RandomForestRegressor(
            n_estimators=400,
            random_state=seed,
            n_jobs=threads,
            max_features="sqrt",
            min_samples_leaf=1,
        )
    if name == "xgb_hist":
        return XGBRegressor(
            n_estimators=700,
            max_depth=6,
            learning_rate=0.03,
            subsample=0.85,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            objective="reg:squarederror",
            eval_metric="rmse",
            tree_method="hist",
            n_jobs=threads,
            random_state=seed,
        )
    if name == "xgb_gpu":
        device = f"cuda:{gpu_id}" if gpu_id is not None else "cuda"
        return XGBRegressor(
            n_estimators=900,
            max_depth=6,
            learning_rate=0.025,
            subsample=0.85,
            colsample_bytree=0.8,
            reg_lambda=1.0,
            objective="reg:squarederror",
            eval_metric="rmse",
            tree_method="hist",
            device=device,
            n_jobs=threads,
            random_state=seed,
        )
    raise ValueError(f"Unknown model: {name}")


def load_one_dataset(name: str, data_dir: str):
    from skfp.datasets.moleculeace import load_moleculeace_dataset, load_moleculeace_splits

    smiles, y = load_moleculeace_dataset(name, data_dir=data_dir)
    train_idx, test_idx = load_moleculeace_splits(
        name,
        split_type="activity_cliff",
        data_dir=data_dir,
    )
    y = np.asarray(y, dtype=np.float32)
    return list(smiles), y, np.asarray(train_idx), np.asarray(test_idx)


def run_one(
    dataset: str,
    model_name: str,
    seed: int,
    data_dir: str,
    model_threads: int,
    tanimoto_threshold: float,
    cliff_delta: float,
    smooth_delta: float,
    gpu_id: int | None,
) -> dict:
    started = time.time()
    smiles, y, train_idx, test_idx = load_one_dataset(dataset, data_dir)
    x, invalid = featurize_smiles(smiles, MorganConfig(radius=2, n_bits=2048))

    model = make_model(model_name, seed, model_threads, gpu_id)
    model.fit(x[train_idx], y[train_idx])
    pred = model.predict(x[test_idx])

    from pccfrl.metrics import CliffMetricConfig

    row = {
        "dataset": dataset,
        "model": model_name,
        "seed": seed,
        "n_total": int(len(y)),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "invalid_smiles": len(invalid),
        "seconds": round(time.time() - started, 3),
        "gpu_id": gpu_id,
    }
    row.update(regression_metrics(y[test_idx], pred))
    row.update(
        cliff_pair_metrics(
            [smiles[i] for i in test_idx],
            y[test_idx],
            pred,
            CliffMetricConfig(
                tanimoto_threshold=tanimoto_threshold,
                cliff_delta=cliff_delta,
                smooth_delta=smooth_delta,
            ),
        )
    )
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", default="all")
    parser.add_argument("--models", default="ridge,rf,xgb_hist")
    parser.add_argument("--seeds", default="0")
    parser.add_argument("--data-dir", default="data/raw/skfp_cache")
    parser.add_argument("--output", default="results/moleculeace_baselines.jsonl")
    parser.add_argument("--max-workers", type=int, default=8)
    parser.add_argument("--model-threads", type=int, default=4)
    parser.add_argument("--tanimoto-threshold", type=float, default=0.7)
    parser.add_argument("--cliff-delta", type=float, default=1.0)
    parser.add_argument("--smooth-delta", type=float, default=0.3)
    parser.add_argument("--gpu-ids", default="")
    args = parser.parse_args()

    os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

    datasets = parse_csv_arg(args.datasets)
    models = parse_csv_arg(args.models)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    raw_gpu_ids = [x.strip() for x in args.gpu_ids.split(",") if x.strip()]
    gpu_ids = [int(x) for x in raw_gpu_ids]

    jobs = [
        (dataset, model, seed, gpu_ids[i % len(gpu_ids)] if gpu_ids else None)
        for i, (dataset, model, seed) in enumerate(
            (dataset, model, seed)
            for dataset in datasets
            for model in models
            for seed in seeds
        )
    ]
    print(f"Running {len(jobs)} jobs -> {output}", flush=True)

    with output.open("a", encoding="utf-8") as fh:
        with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
            futures = {
                pool.submit(
                    run_one,
                    dataset,
                    model,
                    seed,
                    args.data_dir,
                    args.model_threads,
                    args.tanimoto_threshold,
                    args.cliff_delta,
                    args.smooth_delta,
                    gpu_id,
                ): (dataset, model, seed, gpu_id)
                for dataset, model, seed, gpu_id in jobs
            }
            for fut in as_completed(futures):
                dataset, model, seed, gpu_id = futures[fut]
                try:
                    row = fut.result()
                except Exception as exc:
                    row = {
                        "dataset": dataset,
                        "model": model,
                        "seed": seed,
                        "gpu_id": gpu_id,
                        "error": repr(exc),
                    }
                fh.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
                fh.flush()
                status = "ERROR" if "error" in row else f"RMSE={row['rmse']:.4f}"
                print(f"[{dataset} | {model} | seed={seed}] {status}", flush=True)


if __name__ == "__main__":
    main()
