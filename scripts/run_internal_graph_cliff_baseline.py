from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import numpy as np
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, balanced_accuracy_score, f1_score, matthews_corrcoef, roc_auc_score

import run_chemical_ood_cliff_probe as chem_ood
import run_pair_target_cliff_classifier_blend_probe as cliff
import run_pair_target_cold_target_blend_probe as blend
import run_pair_target_cold_target_probe as cold


RECORDS = None
GRAPH_X: np.ndarray | None = None
GRAPH_A: np.ndarray | None = None
GRAPH_M: np.ndarray | None = None
GRAPH_ID_A: np.ndarray | None = None
GRAPH_ID_B: np.ndarray | None = None
FEAT_DIM = 0


def parse_csv_arg(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(",") if x.strip()]


def parse_int_csv_arg(value: str) -> set[int]:
    return {int(x) for x in parse_csv_arg(value)}


def task_key(row: dict) -> tuple[str, int, int, str]:
    return (str(row["split_mode"]), int(row["seed"]), int(row["fold"]), str(row["variant"]))


def completed_task_keys(path: Path) -> set[tuple[str, int, int, str]]:
    keys = set()
    if not path.exists():
        return keys
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            try:
                keys.add(task_key(json.loads(line)))
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                continue
    return keys


def atom_features(atom: Chem.Atom) -> np.ndarray:
    atomic_options = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]
    degree_options = [0, 1, 2, 3, 4, 5]
    charge_options = [-1, 0, 1]
    feats = []
    z = atom.GetAtomicNum()
    feats.extend([1.0 if z == v else 0.0 for v in atomic_options])
    feats.append(1.0 if z not in atomic_options else 0.0)
    deg = atom.GetDegree()
    feats.extend([1.0 if deg == v else 0.0 for v in degree_options])
    charge = atom.GetFormalCharge()
    feats.extend([1.0 if charge == v else 0.0 for v in charge_options])
    feats.append(1.0 if charge not in charge_options else 0.0)
    feats.append(1.0 if atom.GetIsAromatic() else 0.0)
    feats.append(float(atom.GetTotalNumHs()) / 4.0)
    return np.asarray(feats, dtype=np.float32)


def graph_arrays_for_smiles(smiles: str, max_atoms: int, feat_dim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mol = Chem.MolFromSmiles(smiles)
    x = np.zeros((max_atoms, feat_dim), dtype=np.float32)
    a = np.zeros((max_atoms, max_atoms), dtype=np.float16)
    m = np.zeros(max_atoms, dtype=np.float32)
    if mol is None:
        return x, a, m
    n = min(mol.GetNumAtoms(), max_atoms)
    for i, atom in enumerate(mol.GetAtoms()):
        if i >= max_atoms:
            break
        x[i] = atom_features(atom)
        m[i] = 1.0
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i < max_atoms and j < max_atoms:
            a[i, j] = 1.0
            a[j, i] = 1.0
    if n:
        a[np.arange(n), np.arange(n)] = 1.0
        deg = np.asarray(a[:n, :n].sum(axis=1), dtype=np.float32)
        deg_inv = np.power(np.maximum(deg, 1.0), -0.5)
        a[:n, :n] = (deg_inv[:, None] * np.asarray(a[:n, :n], dtype=np.float32) * deg_inv[None, :]).astype(np.float16)
    return x, a, m


def build_graph_tensors(records, max_atoms_arg: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    smiles = sorted(set(records["smiles_a"]).union(set(records["smiles_b"])))
    counts = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            counts.append(mol.GetNumAtoms())
    max_atoms = max_atoms_arg if max_atoms_arg > 0 else max(counts)
    max_atoms = int(max_atoms)
    feat_dim = len(atom_features(Chem.MolFromSmiles("CC").GetAtomWithIdx(0)))
    x_all = np.zeros((len(smiles), max_atoms, feat_dim), dtype=np.float32)
    a_all = np.zeros((len(smiles), max_atoms, max_atoms), dtype=np.float16)
    m_all = np.zeros((len(smiles), max_atoms), dtype=np.float32)
    for i, smi in enumerate(smiles):
        x_all[i], a_all[i], m_all[i] = graph_arrays_for_smiles(smi, max_atoms, feat_dim)
    index = {smi: i for i, smi in enumerate(smiles)}
    id_a = records["smiles_a"].map(index).to_numpy(dtype=np.int64)
    id_b = records["smiles_b"].map(index).to_numpy(dtype=np.int64)
    return x_all, a_all, m_all, id_a, id_b


def make_split(records, split_mode: str, seed: int, fold: int, n_folds: int, val_fraction: float):
    if split_mode in {"random", "target_cluster", "target_cluster_balanced"}:
        folds = cliff.make_folds(records, seed, n_folds, split_mode)
        heldout = folds[fold]
        outer_train_mask, test_mask = cold.train_test_masks(records, heldout)
        val_targets = blend.validation_targets(records, heldout, seed, fold, val_fraction)
        val_mask = records["dataset_norm"].isin(val_targets).to_numpy()
        return outer_train_mask, val_mask, test_mask, ";".join(heldout), ";".join(val_targets)
    if split_mode in {"cold_pair", "pair_scaffold", "ligand_component"}:
        folds = chem_ood.balanced_group_folds(records, split_mode, seed, n_folds)
        test_groups = folds[fold]
        all_groups = set(records[split_mode].astype(str).unique())
        train_groups = all_groups - test_groups
        val_groups = chem_ood.validation_groups(records, split_mode, train_groups, seed, fold, val_fraction)
        outer_train_mask = records[split_mode].astype(str).isin(train_groups).to_numpy()
        test_mask = records[split_mode].astype(str).isin(test_groups).to_numpy()
        val_mask = records[split_mode].astype(str).isin(val_groups).to_numpy()
        return outer_train_mask, val_mask, test_mask, ";".join(sorted(test_groups)), ";".join(sorted(val_groups))
    raise ValueError(split_mode)


def n_folds_for(records, split_mode: str, seed: int, n_folds: int) -> int:
    if split_mode in {"random", "target_cluster", "target_cluster_balanced"}:
        return len(cliff.make_folds(records, seed, n_folds, split_mode))
    return len(chem_ood.balanced_group_folds(records, split_mode, seed, n_folds))


def row_metrics(y_delta: np.ndarray, prob: np.ndarray, threshold: float) -> dict[str, float]:
    y = (np.abs(y_delta) >= threshold).astype(np.int32)
    pred = (prob >= 0.5).astype(np.int32)
    return {
        "roc_auc": float(roc_auc_score(y, prob)) if len(np.unique(y)) > 1 else float("nan"),
        "pr_auc": float(average_precision_score(y, prob)) if len(np.unique(y)) > 1 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, pred)) if len(np.unique(pred)) > 1 else 0.0,
        "abs_delta_spearman": float(spearmanr(np.abs(y_delta), prob).statistic) if np.std(prob) > 0 else float("nan"),
        "score_mean": float(np.mean(prob)),
        "score_std": float(np.std(prob)),
    }


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    import torch

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def build_model(feat_dim: int, hidden: int, layers: int, dropout: float):
    import torch
    from torch import nn

    class GraphEncoder(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.input = nn.Linear(feat_dim, hidden)
            self.layers = nn.ModuleList([nn.Sequential(nn.Linear(hidden, hidden), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden, hidden)) for _ in range(layers)])
            self.norms = nn.ModuleList([nn.LayerNorm(hidden) for _ in range(layers)])

        def forward(self, x, adj, mask):
            h = self.input(x)
            for layer, norm in zip(self.layers, self.norms):
                msg = torch.bmm(adj, h)
                h = norm(h + layer(msg))
                h = torch.nn.functional.silu(h)
                h = h * mask.unsqueeze(-1)
            denom = mask.sum(dim=1, keepdim=True).clamp_min(1.0)
            return (h * mask.unsqueeze(-1)).sum(dim=1) / denom

    class PairGraphModel(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.encoder = GraphEncoder()
            self.head = nn.Sequential(
                nn.Linear(hidden * 3, hidden),
                nn.SiLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden, hidden // 2),
                nn.SiLU(),
                nn.Linear(hidden // 2, 1),
            )

        def forward(self, xa, aa, ma, xb, ab, mb):
            ha = self.encoder(xa, aa, ma)
            hb = self.encoder(xb, ab, mb)
            feat = torch.cat([torch.abs(ha - hb), ha * hb, (ha + hb) * 0.5], dim=1)
            return self.head(feat).squeeze(1)

    return PairGraphModel()


def batch_graph(ids: np.ndarray, device):
    assert GRAPH_X is not None and GRAPH_A is not None and GRAPH_M is not None
    import torch

    x = torch.from_numpy(GRAPH_X[ids]).to(device)
    a = torch.from_numpy(np.asarray(GRAPH_A[ids], dtype=np.float32)).to(device)
    m = torch.from_numpy(GRAPH_M[ids]).to(device)
    return x, a, m


def train_predict(task: dict, train_mask: np.ndarray, val_mask: np.ndarray, test_mask: np.ndarray) -> np.ndarray:
    assert RECORDS is not None and GRAPH_ID_A is not None and GRAPH_ID_B is not None
    import torch

    set_seed(int(task["seed"]))
    device = torch.device(f"cuda:{task['gpu_id']}" if task["gpu_id"] is not None and torch.cuda.is_available() else "cpu")
    model = build_model(FEAT_DIM, int(task["hidden"]), int(task["layers"]), float(task["dropout"])).to(device)
    y_all = (np.abs(RECORDS["delta"].to_numpy(dtype=np.float32)) >= task["cliff_threshold"]).astype(np.float32)
    y_train = y_all[train_mask & ~val_mask]
    pos = float(y_train.sum())
    neg = float(len(y_train) - pos)
    loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor([neg / max(pos, 1.0)], device=device))
    opt = torch.optim.AdamW(model.parameters(), lr=task["lr"], weight_decay=task["weight_decay"])
    y_t = torch.from_numpy(y_all)
    train_idx = np.flatnonzero(train_mask & ~val_mask)
    val_idx = np.flatnonzero(val_mask)
    best_state = None
    best_val = float("inf")
    stale = 0
    rng = np.random.default_rng(71_003 + int(task["seed"]) * 1009 + int(task["fold"]))
    batch_size = int(task["batch_size"])
    for _epoch in range(int(task["epochs"])):
        rng.shuffle(train_idx)
        model.train()
        for start in range(0, len(train_idx), batch_size):
            idx = train_idx[start : start + batch_size]
            xa, aa, ma = batch_graph(GRAPH_ID_A[idx], device)
            xb, ab, mb = batch_graph(GRAPH_ID_B[idx], device)
            logits = model(xa, aa, ma, xb, ab, mb)
            loss = loss_fn(logits, y_t[idx].to(device))
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
        if len(val_idx):
            losses = []
            model.eval()
            with torch.no_grad():
                for start in range(0, len(val_idx), batch_size):
                    idx = val_idx[start : start + batch_size]
                    xa, aa, ma = batch_graph(GRAPH_ID_A[idx], device)
                    xb, ab, mb = batch_graph(GRAPH_ID_B[idx], device)
                    logits = model(xa, aa, ma, xb, ab, mb)
                    losses.append(float(loss_fn(logits, y_t[idx].to(device)).detach().cpu()))
            val_loss = float(np.mean(losses))
            if val_loss < best_val:
                best_val = val_loss
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                stale = 0
            else:
                stale += 1
            if stale >= int(task["patience"]):
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    probs = []
    test_idx = np.flatnonzero(test_mask)
    model.eval()
    with torch.no_grad():
        for start in range(0, len(test_idx), batch_size):
            idx = test_idx[start : start + batch_size]
            xa, aa, ma = batch_graph(GRAPH_ID_A[idx], device)
            xb, ab, mb = batch_graph(GRAPH_ID_B[idx], device)
            logits = model(xa, aa, ma, xb, ab, mb)
            probs.append(torch.sigmoid(logits).detach().cpu().numpy())
    return np.concatenate(probs).astype(np.float32)


def run_task(task: dict) -> dict:
    assert RECORDS is not None
    train_mask, val_mask, test_mask, test_groups, val_groups = make_split(
        RECORDS, task["split_mode"], int(task["seed"]), int(task["fold"]), int(task["folds"]), float(task["val_fraction"])
    )
    prob = train_predict(task, train_mask, val_mask, test_mask)
    y_test_delta = RECORDS.loc[test_mask, "delta"].to_numpy(dtype=np.float32)
    row = {
        "benchmark": "PC-CFRL-internal",
        "split_mode": task["split_mode"],
        "variant": "graph_gcn_pair",
        "seed": int(task["seed"]),
        "fold": int(task["fold"]),
        "heldout_groups": test_groups,
        "val_groups": val_groups,
        "n_train": int(train_mask.sum()),
        "n_val": int(val_mask.sum()),
        "n_test": int(test_mask.sum()),
        "test_positive_rate": float((np.abs(y_test_delta) >= task["cliff_threshold"]).mean()),
        "train_positive_rate": float((np.abs(RECORDS.loc[train_mask, "delta"].to_numpy(dtype=np.float32)) >= task["cliff_threshold"]).mean()),
        "cliff_threshold": float(task["cliff_threshold"]),
        "gpu_id": task["gpu_id"],
    }
    if "pair_id" in RECORDS:
        train_pairs = set(RECORDS.loc[train_mask, "pair_id"].astype(int))
        row["pair_overlap_rate"] = float(RECORDS.loc[test_mask, "pair_id"].astype(int).isin(train_pairs).mean())
    if "canon_a" in RECORDS:
        train_ligands = set(RECORDS.loc[train_mask, "canon_a"]).union(set(RECORDS.loc[train_mask, "canon_b"]))
        row["any_ligand_overlap_rate"] = float((RECORDS.loc[test_mask, "canon_a"].isin(train_ligands) | RECORDS.loc[test_mask, "canon_b"].isin(train_ligands)).mean())
    row.update(row_metrics(y_test_delta, prob, float(task["cliff_threshold"])))
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", default="results/cross_target_analog_ambiguity_t07_pairs.csv")
    parser.add_argument("--saprot-features", default="data/processed/saprot_target_features.csv")
    parser.add_argument("--saprot-pca", type=int, default=32)
    parser.add_argument("--split-modes", default="target_cluster_balanced,ligand_component")
    parser.add_argument("--seeds", default="0,1,2,3,4")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--fold-filter", default="")
    parser.add_argument("--val-fraction", type=float, default=0.15)
    parser.add_argument("--cliff-threshold", type=float, default=1.0)
    parser.add_argument("--max-atoms", type=int, default=0)
    parser.add_argument("--hidden", type=int, default=64)
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.15)
    parser.add_argument("--epochs", type=int, default=8)
    parser.add_argument("--patience", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--gpu-id", type=int, default=-1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", default="results/internal_graph_cliff_baseline_seed0_4.jsonl")
    args = parser.parse_args()

    records = chem_ood.add_ood_columns(cold.make_records(args.pairs))
    cold.RECORDS = records
    cold.SAPROT_MAP, cold.SAPROT_DIM = cold.load_saprot_features(args.saprot_features, args.saprot_pca)
    global RECORDS, GRAPH_X, GRAPH_A, GRAPH_M, GRAPH_ID_A, GRAPH_ID_B, FEAT_DIM
    RECORDS = records
    GRAPH_X, GRAPH_A, GRAPH_M, GRAPH_ID_A, GRAPH_ID_B = build_graph_tensors(records, args.max_atoms)
    FEAT_DIM = int(GRAPH_X.shape[-1])

    split_modes = parse_csv_arg(args.split_modes)
    seeds = [int(x) for x in parse_csv_arg(args.seeds)]
    fold_filter = parse_int_csv_arg(args.fold_filter) if args.fold_filter else None
    tasks = []
    for split_mode in split_modes:
        for seed in seeds:
            n = n_folds_for(records, split_mode, seed, args.folds)
            for fold in range(n):
                if fold_filter is not None and fold not in fold_filter:
                    continue
                tasks.append(
                    {
                        "split_mode": split_mode,
                        "seed": seed,
                        "fold": fold,
                        "folds": args.folds,
                        "val_fraction": args.val_fraction,
                        "cliff_threshold": args.cliff_threshold,
                        "hidden": args.hidden,
                        "layers": args.layers,
                        "dropout": args.dropout,
                        "epochs": args.epochs,
                        "patience": args.patience,
                        "batch_size": args.batch_size,
                        "lr": args.lr,
                        "weight_decay": args.weight_decay,
                        "gpu_id": args.gpu_id if args.gpu_id >= 0 else None,
                        "variant": "graph_gcn_pair",
                    }
                )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if args.resume:
        done = completed_task_keys(out)
        tasks = [task for task in tasks if task_key(task) not in done]
    print(
        f"records={len(records)} targets={records['dataset_norm'].nunique()} unique_graphs={len(GRAPH_X)} "
        f"max_atoms={GRAPH_X.shape[1]} feat_dim={FEAT_DIM} tasks={len(tasks)}",
        flush=True,
    )
    with out.open("a", encoding="utf-8") as fh:
        for task in tasks:
            row = run_task(task)
            fh.write(json.dumps(row, sort_keys=True) + "\n")
            fh.flush()
            print(
                f"[{row['split_mode']} seed={row['seed']} fold={row['fold']} graph_gcn_pair] "
                f"roc={row['roc_auc']:.4f} pr={row['pr_auc']:.4f}",
                flush=True,
            )


if __name__ == "__main__":
    main()
