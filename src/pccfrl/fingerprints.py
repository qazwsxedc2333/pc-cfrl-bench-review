from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


@dataclass(frozen=True)
class MorganConfig:
    radius: int = 2
    n_bits: int = 2048


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return mol


def morgan_bitvect(smiles: str, config: MorganConfig = MorganConfig()):
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=config.radius,
        fpSize=config.n_bits,
    )
    return generator.GetFingerprint(mol)


def morgan_numpy(smiles: str, config: MorganConfig = MorganConfig()) -> np.ndarray:
    fp = morgan_bitvect(smiles, config)
    arr = np.zeros((config.n_bits,), dtype=np.float32)
    if fp is None:
        return arr
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def featurize_smiles(
    smiles: Iterable[str],
    config: MorganConfig = MorganConfig(),
) -> tuple[np.ndarray, list[int]]:
    features: list[np.ndarray] = []
    invalid: list[int] = []
    for idx, smi in enumerate(smiles):
        fp = morgan_numpy(smi, config)
        if not fp.any():
            invalid.append(idx)
        features.append(fp)
    return np.asarray(features, dtype=np.float32), invalid


def bitvectors_for_similarity(
    smiles: Iterable[str],
    config: MorganConfig = MorganConfig(),
):
    return [morgan_bitvect(smi, config) for smi in smiles]


def tanimoto(fp_a, fp_b) -> float:
    if fp_a is None or fp_b is None:
        return 0.0
    return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))

