import json
import numpy as np
import pytest
import torch
from pathlib import Path

from pipeline.cnn.experiments import ExperimentSpec, run_experiment, _save_history, build_experiments
from pipeline.cnn.model import StandardCNN
from pipeline.cnn.training import TrainingConfig


def _make_fake_dataset(tmp_path: Path, k_mer: int, seq_size: str) -> None:
    img_size = 2 ** k_mer
    for label in (0, 1):
        label_dir = tmp_path / "fcgr_images" / str(k_mer) / seq_size / "diff" / str(label)
        label_dir.mkdir(parents=True)
        for i in range(6):
            np.save(label_dir / f"NC_test_{i}.npy", np.random.rand(img_size, img_size).astype(np.float32))


def test_experiment_spec_has_k_mer_field():
    spec = ExperimentSpec(
        name="test",
        model_cls=StandardCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2},
        config=TrainingConfig(epochs=1, batch_size=4),
        seq_size="500",
        undersample=False,
        k_mer=6,
    )
    assert spec.k_mer == 6


def test_run_experiment_returns_metrics_and_history(tmp_path):
    _make_fake_dataset(tmp_path, k_mer=6, seq_size="500")
    spec = ExperimentSpec(
        name="test_exp",
        model_cls=StandardCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2},
        config=TrainingConfig(epochs=2, batch_size=4),
        seq_size="500",
        undersample=False,
        k_mer=6,
    )
    result = run_experiment(spec, tmp_path, torch.device("cpu"))
    assert isinstance(result, tuple) and len(result) == 2
    metrics, history = result
    assert hasattr(metrics, "accuracy")
    assert len(history) == 2
    assert set(history[0].keys()) == {"epoch", "train_loss", "val_loss", "val_acc"}


def test_run_experiment_uses_spec_k_mer(tmp_path):
    # Only k_mer=7 data exists; spec.k_mer=7 should load it successfully
    _make_fake_dataset(tmp_path, k_mer=7, seq_size="500")
    spec = ExperimentSpec(
        name="test_kmer7",
        model_cls=StandardCNN,
        model_kwargs={"k_mer": 7, "num_classes": 2},
        config=TrainingConfig(epochs=1, batch_size=4),
        seq_size="500",
        undersample=False,
        k_mer=7,
    )
    metrics, history = run_experiment(spec, tmp_path, torch.device("cpu"))
    assert len(history) == 1


def test_save_history_writes_json(tmp_path):
    history = [
        {"epoch": 1, "train_loss": 0.8, "val_loss": 0.9, "val_acc": 0.6},
        {"epoch": 2, "train_loss": 0.5, "val_loss": 0.6, "val_acc": 0.7},
    ]
    _save_history("my_experiment", history, tmp_path)
    out = json.loads((tmp_path / "my_experiment.json").read_text())
    assert out["name"] == "my_experiment"
    assert len(out["epochs"]) == 2
    assert out["epochs"][0] == {"epoch": 1, "train_loss": 0.8, "val_loss": 0.9, "val_acc": 0.6}


def test_save_history_creates_dir(tmp_path):
    histories_dir = tmp_path / "results" / "histories"
    _save_history("exp", [{"epoch": 1, "train_loss": 0.5, "val_loss": 0.5, "val_acc": 0.5}], histories_dir)
    assert histories_dir.exists()
    assert (histories_dir / "exp.json").exists()


def test_build_experiments_contains_kmer_comparison():
    specs = build_experiments()
    kmer_specs = [s for s in specs if "_kmer" in s.name]
    assert len(kmer_specs) == 36  # 3 models × 3 k_mers × 4 seq_sizes


def test_kmer_specs_have_correct_k_mer_field():
    specs = build_experiments()
    kmer_specs = {s.name: s for s in specs if "_kmer" in s.name}

    spec = kmer_specs["ResNetCNN_kmer7_seq1000_undersample"]
    assert spec.k_mer == 7
    assert spec.seq_size == "1000"
    assert spec.model_kwargs["k_mer"] == 7

    spec8 = kmer_specs["TwoDRPDWT_kmer8_seq2000_undersample"]
    assert spec8.k_mer == 8
    assert spec8.seq_size == "2000"


def test_kmer_specs_dont_clash_with_original_names():
    specs = build_experiments()
    names = [s.name for s in specs]
    # Original seq_size sweep: no _kmer token
    assert "ResNetCNN_seq1000_undersample" in names
    # k_mer comparison: has _kmer token — different names
    assert "ResNetCNN_kmer6_seq1000_undersample" in names
