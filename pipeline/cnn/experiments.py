"""Experiment runner: fits multiple model/config combos and writes RESULTS.md."""

import json
import time
from dataclasses import dataclass
from pathlib import Path

import torch

from pipeline.config import get_settings
from pipeline.cnn.dataset import FCGRDiffDataset
from pipeline.cnn.model import TwoDSquaredRPDWTCNN, StandardCNN, ResNetCNN
from pipeline.cnn.training import TrainingConfig, EvalMetrics, fit, evaluate


@dataclass
class ExperimentSpec:
    name: str
    model_cls: type
    model_kwargs: dict
    config: TrainingConfig
    seq_size: str | None  # None = all sizes combined
    undersample: bool
    k_mer: int = 6


def _save_history(name: str, history: list[dict], histories_dir: Path) -> None:
    histories_dir.mkdir(parents=True, exist_ok=True)
    payload = {"name": name, "epochs": history}
    (histories_dir / f"{name}.json").write_text(json.dumps(payload, indent=2))


def _parse_results_md(path: Path) -> dict[str, "EvalMetrics"]:
    """Load metrics from an existing RESULTS.md to enable skipping completed runs."""
    if not path.exists():
        return {}
    text = path.read_text()
    metrics: dict[str, dict] = {}

    in_table = False
    for line in text.splitlines():
        if "| Experiment |" in line:
            in_table = True
            continue
        if in_table and line.startswith("|---"):
            continue
        if in_table and line.startswith("| "):
            parts = [p.strip() for p in line.split("|")[1:-1]]
            if len(parts) >= 7:
                name = parts[0]
                try:
                    metrics[name] = {
                        "accuracy": float(parts[3]),
                        "precision": float(parts[4]),
                        "recall": float(parts[5]),
                        "f1": float(parts[6]),
                        "tp": 0, "tn": 0, "fp": 0, "fn": 0,
                    }
                except (ValueError, IndexError):
                    pass
        elif in_table and not line.startswith("|"):
            in_table = False

    current_name: str | None = None
    for line in text.splitlines():
        if line.startswith("### "):
            current_name = line[4:].strip()
        elif current_name and current_name in metrics:
            if "| actual 0 |" in line:
                parts = [p.strip() for p in line.split("|")[1:-1]]
                if len(parts) >= 3:
                    try:
                        metrics[current_name]["tn"] = int(parts[1])
                        metrics[current_name]["fp"] = int(parts[2])
                    except (ValueError, IndexError):
                        pass
            elif "| actual 1 |" in line:
                parts = [p.strip() for p in line.split("|")[1:-1]]
                if len(parts) >= 3:
                    try:
                        metrics[current_name]["fn"] = int(parts[1])
                        metrics[current_name]["tp"] = int(parts[2])
                    except (ValueError, IndexError):
                        pass

    return {name: EvalMetrics(**m) for name, m in metrics.items()}


def run_experiment(spec: ExperimentSpec, data_dir: Path, device: torch.device) -> tuple[EvalMetrics, list[dict]]:
    print(f"\n{'='*60}")
    print(f"EXPERIMENT: {spec.name}")
    print(f"{'='*60}")

    if "k_mer" in spec.model_kwargs and spec.model_kwargs["k_mer"] != spec.k_mer:
        raise ValueError(
            f"spec.k_mer={spec.k_mer} conflicts with model_kwargs['k_mer']={spec.model_kwargs['k_mer']}"
        )

    ds = FCGRDiffDataset(
        k_mer=spec.k_mer,
        seq_size=spec.seq_size,
        data_dir=data_dir,
        undersample=spec.undersample,
    )
    label_counts = {0: 0, 1: 0}
    for _, label in ds.samples:
        label_counts[label] += 1
    print(f"Dataset: {len(ds)} samples | class 0: {label_counts[0]} | class 1: {label_counts[1]}")

    model = spec.model_cls(**spec.model_kwargs)
    param_count = sum(p.numel() for p in model.parameters())
    print(f"Model: {spec.model_cls.__name__} | params: {param_count:,}")

    t0 = time.time()
    model, history = fit(model, ds, config=spec.config, device=device, verbose=True)
    elapsed = time.time() - t0

    metrics = evaluate(model, ds, config=spec.config, device=device, verbose=True)
    print(f"Training time: {elapsed:.1f}s")
    return metrics, history


def build_experiments() -> list[ExperimentSpec]:
    experiments = []

    for seq_size in ["500", "1000", "1500", "2000"]:
        experiments.append(ExperimentSpec(
            name=f"StandardCNN_seq{seq_size}_undersample",
            model_cls=StandardCNN,
            model_kwargs={"k_mer": 6, "num_classes": 2, "dropout": 0.5},
            config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
            seq_size=seq_size,
            undersample=True,
        ))

    experiments.append(ExperimentSpec(
        name="StandardCNN_allseq_undersample",
        model_cls=StandardCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2, "dropout": 0.5},
        config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=256),
        seq_size=None,
        undersample=True,
    ))

    for seq_size in ["500", "1000", "1500", "2000"]:
        experiments.append(ExperimentSpec(
            name=f"ResNetCNN_seq{seq_size}_undersample",
            model_cls=ResNetCNN,
            model_kwargs={"k_mer": 6, "num_classes": 2, "dropout": 0.4},
            config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
            seq_size=seq_size,
            undersample=True,
        ))

    experiments.append(ExperimentSpec(
        name="ResNetCNN_allseq_undersample",
        model_cls=ResNetCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2, "dropout": 0.4},
        config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=256),
        seq_size=None,
        undersample=True,
    ))

    experiments.append(ExperimentSpec(
        name="TwoDRPDWT_seq1000_undersample",
        model_cls=TwoDSquaredRPDWTCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2},
        config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
        seq_size="1000",
        undersample=True,
    ))

    experiments.append(ExperimentSpec(
        name="TwoDRPDWT_allseq_undersample",
        model_cls=TwoDSquaredRPDWTCNN,
        model_kwargs={"k_mer": 6, "num_classes": 2},
        config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=256),
        seq_size=None,
        undersample=True,
    ))

    # k_mer comparison: 3 models × k_mer ∈ {6, 7, 8} × seq_size ∈ {500, 1000, 1500, 2000}
    for k_mer in [6, 7, 8]:
        for seq_size in ["500", "1000", "1500", "2000"]:
            experiments.append(ExperimentSpec(
                name=f"ResNetCNN_kmer{k_mer}_seq{seq_size}_undersample",
                model_cls=ResNetCNN,
                model_kwargs={"k_mer": k_mer, "num_classes": 2, "dropout": 0.4},
                config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
                seq_size=seq_size,
                undersample=True,
                k_mer=k_mer,
            ))
            experiments.append(ExperimentSpec(
                name=f"StandardCNN_kmer{k_mer}_seq{seq_size}_undersample",
                model_cls=StandardCNN,
                model_kwargs={"k_mer": k_mer, "num_classes": 2, "dropout": 0.5},
                config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
                seq_size=seq_size,
                undersample=True,
                k_mer=k_mer,
            ))
            experiments.append(ExperimentSpec(
                name=f"TwoDRPDWT_kmer{k_mer}_seq{seq_size}_undersample",
                model_cls=TwoDSquaredRPDWTCNN,
                model_kwargs={"k_mer": k_mer, "num_classes": 2},
                config=TrainingConfig(epochs=40, lr=1e-3, optimizer="adam", batch_size=128),
                seq_size=seq_size,
                undersample=True,
                k_mer=k_mer,
            ))

    return experiments


def write_results(results: list[tuple[ExperimentSpec, EvalMetrics]], output_path: Path) -> None:
    lines = [
        "# CNN Experiment Results",
        "",
        "**Goal:** Accuracy > 80% on FCGR diff images for variant pathogenicity classification.",
        "",
        "**Setup:** k-mer=6, images 64×64. Class imbalance handled by undersampling class 0 to match class 1 (~7 883 each). "
        "Adam optimizer with cosine annealing LR scheduler. Class-weighted CrossEntropy loss.",
        "",
        "## Results",
        "",
        "| Experiment | Seq Size | Samples | Accuracy | Precision | Recall | F1 | >80%? |",
        "|---|---|---|---|---|---|---|---|",
    ]

    for spec, m in results:
        seq_label = spec.seq_size if spec.seq_size else "all"
        samples = "balanced"
        above = "✅" if m.accuracy >= 0.80 else "❌"
        lines.append(
            f"| {spec.name} | {seq_label} | {samples} | "
            f"{m.accuracy:.4f} | {m.precision:.4f} | {m.recall:.4f} | {m.f1:.4f} | {above} |"
        )

    lines += [
        "",
        "## Confusion Matrices",
        "",
    ]

    for spec, m in results:
        lines += [
            f"### {spec.name}",
            "",
            f"|  | pred 0 | pred 1 |",
            f"|---|---|---|",
            f"| actual 0 | {m.tn} | {m.fp} |",
            f"| actual 1 | {m.fn} | {m.tp} |",
            "",
        ]

    lines += [
        "## Architecture Notes",
        "",
        "- **StandardCNN**: 3× (Conv→BN→ReLU→MaxPool) + FC(256)+Dropout+FC(2). Simple baseline with regularization.",
        "- **ResNetCNN**: Stem + 3 stages with residual blocks + AdaptiveAvgPool. Better gradient flow.",
        "- **TwoDSquaredRPDWTCNN**: Paper architecture (Alzu'bi et al. CMC 2022) with Haar DWT + random projection downsampling.",
        "",
        "## Observations",
        "",
        "- All seq_sizes (500/1000/1500/2000) contain identical class counts (data is replicated across sizes).",
        "- True signal: whether FCGR k-mer frequency differences can separate benign from pathogenic variants.",
        "- Undersampling balances to ~7 883 per class (~15 766 total).",
    ]

    output_path.write_text("\n".join(lines) + "\n")
    print(f"\nResults written to {output_path}")


def _append_results(new_results: list[tuple[ExperimentSpec, EvalMetrics]], output_path: Path) -> None:
    """Append new experiment rows to existing RESULTS.md (preserving custom content)."""
    if not output_path.exists():
        write_results(new_results, output_path)
        return

    rows = ["", f"## k-mer=8 Results ({len(new_results)} new experiments)", ""]
    rows += [
        "| Experiment | Seq Size | Samples | Accuracy | Precision | Recall | F1 | >80%? |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for spec, m in new_results:
        seq_label = spec.seq_size if spec.seq_size else "all"
        above = "✅" if m.accuracy >= 0.80 else "❌"
        rows.append(
            f"| {spec.name} | {seq_label} | balanced | "
            f"{m.accuracy:.4f} | {m.precision:.4f} | {m.recall:.4f} | {m.f1:.4f} | {above} |"
        )

    rows += ["", "### Confusion Matrices", ""]
    for spec, m in new_results:
        rows += [
            f"#### {spec.name}", "",
            "|  | pred 0 | pred 1 |", "|---|---|---|",
            f"| actual 0 | {m.tn} | {m.fp} |",
            f"| actual 1 | {m.fn} | {m.tp} |", "",
        ]

    existing = output_path.read_text().rstrip()
    output_path.write_text(existing + "\n" + "\n".join(rows) + "\n")
    print(f"\nAppended {len(new_results)} new results to {output_path}")


def main() -> None:
    settings = get_settings()
    device = torch.device(
        "cuda" if torch.cuda.is_available()
        else "mps" if torch.backends.mps.is_available()
        else "cpu"
    )
    print(f"Device: {device}")

    experiments = build_experiments()
    histories_dir = Path("results/histories")
    results_md = Path("RESULTS.md")

    new_results: list[tuple[ExperimentSpec, EvalMetrics]] = []

    for spec in experiments:
        json_path = histories_dir / f"{spec.name}.json"
        if json_path.exists():
            print(f"[SKIP] {spec.name}")
            continue

        try:
            metrics, history = run_experiment(spec, settings.data_dir, device)
            new_results.append((spec, metrics))
            try:
                _save_history(spec.name, history, histories_dir)
            except Exception as save_err:
                print(f"[WARN] Failed to save history for {spec.name}: {save_err}")
            if metrics.accuracy >= 0.80:
                print(f"\n*** TARGET REACHED: {spec.name} achieved {metrics.accuracy:.4f} accuracy ***")
        except Exception as e:
            print(f"[ERROR] {spec.name} failed: {e}")

    if new_results:
        _append_results(new_results, results_md)


if __name__ == "__main__":
    main()
