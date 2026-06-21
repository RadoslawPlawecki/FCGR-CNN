import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

MODEL_COLORS = {
    "ResNetCNN": "#2196F3",
    "StandardCNN": "#FF9800",
    "TwoDRPDWT": "#4CAF50",
}


def _load_histories(histories_dir: Path) -> list[dict]:
    if not histories_dir.exists():
        raise FileNotFoundError(f"Histories directory not found: {histories_dir}")
    results = []
    for path in sorted(histories_dir.glob("*.json")):
        try:
            data = json.loads(path.read_text())
            if "name" not in data or "epochs" not in data:
                print(f"[WARN] Skipping {path.name}: missing 'name' or 'epochs'")
                continue
            results.append(data)
        except (json.JSONDecodeError, KeyError) as exc:
            print(f"[WARN] Skipping {path.name}: {exc}")
    return results


def _seq_size_from_name(name: str) -> int:
    for part in name.split("_"):
        if part.startswith("seq"):
            try:
                return int(part[3:])
            except ValueError:
                pass
    return 0


def plot_training_curves(histories_dir: Path, figures_dir: Path) -> None:
    """One figure per model: 2 rows × N_seqsizes columns (loss curves + val_acc)."""
    all_histories = _load_histories(histories_dir)

    # Original seq_size sweep only — exclude k_mer comparison runs
    sweep = [h for h in all_histories if "_kmer" not in h["name"]]

    # For models with no non-kmer runs, fall back to kmer6 as the seq-size sweep
    represented = {h["name"].split("_")[0] for h in sweep}
    kmer6_fallback = [
        h for h in all_histories
        if "_kmer6_" in h["name"] and h["name"].split("_")[0] not in represented
    ]
    sweep += kmer6_fallback

    by_model: dict[str, list[dict]] = defaultdict(list)
    for h in sweep:
        model = h["name"].split("_")[0]
        by_model[model].append(h)

    figures_dir.mkdir(parents=True, exist_ok=True)

    for model_name, histories in by_model.items():
        histories = sorted(histories, key=lambda h: _seq_size_from_name(h["name"]))
        n_cols = len(histories)
        color = MODEL_COLORS.get(model_name, "#333333")

        with plt.style.context("seaborn-v0_8-whitegrid"):
            fig, axes = plt.subplots(2, n_cols, figsize=(4 * n_cols, 6), squeeze=False)
            fig.suptitle(model_name, fontsize=13, fontweight="bold")

            for col, h in enumerate(histories):
                epochs = [e["epoch"] for e in h["epochs"]]
                train_loss = [e["train_loss"] for e in h["epochs"]]
                val_loss = [e["val_loss"] for e in h["epochs"]]
                val_acc = [e["val_acc"] for e in h["epochs"]]

                seq_size = _seq_size_from_name(h["name"])
                col_title = f"seq_size = {seq_size}" if seq_size else h["name"]

                ax0 = axes[0][col]
                ax0.plot(epochs, train_loss, color=color, linestyle="-", label="train loss")
                ax0.plot(epochs, val_loss, color=color, linestyle="--", alpha=0.7, label="val loss")
                ax0.set_title(col_title, fontsize=11)
                ax0.set_xlabel("Epoch")
                ax0.set_ylabel("Loss")
                if col == 0:
                    ax0.legend(fontsize=9)

                ax1 = axes[1][col]
                ax1.plot(epochs, val_acc, color=color, linestyle="-")
                ax1.set_xlabel("Epoch")
                ax1.set_ylabel("Validation Accuracy")
                ax1.set_ylim(0, 1)

            plt.tight_layout()
            stem = figures_dir / f"curves_{model_name}"
            fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
            fig.savefig(f"{stem}.pdf", dpi=300, bbox_inches="tight")
            plt.close(fig)

    print(f"[INFO] Training curve figures saved to {figures_dir}")


def _kmer_from_name(name: str) -> int | None:
    for part in name.split("_"):
        if part.startswith("kmer"):
            try:
                return int(part[4:])
            except ValueError:
                pass
    return None


def plot_kmer_comparison(histories_dir: Path, figures_dir: Path) -> None:
    """Single figure: mean ± std val_acc and val_loss vs k_mer, one line per model."""
    all_histories = _load_histories(histories_dir)

    kmer_histories = [h for h in all_histories if "_kmer" in h["name"]]
    if not kmer_histories:
        print("[WARN] No k_mer comparison histories found — skipping kmer_comparison figure.")
        return

    data: dict[tuple[str, int], list[tuple[float, float]]] = defaultdict(list)
    for h in kmer_histories:
        model = h["name"].split("_")[0]
        k_mer = _kmer_from_name(h["name"])
        if k_mer is None:
            continue
        last = h["epochs"][-1]
        data[(model, k_mer)].append((last["val_acc"], last["val_loss"]))

    models = sorted({m for m, _ in data})
    k_mers = sorted({k for _, k in data})

    figures_dir.mkdir(parents=True, exist_ok=True)

    with plt.style.context("seaborn-v0_8-whitegrid"):
        fig, (ax_acc, ax_loss) = plt.subplots(1, 2, figsize=(10, 4))

        for model in models:
            color = MODEL_COLORS.get(model, "#333333")
            mean_accs, std_accs, mean_losses, std_losses = [], [], [], []

            for k in k_mers:
                vals = data.get((model, k), [])
                accs = [v[0] for v in vals]
                losses = [v[1] for v in vals]

                mean_accs.append(float(np.mean(accs)) if accs else float("nan"))
                std_accs.append(float(np.std(accs)) if accs else 0.0)
                mean_losses.append(float(np.mean(losses)) if losses else float("nan"))
                std_losses.append(float(np.std(losses)) if losses else 0.0)

                for acc in accs:
                    ax_acc.scatter(k, acc, color=color, alpha=0.3, s=20, zorder=3)
                for loss in losses:
                    ax_loss.scatter(k, loss, color=color, alpha=0.3, s=20, zorder=3)

            ma = np.array(mean_accs)
            sa = np.array(std_accs)
            ml = np.array(mean_losses)
            sl = np.array(std_losses)

            ax_acc.plot(k_mers, ma, color=color, marker="o", label=model)
            ax_acc.fill_between(k_mers, ma - sa, ma + sa, color=color, alpha=0.15)

            ax_loss.plot(k_mers, ml, color=color, marker="o", label=model)
            ax_loss.fill_between(k_mers, ml - sl, ml + sl, color=color, alpha=0.15)

        for ax, ylabel, title in [
            (ax_acc, "Validation Accuracy", "Accuracy vs k-mer size"),
            (ax_loss, "Final Validation Loss", "Validation Loss vs k-mer size"),
        ]:
            ax.set_xlabel("k-mer size")
            ax.set_ylabel(ylabel)
            ax.set_xticks(k_mers)
            ax.set_title(title)
            ax.legend(fontsize=9)

        plt.tight_layout()
        fig.savefig(figures_dir / "kmer_comparison.png", dpi=300, bbox_inches="tight")
        fig.savefig(figures_dir / "kmer_comparison.pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"[INFO] k-mer comparison figure saved to {figures_dir}")


def main() -> None:
    hist_dir = Path("results/histories")
    fig_dir = Path("results/figures")
    plot_training_curves(hist_dir, fig_dir)
    plot_kmer_comparison(hist_dir, fig_dir)


if __name__ == "__main__":
    main()
