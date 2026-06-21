import json
import pytest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

from pipeline.cnn.plot import plot_training_curves


def _write_history(path: Path, name: str, n_epochs: int = 5) -> None:
    data = {
        "name": name,
        "epochs": [
            {
                "epoch": i + 1,
                "train_loss": 1.0 - i * 0.1,
                "val_loss": 1.1 - i * 0.09,
                "val_acc": 0.5 + i * 0.05,
            }
            for i in range(n_epochs)
        ],
    }
    path.write_text(json.dumps(data))


def test_plot_training_curves_creates_png_and_pdf(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    for seq_size in ["500", "1000"]:
        _write_history(hist_dir / f"ResNetCNN_seq{seq_size}_undersample.json", f"ResNetCNN_seq{seq_size}_undersample")

    plot_training_curves(hist_dir, fig_dir)

    assert (fig_dir / "curves_ResNetCNN.png").exists()
    assert (fig_dir / "curves_ResNetCNN.pdf").exists()
    assert (fig_dir / "curves_ResNetCNN.png").stat().st_size > 0


def test_plot_training_curves_excludes_kmer_files(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    # Only a kmer file — should produce no figure
    _write_history(hist_dir / "ResNetCNN_kmer6_seq500_undersample.json", "ResNetCNN_kmer6_seq500_undersample")

    plot_training_curves(hist_dir, fig_dir)

    assert not (fig_dir / "curves_ResNetCNN.png").exists()


def test_plot_training_curves_one_figure_per_model(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    for model in ["ResNetCNN", "StandardCNN"]:
        _write_history(hist_dir / f"{model}_seq500_undersample.json", f"{model}_seq500_undersample")

    plot_training_curves(hist_dir, fig_dir)

    assert (fig_dir / "curves_ResNetCNN.png").exists()
    assert (fig_dir / "curves_StandardCNN.png").exists()


def test_plot_training_curves_missing_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        plot_training_curves(tmp_path / "nonexistent", tmp_path / "figures")


def test_plot_training_curves_skips_malformed_json(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    (hist_dir / "bad.json").write_text("not json{{{")
    _write_history(hist_dir / "ResNetCNN_seq500_undersample.json", "ResNetCNN_seq500_undersample")

    # Should not raise — bad file is skipped, good file produces figure
    plot_training_curves(hist_dir, fig_dir)
    assert (fig_dir / "curves_ResNetCNN.png").exists()


from pipeline.cnn.plot import plot_kmer_comparison


def test_plot_kmer_comparison_creates_png_and_pdf(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    for k_mer in [6, 7, 8]:
        for seq_size in ["500", "1000", "1500", "2000"]:
            name = f"ResNetCNN_kmer{k_mer}_seq{seq_size}_undersample"
            _write_history(hist_dir / f"{name}.json", name)

    plot_kmer_comparison(hist_dir, fig_dir)

    assert (fig_dir / "kmer_comparison.png").exists()
    assert (fig_dir / "kmer_comparison.pdf").exists()
    assert (fig_dir / "kmer_comparison.png").stat().st_size > 0


def test_plot_kmer_comparison_excludes_non_kmer_files(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    # Only non-kmer files — should produce no figure
    _write_history(hist_dir / "ResNetCNN_seq500_undersample.json", "ResNetCNN_seq500_undersample")

    plot_kmer_comparison(hist_dir, fig_dir)

    assert not (fig_dir / "kmer_comparison.png").exists()


def test_plot_kmer_comparison_missing_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        plot_kmer_comparison(tmp_path / "nonexistent", tmp_path / "figures")


def test_plot_kmer_comparison_multiple_models(tmp_path):
    hist_dir = tmp_path / "histories"
    hist_dir.mkdir()
    fig_dir = tmp_path / "figures"

    for model in ["ResNetCNN", "StandardCNN"]:
        for k_mer in [6, 7, 8]:
            for seq_size in ["500", "1000"]:
                name = f"{model}_kmer{k_mer}_seq{seq_size}_undersample"
                _write_history(hist_dir / f"{name}.json", name)

    plot_kmer_comparison(hist_dir, fig_dir)

    assert (fig_dir / "kmer_comparison.png").exists()
