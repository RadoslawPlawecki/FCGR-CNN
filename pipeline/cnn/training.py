from dataclasses import dataclass, field
from typing import Literal

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, random_split


@dataclass
class TrainingConfig:
    """Training hyperparameters."""
    batch_size: int = 64
    epochs: int = 30
    lr: float = 1e-3
    momentum: float = 0.9
    val_split: float = 0.2
    optimizer: Literal["sgd", "adam"] = "adam"
    weight_decay: float = 1e-4

    def build_optimizer(self, model: nn.Module) -> optim.Optimizer:
        if self.optimizer == "adam":
            return optim.Adam(model.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        return optim.SGD(model.parameters(), lr=self.lr, momentum=self.momentum, weight_decay=self.weight_decay)


@dataclass
class EvalMetrics:
    accuracy: float
    precision: float
    recall: float
    f1: float
    tp: int
    tn: int
    fp: int
    fn: int


def fit(
    model: nn.Module,
    dataset: Dataset,
    config: TrainingConfig | None = None,
    device: torch.device | None = None,
    verbose: bool = True,
) -> tuple[nn.Module, list[dict]]:
    """Train model; returns (trained_model, epoch_history)."""
    if config is None:
        config = TrainingConfig()
    if device is None:
        device = torch.device(
            "cuda" if torch.cuda.is_available()
            else "mps" if torch.backends.mps.is_available()
            else "cpu"
        )

    val_size = int(len(dataset) * config.val_split)  # type: ignore[arg-type]
    train_size = len(dataset) - val_size  # type: ignore[arg-type]
    train_ds, val_ds = random_split(dataset, [train_size, val_size], generator=torch.Generator().manual_seed(42))

    train_loader = DataLoader(train_ds, batch_size=config.batch_size, shuffle=True, num_workers=4, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size, num_workers=4, pin_memory=True)

    train_labels = torch.tensor([train_ds[i][1] for i in range(len(train_ds))])
    num_classes = 2
    class_weights = torch.zeros(num_classes, device=device)
    for c in range(num_classes):
        count = (train_labels == c).sum().item()
        class_weights[c] = len(train_labels) / (num_classes * count) if count > 0 else 1.0

    model = model.to(device)
    optimizer = config.build_optimizer(model)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=config.epochs)

    if verbose:
        print(f"[INFO] Device={device} | train={train_size} val={val_size} | opt={config.optimizer} lr={config.lr}")
        print(f"[INFO] Class weights: 0={class_weights[0]:.3f}  1={class_weights[1]:.3f}")

    history = []
    for epoch in range(1, config.epochs + 1):
        model.train()
        train_loss = 0.0
        for X, y in train_loader:
            X, y = X.to(device), y.to(device)
            optimizer.zero_grad()
            loss = criterion(model(X), y)
            loss.backward()
            optimizer.step()
            train_loss += loss.item() * len(X)
        scheduler.step()

        model.eval()
        val_loss, correct = 0.0, 0
        with torch.no_grad():
            for X, y in val_loader:
                X, y = X.to(device), y.to(device)
                out = model(X)
                val_loss += criterion(out, y).item() * len(X)
                correct += (out.argmax(1) == y).sum().item()

        val_acc = correct / val_size
        entry = {
            "epoch": epoch,
            "train_loss": train_loss / train_size,
            "val_loss": val_loss / val_size,
            "val_acc": val_acc,
        }
        history.append(entry)
        if verbose:
            print(
                f"Epoch {epoch:>3}/{config.epochs} | "
                f"train_loss={entry['train_loss']:.4f} | "
                f"val_loss={entry['val_loss']:.4f} | "
                f"val_acc={val_acc:.4f}"
            )

    return model, history


def evaluate(
    model: nn.Module,
    dataset: Dataset,
    config: TrainingConfig | None = None,
    device: torch.device | None = None,
    verbose: bool = True,
) -> EvalMetrics:
    """Evaluate model on full dataset; returns EvalMetrics."""
    if config is None:
        config = TrainingConfig()
    if device is None:
        device = torch.device(
            "cuda" if torch.cuda.is_available()
            else "mps" if torch.backends.mps.is_available()
            else "cpu"
        )

    loader = DataLoader(dataset, batch_size=config.batch_size, num_workers=4, pin_memory=True)
    model.eval().to(device)

    all_preds: list[int] = []
    all_labels: list[int] = []
    with torch.no_grad():
        for X, y in loader:
            preds = model(X.to(device)).argmax(1).cpu().tolist()
            all_preds.extend(preds)
            all_labels.extend(y.tolist())

    preds_t = torch.tensor(all_preds)
    labels_t = torch.tensor(all_labels)

    tp = int(((preds_t == 1) & (labels_t == 1)).sum())
    tn = int(((preds_t == 0) & (labels_t == 0)).sum())
    fp = int(((preds_t == 1) & (labels_t == 0)).sum())
    fn = int(((preds_t == 0) & (labels_t == 1)).sum())
    total = len(all_labels)

    accuracy  = (tp + tn) / total
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1        = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

    metrics = EvalMetrics(accuracy=accuracy, precision=precision, recall=recall, f1=f1,
                          tp=tp, tn=tn, fp=fp, fn=fn)

    if verbose:
        print("\n=== Evaluation Results ===")
        print(f"  Accuracy:  {accuracy:.4f}")
        print(f"  Precision: {precision:.4f}")
        print(f"  Recall:    {recall:.4f}")
        print(f"  F1:        {f1:.4f}")
        print(f"\nConfusion matrix (rows=actual, cols=predicted):")
        print(f"              pred 0   pred 1")
        print(f"  actual 0    {tn:>6}   {fp:>6}")
        print(f"  actual 1    {fn:>6}   {tp:>6}")

    return metrics
