import random
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import Dataset


class FCGRDiffDataset(Dataset[tuple[torch.Tensor, int]]):
    """Loads diff FCGR images from data/fcgr_images/{k_mer}/{seq_size}/diff/{0,1}/.

    Args:
        k_mer:       FCGR k-mer size (determines subdirectory).
        seq_size:    Sequence window size string, e.g. "500", "1000", "1500", "2000".
                     Pass None to combine all available seq_sizes.
        data_dir:    Root data directory (settings.data_dir).
        undersample: If True, downsample the majority class to match the minority class.
        seed:        Random seed for undersampling reproducibility.
    """

    def __init__(
        self,
        k_mer: int,
        seq_size: str | None,
        data_dir: Path,
        undersample: bool = False,
        seed: int = 42,
    ) -> None:
        by_label: dict[int, list[Path]] = {0: [], 1: []}

        seq_sizes = (
            [seq_size]
            if seq_size is not None
            else [d.name for d in (data_dir / "fcgr_images" / str(k_mer)).iterdir() if d.is_dir()]
        )

        for ss in seq_sizes:
            for label in (0, 1):
                label_dir = data_dir / "fcgr_images" / str(k_mer) / ss / "diff" / str(label)
                if not label_dir.exists():
                    continue
                by_label[label].extend(sorted(label_dir.glob("*.npy")))

        if undersample:
            min_count = min(len(by_label[0]), len(by_label[1]))
            rng = random.Random(seed)
            for label in (0, 1):
                if len(by_label[label]) > min_count:
                    by_label[label] = rng.sample(by_label[label], min_count)

        self.samples: list[tuple[Path, int]] = []
        for label, paths in by_label.items():
            for p in paths:
                self.samples.append((p, label))

        if not self.samples:
            raise FileNotFoundError(
                f"No .npy files found under {data_dir}/fcgr_images/{k_mer}/<seq_size>/diff/{{0,1}}"
            )

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> tuple[torch.Tensor, int]:
        path, label = self.samples[idx]
        img = np.load(path).astype(np.float32)
        return torch.from_numpy(img).unsqueeze(0), label  # (1, H, W), int
