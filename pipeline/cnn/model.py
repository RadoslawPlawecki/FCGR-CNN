import math

import torch
import torch.nn as nn


class BaseCNN(nn.Module):
    def __init__(self):
        super(BaseCNN, self).__init__()

    def forward(self, x):
        raise NotImplementedError("Subclasses should implement this method.")


class _HaarDWT(nn.Module):
    """Haar low-pass (approximation) coefficients along the H dimension, stride 2.

    Implements left multiplication by the Haar DWT matrix R:
      a[i] = (x[2i] + x[2i+1]) / sqrt(2)
    Output H is floor(H_in / 2).
    """

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        n = x.size(2) // 2
        return (x[:, :, :2 * n:2, :] + x[:, :, 1:2 * n + 1:2, :]) / math.sqrt(2)


class _RandomProjection(nn.Module):
    """Fixed (non-learnable) random projection along the W dimension.

    Implements right multiplication by matrix C^RP of shape (W_in, W_out).
    Entries ~ N(0, 1/W_in) following the Johnson-Lindenstrauss convention.
    """

    def __init__(self, in_dim: int, out_dim: int) -> None:
        super().__init__()
        weight = torch.randn(in_dim, out_dim) / math.sqrt(in_dim)
        self.register_buffer("weight", weight)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, C, H, W) → (B, C, H, out_dim)
        return x @ self.weight  # type: ignore[operator]


class _TwoDSquaredRPDWT(nn.Module):
    """(2D)^2 RP DWT downsampling layer.

    Implements Y = R^{DWT}_{k×d} · X_{d×n} · C^{RP}_{n×h} from:
      Alzu'bi et al., CMC 2022, https://www.techscience.com/cmc/v70n3/44931

    - R^DWT: Haar low-pass transform in H direction  → H // 2
    - C^RP:  Random projection in W direction         → W // 2
    Both spatial dimensions are halved.
    """

    def __init__(self, in_w: int) -> None:
        super().__init__()
        self.dwt = _HaarDWT()
        self.rp = _RandomProjection(in_w, in_w // 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.dwt(x)  # (B, C, H//2, W)
        x = self.rp(x)   # (B, C, H//2, W//2)
        return x


class TwoDSquaredRPDWTCNN(BaseCNN):
    """(2D)^2 RP DWT-CNN from Alzu'bi et al. (CMC 2022, v70n3/44931), Table 2.

    Layers:
      1  Conv(10,  5×5, pad=2) + ReLU + MaxPool(2×2)
      2  Conv(15,  5×5, pad=1) + ReLU + MaxPool(2×2)
      3  Conv(20,  5×5)        + ReLU + (2D)^2 RP DWT
      4  Flatten → Linear(num_classes)

    Args:
        k_mer:       FCGR k-mer size; input image is (2^k × 2^k).
        num_classes: Number of output classes.
        in_channels: Input channels (1 for diff images).
    """

    def __init__(self, k_mer: int = 6, num_classes: int = 2, in_channels: int = 1) -> None:
        super().__init__()

        img = 2 ** k_mer
        h = img + 2 * 2 - 5 + 1   # Conv1 pad=2  (same: h unchanged)
        h = h // 2                  # MaxPool
        h = h + 2 * 1 - 5 + 1      # Conv2 pad=1  (h - 2)
        h = h // 2                  # MaxPool
        h_conv3 = h - 5 + 1         # Conv3 valid  (h - 4)
        h_out = h_conv3 // 2        # DWT: halves H
        w_out = h_conv3 // 2        # RP:  halves W (feature maps are square)
        fc_in = 20 * h_out * w_out

        self.conv_layers = nn.Sequential(
            nn.Conv2d(in_channels, 10, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),
            nn.Conv2d(10, 15, kernel_size=5, padding=1),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),
            nn.Conv2d(15, 20, kernel_size=5),
            nn.ReLU(),
        )
        self.downsample = _TwoDSquaredRPDWT(in_w=h_conv3)
        self.classifier = nn.Linear(fc_in, num_classes)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.conv_layers(x)
        x = self.downsample(x)
        x = torch.flatten(x, 1)
        return self.classifier(x)


class StandardCNN(BaseCNN):
    """Plain CNN with BatchNorm and Dropout.

    Architecture:
      3× (Conv → BN → ReLU → MaxPool)
      Flatten → FC(256) → Dropout → FC(num_classes)

    Args:
        k_mer:       FCGR k-mer size; input is (2^k × 2^k).
        num_classes: Output classes.
        in_channels: Input channels (1 for diff images).
        dropout:     Dropout probability before final FC.
    """

    def __init__(self, k_mer: int = 6, num_classes: int = 2, in_channels: int = 1, dropout: float = 0.5) -> None:
        super().__init__()
        img = 2 ** k_mer  # 64 for k=6
        self.features = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(), nn.MaxPool2d(2),  # 32x32
            nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(), nn.MaxPool2d(2),           # 16x16
            nn.Conv2d(64, 128, 3, padding=1), nn.BatchNorm2d(128), nn.ReLU(), nn.MaxPool2d(2),         # 8x8
        )
        spatial = img // 8
        self.classifier = nn.Sequential(
            nn.Flatten(),
            nn.Linear(128 * spatial * spatial, 256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_classes),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.classifier(self.features(x))


class _ResBlock(nn.Module):
    def __init__(self, channels: int) -> None:
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv2d(channels, channels, 3, padding=1, bias=False),
            nn.BatchNorm2d(channels),
            nn.ReLU(),
            nn.Conv2d(channels, channels, 3, padding=1, bias=False),
            nn.BatchNorm2d(channels),
        )
        self.relu = nn.ReLU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.relu(self.block(x) + x)


class ResNetCNN(BaseCNN):
    """Small ResNet-style CNN for 64×64 FCGR diff images.

    Stem: Conv(32) → BN → ReLU → MaxPool
    Stage 1: ResBlock(32) × 2 → MaxPool
    Stage 2: ResBlock(64) × 2 → MaxPool
    Stage 3: ResBlock(128) × 2 → AdaptiveAvgPool(4×4)
    Head: Flatten → Dropout → Linear(num_classes)
    """

    def __init__(self, k_mer: int = 6, num_classes: int = 2, in_channels: int = 1, dropout: float = 0.4) -> None:  # noqa: ARG002
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1, bias=False),
            nn.BatchNorm2d(32), nn.ReLU(),
            nn.MaxPool2d(2),  # 32×32
        )
        self.stage1 = nn.Sequential(_ResBlock(32), _ResBlock(32), nn.MaxPool2d(2))  # 16×16
        self.proj1 = nn.Sequential(nn.Conv2d(32, 64, 1, bias=False), nn.BatchNorm2d(64))
        self.stage2 = nn.Sequential(_ResBlock(64), _ResBlock(64), nn.MaxPool2d(2))  # 8×8
        self.proj2 = nn.Sequential(nn.Conv2d(64, 128, 1, bias=False), nn.BatchNorm2d(128))
        self.stage3 = nn.Sequential(_ResBlock(128), _ResBlock(128), nn.AdaptiveAvgPool2d(4))  # 4×4
        self.head = nn.Sequential(nn.Flatten(), nn.Dropout(dropout), nn.Linear(128 * 4 * 4, num_classes))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.stem(x)
        x = self.stage1(x)
        x = self.proj1(x)
        x = self.stage2(x)
        x = self.proj2(x)
        x = self.stage3(x)
        return self.head(x)
