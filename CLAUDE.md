# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
uv sync                          # install dependencies
uv run pytest                    # run all tests
uv run pytest tests/path/test.py # run single test file
python -m pipeline.image_creation.fcgr_utils   # generate FCGR diff images
python -m pipeline.cnn.run       # train and evaluate the CNN
```

## Architecture

This project classifies DNA variants as benign/pathogenic using Frequency Chaos Game Representation (FCGR) images fed into a CNN.

**Data flow:**

1. `data/fcgr.csv` — raw variant data (reference + mutated sequences with ClinVar labels)
2. `pipeline/image_creation/fcgr_utils.py` — converts sequences into FCGR difference images saved as `.npy` files
3. `data/fcgr_images/{k_mer}/{seq_size}/diff/{0,1}/` — label-partitioned dataset (0=benign, 1=pathogenic)
4. `pipeline/cnn/` — loads images, trains, evaluates

**CNN pipeline (`pipeline/cnn/`):**

- `dataset.py` — `FCGRDiffDataset`: loads `.npy` files from the label dirs; path template `data/fcgr_images/{k_mer}/{seq_size}/diff/{label}/`
- `model.py` — `TwoDSquaredRPDWTCNN`: implements Alzu'bi et al. (CMC 2022) — 3 conv layers then a custom (2D)² RP-DWT downsampling layer (Haar DWT on H axis + random projection on W axis) instead of a final pooling layer
- `training.py` — `fit()` uses SGD+momentum with inverse-frequency class weights baked into CrossEntropyLoss; `evaluate()` prints accuracy/precision/recall/F1 + confusion matrix
- `run.py` — entry point wiring dataset → model → fit → evaluate

**Config:** `pipeline/config.py` uses `pydantic-settings`; override `data_dir` via env var.

**Current data:** only k_mer=6 is generated; seq_sizes 500/1000/1500/2000 each have ~17 181 class-0 and ~7 883 class-1 samples (~2.2:1 imbalance).

**Model reference:** Alzu'bi et al., *CMC* 2022, v70n3/44931 — the (2D)² RP-DWT-CNN architecture. Input image size is `2^k_mer × 2^k_mer` (64×64 for k=6).
