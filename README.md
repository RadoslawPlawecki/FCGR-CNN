# FCGR-CNN

Frequency Chaos Game Representation (FCGR) image pipeline for CNN-based DNA variant classification.

## Setup

Requires Python 3.13+.

### uv

```bash
uv sync
```

### pip

```bash
pip install -r requirements.txt
```

### conda

```bash
conda create -n fcgr-cnn python=3.13
conda activate fcgr-cnn
pip install -r requirements.txt
```

## Data

Before running the image creation pipeline, download the FCGR sequence CSV file:

**[Download fcgr.csv](https://drive.google.com/file/d/1Ngn8EO_pBfItN6VPO0X0wpUaWY_HuWmR/view?usp=sharing)**

Place the file at:

```
data/fcgr.csv
```

## Usage

### Generate FCGR images

```bash
python -m pipeline.image_creation.fcgr_utils
```

This will create FCGR difference images under `data/fcgr_images/k_<k>/diff/`, which are used as CNN input.

Optional parameters in `create_fcgr_img()`:

| Parameter | Default | Description |
|---|---|---|
| `k_mer` | `6` | k-mer length (2–6) |
| `stop_after` | `None` | Process only first N rows (for testing) |
| `save_original` | `False` | Also save ref and mut images |
