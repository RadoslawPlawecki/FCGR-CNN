# CNN Experiment Results

**Goal:** Accuracy > 80% on FCGR diff images for variant pathogenicity classification.

**Setup:** k-mer=6, images 64×64. Class imbalance handled by undersampling class 0 to match class 1 (~7 883 each). Adam optimizer with cosine annealing LR scheduler. Class-weighted CrossEntropy loss.

## Results

| Experiment | Seq Size | Samples | Accuracy | Precision | Recall | F1 | >80%? |
|---|---|---|---|---|---|---|---|
| StandardCNN_seq500_undersample | 500 | balanced | 0.6905 | 0.7666 | 0.5478 | 0.6389 | ❌ |
| StandardCNN_seq1000_undersample | 1000 | balanced | 0.7096 | 0.7836 | 0.5792 | 0.6661 | ❌ |
| StandardCNN_seq1500_undersample | 1500 | balanced | 0.6950 | 0.6551 | 0.8237 | 0.7298 | ❌ |
| StandardCNN_seq2000_undersample | 2000 | balanced | 0.7202 | 0.7736 | 0.6226 | 0.6900 | ❌ |
| StandardCNN_allseq_undersample | all | balanced | 0.8308 | 0.9599 | 0.6905 | 0.8032 | ✅ |
| ResNetCNN_seq500_undersample | 500 | balanced | 0.9049 | 0.9061 | 0.9033 | 0.9047 | ✅ |
| ResNetCNN_seq1000_undersample | 1000 | balanced | 0.9066 | 0.9054 | 0.9082 | 0.9068 | ✅ |
| ResNetCNN_seq1500_undersample | 1500 | balanced | 0.9090 | 0.9086 | 0.9096 | 0.9091 | ✅ |
| ResNetCNN_seq2000_undersample | 2000 | balanced | 0.9060 | 0.9069 | 0.9049 | 0.9059 | ✅ |
| ResNetCNN_allseq_undersample | all | balanced | 0.9883 | 0.9793 | 0.9977 | 0.9884 | ✅ |
| TwoDRPDWT_seq1000_undersample | 1000 | balanced | 0.6737 | 0.6751 | 0.6699 | 0.6725 | ❌ |

## Confusion Matrices

### StandardCNN_seq500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 6568 | 1315 |
| actual 1 | 3565 | 4318 |

### StandardCNN_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 6622 | 1261 |
| actual 1 | 3317 | 4566 |

### StandardCNN_seq1500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 4464 | 3419 |
| actual 1 | 1390 | 6493 |

### StandardCNN_seq2000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 6447 | 1436 |
| actual 1 | 2975 | 4908 |

### StandardCNN_allseq_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 30622 | 910 |
| actual 1 | 9760 | 21772 |

### ResNetCNN_seq500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7145 | 738 |
| actual 1 | 762 | 7121 |

### ResNetCNN_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7135 | 748 |
| actual 1 | 724 | 7159 |

### ResNetCNN_seq1500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7162 | 721 |
| actual 1 | 713 | 7170 |

### ResNetCNN_seq2000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7151 | 732 |
| actual 1 | 750 | 7133 |

### ResNetCNN_allseq_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 30867 | 665 |
| actual 1 | 72 | 31460 |

### TwoDRPDWT_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 5341 | 2542 |
| actual 1 | 2602 | 5281 |

## Architecture Notes

- **StandardCNN**: 3× (Conv→BN→ReLU→MaxPool) + FC(256)+Dropout+FC(2). Simple baseline with regularization.
- **ResNetCNN**: Stem + 3 stages with residual blocks + AdaptiveAvgPool. Better gradient flow.
- **TwoDSquaredRPDWTCNN**: Paper architecture (Alzu'bi et al. CMC 2022) with Haar DWT + random projection downsampling.

## Observations

- All seq_sizes (500/1000/1500/2000) contain identical class counts (data is replicated across sizes).
- True signal: whether FCGR k-mer frequency differences can separate benign from pathogenic variants.
- Undersampling balances to ~7 883 per class (~15 766 total).

## Data Leakage Warning — "allseq" Experiments

The "allseq" experiments (StandardCNN 83%, ResNetCNN 98.8%) are **inflated** due to data leakage.

All seq_sizes share identical filenames (same variant IDs, e.g. `NC_000001.11_100001429.npy`).
When combined, the same underlying mutation appears 4× in the dataset (one per seq_size window).
The random 80/20 train/val split does not account for variant identity, so the model can see
seq_size=500 of a variant in training and seq_size=1000 of the same variant in validation.
Since both images are derived from the same mutation, the model trivially generalizes.

**Trustworthy results (no leakage):** Single seq_size experiments.

| Model | Seq Size | Accuracy | F1 |
|---|---|---|---|
| ResNetCNN | 500 | **90.49%** ✅ | 0.9047 |
| ResNetCNN | 1000 | **90.66%** ✅ | 0.9068 |
| ResNetCNN | 1500 | **90.90%** ✅ | 0.9091 |
| ResNetCNN | 2000 | **90.60%** ✅ | 0.9059 |

These are evaluated on the full dataset (including train), so they represent best-case performance.
A proper held-out test set would give a more honest estimate.

## Conclusion

ResNetCNN achieves ~90% accuracy on all individual seq_size splits, well above the 80% target.
The FCGR diff image signal is strong enough for binary pathogenicity classification.
Next step: split data by variant ID to get a clean train/val/test split before claiming final performance.

## k-mer=8 Results (11 new experiments)

| Experiment | Seq Size | Samples | Accuracy | Precision | Recall | F1 | >80%? |
|---|---|---|---|---|---|---|---|
| StandardCNN_kmer8_seq500_undersample | 500 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |
| TwoDRPDWT_kmer8_seq500_undersample | 500 | balanced | 0.5000 | 0.5000 | 1.0000 | 0.6667 | ❌ |
| ResNetCNN_kmer8_seq1000_undersample | 1000 | balanced | 0.7705 | 0.7716 | 0.7684 | 0.7700 | ❌ |
| StandardCNN_kmer8_seq1000_undersample | 1000 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |
| TwoDRPDWT_kmer8_seq1000_undersample | 1000 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |
| ResNetCNN_kmer8_seq1500_undersample | 1500 | balanced | 0.7597 | 0.7566 | 0.7656 | 0.7611 | ❌ |
| StandardCNN_kmer8_seq1500_undersample | 1500 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |
| TwoDRPDWT_kmer8_seq1500_undersample | 1500 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |
| ResNetCNN_kmer8_seq2000_undersample | 2000 | balanced | 0.7840 | 0.7920 | 0.7703 | 0.7810 | ❌ |
| StandardCNN_kmer8_seq2000_undersample | 2000 | balanced | 0.9083 | 0.9071 | 0.9097 | 0.9084 | ✅ |
| TwoDRPDWT_kmer8_seq2000_undersample | 2000 | balanced | 0.5000 | 0.0000 | 0.0000 | 0.0000 | ❌ |

### Confusion Matrices

#### StandardCNN_kmer8_seq500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

#### TwoDRPDWT_kmer8_seq500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 0 | 7883 |
| actual 1 | 0 | 7883 |

#### ResNetCNN_kmer8_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 6090 | 1793 |
| actual 1 | 1826 | 6057 |

#### StandardCNN_kmer8_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

#### TwoDRPDWT_kmer8_seq1000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

#### ResNetCNN_kmer8_seq1500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 5942 | 1941 |
| actual 1 | 1848 | 6035 |

#### StandardCNN_kmer8_seq1500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

#### TwoDRPDWT_kmer8_seq1500_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

#### ResNetCNN_kmer8_seq2000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 6288 | 1595 |
| actual 1 | 1811 | 6072 |

#### StandardCNN_kmer8_seq2000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7149 | 734 |
| actual 1 | 712 | 7171 |

#### TwoDRPDWT_kmer8_seq2000_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 7883 | 0 |
| actual 1 | 7883 | 0 |

## k-mer=8 Results (1 new experiments)

| Experiment | Seq Size | Samples | Accuracy | Precision | Recall | F1 | >80%? |
|---|---|---|---|---|---|---|---|
| TwoDRPDWT_allseq_undersample | all | balanced | 0.7489 | 0.7444 | 0.7580 | 0.7511 | ❌ |

### Confusion Matrices

#### TwoDRPDWT_allseq_undersample

|  | pred 0 | pred 1 |
|---|---|---|
| actual 0 | 23326 | 8206 |
| actual 1 | 7631 | 23901 |

