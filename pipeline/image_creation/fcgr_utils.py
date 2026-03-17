from tqdm import tqdm
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from pipeline.config import get_settings
from pipeline.utils import load_fcgr_seq
from pipeline.FCGR import FCGR

settings = get_settings()

def _process_row(row: dict, k_mer: int, ref_dir: Path, mut_dir: Path, diff_dir: Path, save_original: bool) -> None:
    chromosome = row["accession"]
    location = row["loc"]
    label = row["label"]

    for ref in ["ref_500", "ref_1000", "ref_1500", "ref_2000"]:
        seq_ref = row[ref]
        seq_mut = row[ref.replace("ref", "mut")]
        size = ref.split("_")[1]

        ref_img = FCGR(seq_ref, k_mer=k_mer).fill_matrix()
        mut_img = FCGR(seq_mut, k_mer=k_mer).fill_matrix()
        diff_img = np.abs(ref_img - mut_img)

        if save_original:
            np.save(f"{ref_dir}/{chromosome}_{location}_{label}_{size}_ref.npy", ref_img)
            np.save(f"{mut_dir}/{chromosome}_{location}_{label}_{size}_mut.npy", mut_img)

        np.save(f"{diff_dir}/{chromosome}_{location}_{label}_{size}_diff.npy", diff_img)


def create_fcgr_img(k_mer: int = 6, stop_after: int | None = None, save_original: bool = False) -> None:
    df = load_fcgr_seq(settings.fcgr_file)
    if stop_after is not None:
        df = df.head(stop_after)

    ref_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "ref"
    mut_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "mut"
    diff_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "diff"

    ref_dir.mkdir(parents=True, exist_ok=True)
    mut_dir.mkdir(parents=True, exist_ok=True)
    diff_dir.mkdir(parents=True, exist_ok=True)

    worker = partial(_process_row, k_mer=k_mer, ref_dir=ref_dir, mut_dir=mut_dir,
                     diff_dir=diff_dir, save_original=save_original)

    print(f"[INFO] Starting FCGR image creation for {len(df)} rows with k-mer={k_mer}, workers={settings.workers}...")
    with ProcessPoolExecutor(max_workers=settings.workers) as executor:
        list(tqdm(executor.map(worker, df.iter_rows(named=True)), total=len(df), desc="Creating FCGR images", unit="row"))
        
    print(f"[INFO] Finished creating FCGR images for {len(df)} rows.")
    
def show_fcgr_image(k_mer: int, chromosome: str, location: str, label: str, size: str) -> None:
    ref_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "ref"
    mut_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "mut"
    diff_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "diff"

    ref_img_path = ref_dir / f"{chromosome}_{location}_{label}_{size}_ref.npy"
    mut_img_path = mut_dir / f"{chromosome}_{location}_{label}_{size}_mut.npy"
    diff_img_path = diff_dir / f"{chromosome}_{location}_{label}_{size}_diff.npy"

    images_to_show: list[tuple[str, np.ndarray]] = []

    for title, img_path in [
        ("Reference Image", ref_img_path),
        ("Mutated Image", mut_img_path),
        ("Difference Image", diff_img_path),
    ]:
        if not img_path.exists():
            print(f"File not found: {img_path}")
            continue

        try:
            images_to_show.append((title, np.load(img_path)))
        except Exception as e:
            print(f"Error loading {img_path.name}: {e}")

    if not images_to_show:
        print("No images to display.")
        return
    
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, len(images_to_show), figsize=(5 * len(images_to_show), 5))

    if len(images_to_show) == 1:
        axes = [axes]

    for ax, (title, image) in zip(axes, images_to_show):
        ax.imshow(image, cmap='gray')
        ax.set_title(title)
        ax.axis('off')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    create_fcgr_img()
        