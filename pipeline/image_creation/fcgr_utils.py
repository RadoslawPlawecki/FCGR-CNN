from tqdm import tqdm
import numpy as np

from pipeline.config import get_settings
from pipeline.utils import load_fcgr_seq
from pipeline.FCGR import FCGR

settings = get_settings()

def create_fcgr_img(k_mer: int = 6, stop_after: int | None = None) -> None:
    df = load_fcgr_seq(settings.fcgr_file)
    
    ref_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "ref"
    mut_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "mut"
    diff_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "diff"
    
    if not ref_dir.exists() or not mut_dir.exists() or not diff_dir.exists():
        print(f"[INFO] Creating directories for FCGR images at {settings.data_dir / 'fcgr_images'}...")
    
        ref_dir.mkdir(parents=True, exist_ok=True)
        mut_dir.mkdir(parents=True, exist_ok=True)
        diff_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Directories created successfully.")
    else:
        print(f"[INFO] Directories for FCGR images already exist at {settings.data_dir / 'fcgr_images'}.")

    print(f"[INFO] Starting FCGR image creation for {len(df)} rows with k-mer={k_mer}...")
    for i, row in enumerate(tqdm(df.iter_rows(named=True), desc="Creating FCGR images", unit="row")):
        chromosome = row["accession"]
        location = row["loc"]
        label = row["label"]
        ref_img = None
        
        for ref in ["ref_500", "ref_1000", "ref_1500", "ref_2000"]:
            seq_ref = row[ref]
            seq_mut = row[ref.replace("ref", "mut")]
            size = ref.split("_")[1]
            
            fcgr_ref = FCGR(seq_ref, k_mer=k_mer)
            fcgr_mut = FCGR(seq_mut, k_mer=k_mer)
            fcgr_matrix_ref = fcgr_ref.fill_matrix()
            fcgr_matrix_mut = fcgr_mut.fill_matrix()

            ref_img = fcgr_matrix_ref
            np.save(f"{ref_dir}/{chromosome}_{location}_{label}_{size}_ref.npy", ref_img)

            mut_img = fcgr_matrix_mut
            np.save(f"{mut_dir}/{chromosome}_{location}_{label}_{size}_mut.npy", mut_img)

            diff_img = np.abs(ref_img - mut_img)
            np.save(f"{diff_dir}/{chromosome}_{location}_{label}_{size}_diff.npy", diff_img)
            
        if i == stop_after:
            print(f"[INFO] Stopped after processing {stop_after} rows for testing purposes.")
            break
    print(f"[INFO] Finished creating FCGR images for {len(df)} rows.")
    
def show_fcgr_image(k_mer: int, chromosome: str, location: str, label: str, size: str) -> None:
    ref_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "ref"
    mut_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "mut"
    diff_dir = settings.data_dir / "fcgr_images" / f"k_{k_mer}" / "diff"

    ref_img_path = ref_dir / f"{chromosome}_{location}_{label}_{size}_ref.npy"
    mut_img_path = mut_dir / f"{chromosome}_{location}_{label}_{size}_mut.npy"
    diff_img_path = diff_dir / f"{chromosome}_{location}_{label}_{size}_diff.npy"
    
    try:
        ref_img = np.load(ref_img_path)
        mut_img = np.load(mut_img_path)
        diff_img = np.load(diff_img_path)
    except FileNotFoundError as e:
        print(f"File not found: {e.filename}")
        return
    except Exception as e:
        print(f"Error loading images: {e}")
        return
    
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axes[0].imshow(ref_img, cmap='gray')
    axes[0].set_title('Reference Image')
    axes[0].axis('off')
    
    axes[1].imshow(mut_img, cmap='gray')
    axes[1].set_title('Mutated Image')
    axes[1].axis('off')
    
    axes[2].imshow(diff_img, cmap='gray')
    axes[2].set_title('Difference Image')
    axes[2].axis('off')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    create_fcgr_img()
        