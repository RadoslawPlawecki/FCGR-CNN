"""Single entry point: generate missing FCGR data, run all experiments, plot results."""

from pipeline.config import get_settings
from pipeline.image_creation.fcgr_utils import create_fcgr_img
from pipeline.cnn.experiments import main as run_experiments
from pipeline.cnn.plot import main as run_plots


def main() -> None:
    settings = get_settings()

    for k_mer in [6, 7, 8]:
        if not (settings.data_dir / "fcgr_images" / str(k_mer)).exists():
            print(f"[INFO] Generating FCGR images for k_mer={k_mer}...")
            create_fcgr_img(k_mer=k_mer)

    run_experiments()
    run_plots()


if __name__ == "__main__":
    main()
