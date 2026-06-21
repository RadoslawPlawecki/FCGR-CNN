from pydantic_settings import BaseSettings
from functools import cache
from pathlib import Path
import os

BASE_DIR = Path(__file__).parent.parent

class Settings(BaseSettings):
    data_dir: Path = BASE_DIR / "data"
    fcgr_file: Path = data_dir / "fcgr.csv"
    workers: int = max(1, (os.cpu_count() or 1) - 1)


@cache
def get_settings() -> Settings:
    return Settings()