from pydantic_settings import BaseSettings
from functools import cache
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent

class Settings(BaseSettings):
    data_dir: Path = BASE_DIR / "data"
    fcgr_file: Path = data_dir / "fcgr.csv"


@cache
def get_settings() -> Settings:
    return Settings()