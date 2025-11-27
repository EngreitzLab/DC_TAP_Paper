import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def file_exist(fp: str) -> bool:
    """
    Method to check if a file exists
    """
    return (PROJECT_ROOT / fp).exists()


def mkdir(fp: str):
    """
    Method to make directory within python env and scripts.
    """
    subprocess.run(["mkdir", "-p", (PROJECT_ROOT / fp)])
