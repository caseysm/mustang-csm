from .mustang_align import MustangAnalyzer, mustang_analysis, MustangError, main

__all__ = ['MustangAnalyzer', 'mustang_analysis', 'MustangError', 'main']

# Convenience import for easier access to Path
from pathlib import Path

def get_mustang_path():
    """
    Returns the path to the MUSTANG executable included with this package.

    Returns:
        Path: Path to the MUSTANG executable.
    """
    return Path(__file__).parent / "bin" / "mustang-3.2.4"

