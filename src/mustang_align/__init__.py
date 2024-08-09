"""
MUSTANG Align

This package provides tools for performing MUSTANG analysis on PDB files
and generating pairwise alignments.

Classes:
    MustangAnalyzer: Main class for performing MUSTANG analysis.

Functions:
    mustang_analysis: Convenience function for running MUSTANG analysis.
    main: Entry point for command-line interface.

For command-line usage, use the 'mustang_align' command after installation.
"""

from .mustang_align import MustangAnalyzer, mustang_analysis, MustangError, main

__version__ = "0.1.2"
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


# Example usage
__usage__ = """
# Example usage of mustang_align

from mustang_align import mustang_analysis, Path

input_dir = Path("path/to/pdb/files")
output_dir = Path("path/to/output")

mustang_analysis(input_dir, output_dir)

# For command-line usage:
# mustang_align input_dir output_dir [options]
"""