"""
MUSTANG.py

This module provides the MustangAnalyzer class for performing MUSTANG analysis on PDB files
and generating pairwise alignments.

The MustangAnalyzer uses the MUSTANG tool to align protein structures and creates pairwise alignment files.
"""

import os
import subprocess
import itertools
import logging
from pathlib import Path
from typing import List
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import argparse

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

class MustangError(Exception):
    """Custom exception for MUSTANG related errors."""
    pass

class MustangAnalyzer:
    """
    Performs MUSTANG analysis on PDB files and generates pairwise alignments.
    """

    def __init__(self, input_dir: Path, output_dir: Path, mustang_path: Path = None, alignment_format: str = 'fasta', cpu_percentage: float = 25):
        """
        Initialize the MustangAnalyzer.

        Args:
            input_dir (Path): Directory containing input PDB files.
            output_dir (Path): Directory for output alignment files.
            mustang_path (Path): Path to the MUSTANG executable. If None, defaults to the 'bin' directory.
            alignment_format (str): Output alignment format (default: 'fasta').
            cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        """
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.mustang_path = mustang_path or Path(__file__).resolve().parent / "bin" / "mustang-3.2.4"
        self.alignment_format = alignment_format
        self.output_extension = 'afasta' if alignment_format == 'fasta' else alignment_format
        self.cpu_cores = max(1, int(cpu_count() * (cpu_percentage / 100)))
        self.expected_fastas = 0
        self.actual_fastas = 0

    def run_analysis(self, pairwise_dir: Path = None) -> None:
        """
        Run the MUSTANG analysis pipeline to generate pairwise alignments.

        Args:
            pairwise_dir (Path): Directory to store pairwise alignment files (default: None).
        """
        pairwise_dir = pairwise_dir or self.output_dir / "pairwise_fastas"
        pairwise_dir.mkdir(parents=True, exist_ok=True)

        pdb_files = list(self.input_dir.glob('*.pdb'))
        if not pdb_files:
            logging.warning(f"No PDB files found in {self.input_dir}")
            return

        self.expected_fastas = len(pdb_files) ** 2
        pairs = list(itertools.product(pdb_files, repeat=2))

        with Pool(processes=self.cpu_cores) as pool:
            results = list(tqdm(pool.imap_unordered(self._process_pair_wrapper,
                                [(pdb1, pdb2, pairwise_dir) for pdb1, pdb2 in pairs]),
                                total=len(pairs),
                                desc="Aligning pairs"))

        self.actual_fastas = sum(results)
        if self.expected_fastas != self.actual_fastas:
            logging.warning(f"Mismatch in alignment file count. Expected: {self.expected_fastas}, Actual: {self.actual_fastas}")

    def _process_pair_wrapper(self, args):
        """
        Wrapper function for _process_pair to be used with Pool.imap_unordered.
        """
        return self._process_pair(*args)

    def _process_pair(self, pdb1: Path, pdb2: Path, pairwise_dir: Path) -> bool:
        """
        Process a pair of PDB files with MUSTANG.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            pairwise_dir (Path): Directory to store pairwise alignment files.

        Returns:
            bool: True if the alignment was successful, False otherwise.
        """
        try:
            return self._run_mustang([pdb1, pdb2], pairwise_dir)
        except Exception as e:
            logging.error(f"Error processing {pdb1} and {pdb2}: {e}")
            return False

    def _run_mustang(self, pdb_files: List[Path], pairwise_dir: Path) -> bool:
        """
        Run MUSTANG on a pair of PDB files.

        Args:
            pdb_files (List[Path]): List of paths to the PDB files to align.
            pairwise_dir (Path): Directory to store pairwise alignment files.

        Returns:
            bool: True if the alignment was successful, False otherwise.

        Raises:
            MustangError: If MUSTANG executable is not found or fails to run.
        """
        output_identifier = pairwise_dir / f"{pdb_files[0].stem}_{pdb_files[1].stem}"
        command = [
            str(self.mustang_path),
            '-i'
        ] + [str(pdb) for pdb in pdb_files] + [
            '-o', str(output_identifier),
            '-F', self.alignment_format
        ]

        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            expected_fasta = f"{output_identifier}.{self.output_extension}"
            expected_pdb = f"{output_identifier}.pdb"

            if os.path.exists(expected_fasta):
                if self.output_extension != self.alignment_format:
                    new_name = f"{output_identifier}.{self.alignment_format}"
                    os.rename(expected_fasta, new_name)

                if os.path.exists(expected_pdb):
                    os.remove(expected_pdb)

                return True
            else:
                logging.error(f"MUSTANG did not produce the expected alignment file: {expected_fasta}")
                return False

        except subprocess.CalledProcessError as e:
            logging.error(f"MUSTANG failed for {pdb_files[0]} and {pdb_files[1]}: {e.stderr}")
            return False
        except FileNotFoundError:
            raise MustangError(f"MUSTANG executable not found at {self.mustang_path}")

def mustang_analysis(input_dir: Path, output_dir: Path, mustang_path: Path = None, alignment_format: str = 'fasta', cpu_percentage: float = 25, pairwise_dir: Path = None) -> None:
    """
    Main function to run the MUSTANG analysis for generating pairwise alignments.

    Args:
        input_dir (Path): Directory containing input PDB files.
        output_dir (Path): Directory for output alignment files.
        mustang_path (Path): Path to the MUSTANG executable (default: None).
        alignment_format (str): Output alignment format (default: 'fasta').
        cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        pairwise_dir (Path): Directory to store pairwise alignment files (default: None).
    """
    analyzer = MustangAnalyzer(input_dir, output_dir, mustang_path, alignment_format, cpu_percentage)
    analyzer.run_analysis(pairwise_dir)

def main():
    parser = argparse.ArgumentParser(description="Perform MUSTANG analysis and generate pairwise alignments")
    parser.add_argument("input_dir", type=Path, help="Directory containing input PDB files")
    parser.add_argument("output_dir", type=Path, help="Output directory for alignment results")
    parser.add_argument("--mustang_path", type=Path, help="Path to MUSTANG executable")
    parser.add_argument("--alignment_format", default='fasta', help="Output alignment format (default: fasta)")
    parser.add_argument("--cpu_percentage", type=float, default=25, help="Percentage of CPU cores to use (default: 25)")
    parser.add_argument("--pairwise_dir", type=Path, help="Directory for pairwise alignments (if different from output_dir)")
    args = parser.parse_args()

    mustang_analysis(args.input_dir, args.output_dir, args.mustang_path, args.alignment_format, args.cpu_percentage, args.pairwise_dir)

if __name__ == "__main__":
    main()



