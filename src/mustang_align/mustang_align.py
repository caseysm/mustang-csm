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
import shutil
from pathlib import Path
from typing import List, Union
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

class MustangError(Exception):
    """Custom exception for MUSTANG related errors."""
    pass

class MustangAnalyzer:
    """
    Performs MUSTANG analysis on PDB files and generates pairwise alignments.
    """

    def __init__(self, mustang_path: Path = None, alignment_format: str = 'fasta', cpu_percentage: float = 25):
        """
        Initialize the MustangAnalyzer.

        Args:
            mustang_path (Path): Path to the MUSTANG executable. If None, defaults to the 'bin' directory.
            alignment_format (str): Output alignment format (default: 'fasta').
            cpu_percentage (float): Percentage of CPU cores to use (default: 25).
        """
        self.mustang_path = mustang_path or Path(__file__).resolve().parent / "bin" / "mustang-3.2.4"
        self.alignment_format = alignment_format
        self.output_extension = 'afasta' if alignment_format == 'fasta' else alignment_format
        self.cpu_cores = max(1, int(cpu_count() * (cpu_percentage / 100)))

    def run_all_vs_all(self, input_dir: Path, output_dir: Path) -> None:
        """
        Run all-vs-all MUSTANG analysis on PDB files in the input directory.

        Args:
            input_dir (Path): Directory containing input PDB files.
            output_dir (Path): Directory for output alignment files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_output_dir = output_dir / "pdb_outputs"
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        fasta_output_dir = output_dir / "pairwise_fastas"
        fasta_output_dir.mkdir(parents=True, exist_ok=True)

        pdb_files = list(input_dir.glob('*.pdb'))
        if not pdb_files:
            logging.warning(f"No PDB files found in {input_dir}")
            return

        pairs = list(itertools.product(pdb_files, repeat=2))

        with Pool(processes=self.cpu_cores) as pool:
            results = list(tqdm(pool.imap_unordered(self._process_pair_wrapper,
                                [(pdb1, pdb2, output_dir, pdb_output_dir, fasta_output_dir) for pdb1, pdb2 in pairs]),
                                total=len(pairs),
                                desc="Aligning pairs"))

    def run_pdb_vs_pdb(self, pdb1: Path, pdb2: Path, output_dir: Path) -> None:
        """
        Run MUSTANG analysis on two specific PDB files.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            output_dir (Path): Directory for output alignment files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        pdb_output_dir = output_dir / "pdb_outputs"
        pdb_output_dir.mkdir(parents=True, exist_ok=True)
        fasta_output_dir = output_dir / "pairwise_fastas"
        fasta_output_dir.mkdir(parents=True, exist_ok=True)
        self._process_pair(pdb1, pdb2, output_dir, pdb_output_dir, fasta_output_dir)

    def _process_pair_wrapper(self, args):
        """
        Wrapper function for _process_pair to be used with Pool.imap_unordered.
        """
        return self._process_pair(*args)

    def _process_pair(self, pdb1: Path, pdb2: Path, output_dir: Path, pdb_output_dir: Path, fasta_output_dir: Path) -> bool:
        """
        Process a pair of PDB files with MUSTANG.

        Args:
            pdb1 (Path): Path to the first PDB file.
            pdb2 (Path): Path to the second PDB file.
            output_dir (Path): Directory to store temporary output files.
            pdb_output_dir (Path): Directory to store output PDB files.
            fasta_output_dir (Path): Directory to store output FASTA files.

        Returns:
            bool: True if the alignment was successful, False otherwise.
        """
        try:
            return self._run_mustang([pdb1, pdb2], output_dir, pdb_output_dir, fasta_output_dir)
        except Exception as e:
            logging.error(f"Error processing {pdb1} and {pdb2}: {e}")
            return False

    def _run_mustang(self, pdb_files: List[Path], output_dir: Path, pdb_output_dir: Path, fasta_output_dir: Path) -> bool:
        """
        Run MUSTANG on a pair of PDB files.

        Args:
            pdb_files (List[Path]): List of paths to the PDB files to align.
            output_dir (Path): Directory to store temporary output files.
            pdb_output_dir (Path): Directory to store output PDB files.
            fasta_output_dir (Path): Directory to store output FASTA files.

        Returns:
            bool: True if the alignment was successful, False otherwise.

        Raises:
            MustangError: If MUSTANG executable is not found or fails to run.
        """
        output_identifier = output_dir / f"{pdb_files[0].stem}_{pdb_files[1].stem}"
        command = [
            str(self.mustang_path),
            '-i'
        ] + [str(pdb) for pdb in pdb_files] + [
            '-o', str(output_identifier),
            '-F', self.alignment_format
        ]

        try:
            subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            expected_afasta = f"{output_identifier}.{self.output_extension}"
            expected_pdb = f"{output_identifier}.pdb"

            if os.path.exists(expected_afasta) and os.path.exists(expected_pdb):
                # Move and rename the FASTA file
                fasta_dest = fasta_output_dir / f"{pdb_files[0].stem}_{pdb_files[1].stem}.fasta"
                shutil.move(expected_afasta, fasta_dest)

                # Move the PDB file
                pdb_dest = pdb_output_dir / f"{pdb_files[0].stem}_{pdb_files[1].stem}.pdb"
                shutil.move(expected_pdb, pdb_dest)

                return True
            else:
                logging.error(f"MUSTANG did not produce the expected output files for {output_identifier}")
                return False

        except subprocess.CalledProcessError as e:
            logging.error(f"MUSTANG failed for {pdb_files[0]} and {pdb_files[1]}: {e.stderr}")
            return False
        except FileNotFoundError:
            raise MustangError(f"MUSTANG executable not found at {self.mustang_path}")

def mustang_analysis(mode: str, input_path: Union[Path, List[Path]], output_dir: Path, mustang_path: Path = None, alignment_format: str = 'fasta', cpu_percentage: float = 25) -> None:
    """
    Main function to run the MUSTANG analysis for generating alignments.

    Args:
        mode (str): Analysis mode ('all_vs_all' or 'pdb_vs_pdb').
        input_path (Union[Path, List[Path]]): Input directory (for all_vs_all) or list of two PDB files (for pdb_vs_pdb).
        output_dir (Path): Directory for output alignment files.
        mustang_path (Path): Path to the MUSTANG executable (default: None).
        alignment_format (str): Output alignment format (default: 'fasta').
        cpu_percentage (float): Percentage of CPU cores to use (default: 25).
    """
    analyzer = MustangAnalyzer(mustang_path, alignment_format, cpu_percentage)

    if mode == 'all_vs_all':
        if not isinstance(input_path, Path):
            raise ValueError("For all_vs_all mode, input_path should be a directory Path.")
        analyzer.run_all_vs_all(input_path, output_dir)
    elif mode == 'pdb_vs_pdb':
        if not isinstance(input_path, list) or len(input_path) != 2:
            raise ValueError("For pdb_vs_pdb mode, input_path should be a list of two PDB file Paths.")
        analyzer.run_pdb_vs_pdb(input_path[0], input_path[1], output_dir)
    else:
        raise ValueError("Invalid mode. Choose either 'all_vs_all' or 'pdb_vs_pdb'.")

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Perform MUSTANG analysis and generate alignments")
    parser.add_argument("--mode", choices=['all_vs_all', 'pdb_vs_pdb'], required=True, help="Analysis mode")
    parser.add_argument("--output_dir", type=Path, required=True, help="Output directory for alignment results")
    parser.add_argument("--mustang_path", type=Path, help="Path to MUSTANG executable")
    parser.add_argument("--alignment_format", default='fasta', help="Output alignment format (default: fasta)")
    parser.add_argument("--cpu_percentage", type=float, default=25, help="Percentage of CPU cores to use (default: 25)")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input_dir", type=Path, help="Input directory for all_vs_all mode")
    group.add_argument("--input_pdb1", type=Path, help="First input PDB file for pdb_vs_pdb mode")
    group.add_argument("--input_pdb2", type=Path, help="Second input PDB file for pdb_vs_pdb mode")

    args = parser.parse_args()

    if args.mode == 'all_vs_all' and not args.input_dir:
        parser.error("all_vs_all mode requires --input_dir")
    elif args.mode == 'pdb_vs_pdb' and (not args.input_pdb1 or not args.input_pdb2):
        parser.error("pdb_vs_pdb mode requires both --input_pdb1 and --input_pdb2")

    input_path = args.input_dir if args.mode == 'all_vs_all' else [args.input_pdb1, args.input_pdb2]

    mustang_analysis(args.mode, input_path, args.output_dir, args.mustang_path, args.alignment_format, args.cpu_percentage)

if __name__ == "__main__":
    main()
