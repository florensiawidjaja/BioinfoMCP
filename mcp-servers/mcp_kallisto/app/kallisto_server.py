from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Union
import subprocess

mcp = FastMCP()


@mcp.tool()
def index(
    fasta_files: List[Path],
    index: Path,
    kmer_size: int = 31,
    d_list: Optional[Path] = None,
    make_unique: bool = False,
    aa: bool = False,
    distinguish: bool = False,
    threads: int = 1,
    min_size: Optional[int] = None,
    ec_max_size: Optional[int] = None,
):
    """
    Builds a kallisto index from a FASTA formatted file of target sequences.

    Parameters:
    - fasta_files: List of FASTA files (plaintext or gzipped) containing transcriptome sequences.
    - index: Filename for the kallisto index to be constructed.
    - kmer_size: k-mer (odd) length (default: 31, max: 31).
    - d_list: Path to a FASTA file containing sequences to mask from quantification.
    - make_unique: Replace repeated target names with unique names.
    - aa: Generate index from a FASTA file containing amino acid sequences.
    - distinguish: Generate index where sequences are distinguished by the sequence name.
    - threads: Number of threads to use (default: 1).
    - min_size: Length of minimizers (default: automatically chosen).
    - ec_max_size: Maximum number of targets in an equivalence class (default: no maximum).
    """
    # Validate fasta_files
    if not fasta_files or len(fasta_files) == 0:
        raise ValueError("At least one FASTA file must be provided in fasta_files.")
    for f in fasta_files:
        if not f.exists():
            raise FileNotFoundError(f"FASTA file not found: {f}")

    # Validate index path parent directory exists
    if not index.parent.exists():
        raise FileNotFoundError(f"Index output directory does not exist: {index.parent}")

    # Validate kmer_size
    if kmer_size < 1 or kmer_size > 31 or kmer_size % 2 == 0:
        raise ValueError("kmer_size must be an odd integer between 1 and 31 (inclusive).")

    # Validate threads
    if threads < 1:
        raise ValueError("threads must be >= 1.")

    # Validate min_size if given
    if min_size is not None and min_size < 1:
        raise ValueError("min_size must be >= 1 if specified.")

    # Validate ec_max_size if given
    if ec_max_size is not None and ec_max_size < 1:
        raise ValueError("ec_max_size must be >= 1 if specified.")

    cmd = ["kallisto", "index", "-i", str(index), "-k", str(kmer_size)]
    if d_list:
        if not d_list.exists():
            raise FileNotFoundError(f"d_list FASTA file not found: {d_list}")
        cmd += ["-d", str(d_list)]
    if make_unique:
        cmd.append("--make-unique")
    if aa:
        cmd.append("--aa")
    if distinguish:
        cmd.append("--distinguish")
    if threads != 1:
        cmd += ["-t", str(threads)]
    if min_size is not None:
        cmd += ["-m", str(min_size)]
    if ec_max_size is not None:
        cmd += ["-e", str(ec_max_size)]

    # Add fasta files at the end
    cmd += [str(f) for f in fasta_files]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(index)],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto index failed with exit code {e.returncode}",
        }


@mcp.tool()
def quant(
    fastq_files: List[Path],
    index: Path,
    output_dir: Path,
    bootstrap_samples: int = 0,
    seed: int = 42,
    plaintext: bool = False,
    single: bool = False,
    single_overhang: bool = False,
    fr_stranded: bool = False,
    rf_stranded: bool = False,
    fragment_length: Optional[float] = None,
    sd: Optional[float] = None,
    threads: int = 1,
):
    """
    Runs the quantification algorithm on FASTQ files using a kallisto index.

    Parameters:
    - fastq_files: List of FASTQ files (plaintext or gzipped). For paired-end, provide pairs in order.
    - index: Filename for the kallisto index to be used for quantification.
    - output_dir: Directory to write output to.
    - bootstrap_samples: Number of bootstrap samples (default: 0).
    - seed: Seed for bootstrap sampling (default: 42).
    - plaintext: Output plaintext instead of HDF5.
    - single: Quantify single-end reads.
    - single_overhang: Include reads where unobserved rest of fragment is predicted outside transcript.
    - fr_stranded: Strand specific reads, first read forward.
    - rf_stranded: Strand specific reads, first read reverse.
    - fragment_length: Estimated average fragment length (required if single).
    - sd: Estimated standard deviation of fragment length (required if single).
    - threads: Number of threads to use (default: 1).
    """
    # Validate fastq_files
    if not fastq_files or len(fastq_files) == 0:
        raise ValueError("At least one FASTQ file must be provided in fastq_files.")
    for f in fastq_files:
        if not f.exists():
            raise FileNotFoundError(f"FASTQ file not found: {f}")

    # Validate index file
    if not index.exists():
        raise FileNotFoundError(f"Index file not found: {index}")

    # Validate output_dir exists or create it
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    # Validate bootstrap_samples
    if bootstrap_samples < 0:
        raise ValueError("bootstrap_samples must be >= 0.")

    # Validate seed
    if seed < 0:
        raise ValueError("seed must be >= 0.")

    # Validate threads
    if threads < 1:
        raise ValueError("threads must be >= 1.")

    # Validate single-end parameters
    if single:
        if fragment_length is None or fragment_length <= 0:
            raise ValueError("fragment_length must be > 0 when using single-end mode.")
        if sd is None or sd <= 0:
            raise ValueError("sd must be > 0 when using single-end mode.")
    else:
        # For paired-end, number of fastq files must be even
        if len(fastq_files) % 2 != 0:
            raise ValueError("For paired-end mode, an even number of FASTQ files must be provided.")

    cmd = [
        "kallisto",
        "quant",
        "-i",
        str(index),
        "-o",
        str(output_dir),
        "-t",
        str(threads),
    ]

    if bootstrap_samples != 0:
        cmd += ["-b", str(bootstrap_samples)]
    if seed != 42:
        cmd += ["--seed", str(seed)]
    if plaintext:
        cmd.append("--plaintext")
    if single:
        cmd.append("--single")
    if single_overhang:
        cmd.append("--single-overhang")
    if fr_stranded:
        cmd.append("--fr-stranded")
    if rf_stranded:
        cmd.append("--rf-stranded")
    if single:
        cmd += ["-l", str(fragment_length), "-s", str(sd)]

    # Add fastq files at the end
    cmd += [str(f) for f in fastq_files]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Output files expected:
        # abundance.h5 (unless plaintext), abundance.tsv, run_info.json
        output_files = [
            str(output_dir / "abundance.tsv"),
            str(output_dir / "run_info.json"),
        ]
        if not plaintext:
            output_files.append(str(output_dir / "abundance.h5"))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto quant failed with exit code {e.returncode}",
        }


@mcp.tool()
def quant_tcc(
    tcc_matrix: Path,
    output_dir: Path,
    bootstrap_samples: int = 0,
    seed: int = 42,
    plaintext: bool = False,
    threads: int = 1,
):
    """
    Runs quantification on transcript-compatibility counts (TCC) matrix file.

    Parameters:
    - tcc_matrix: Path to the transcript-compatibility-counts matrix file (MatrixMarket format).
    - output_dir: Directory to write output to.
    - bootstrap_samples: Number of bootstrap samples (default: 0).
    - seed: Seed for bootstrap sampling (default: 42).
    - plaintext: Output plaintext instead of HDF5.
    - threads: Number of threads to use (default: 1).
    """
    if not tcc_matrix.exists():
        raise FileNotFoundError(f"TCC matrix file not found: {tcc_matrix}")

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    if bootstrap_samples < 0:
        raise ValueError("bootstrap_samples must be >= 0.")

    if seed < 0:
        raise ValueError("seed must be >= 0.")

    if threads < 1:
        raise ValueError("threads must be >= 1.")

    cmd = [
        "kallisto",
        "quant-tcc",
        "-t",
        str(threads),
        "-b",
        str(bootstrap_samples),
        "--seed",
        str(seed),
    ]

    if plaintext:
        cmd.append("--plaintext")

    cmd += [str(tcc_matrix)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # quant-tcc output files are not explicitly documented, assume output_dir contains results
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output_dir)],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto quant-tcc failed with exit code {e.returncode}",
        }


@mcp.tool()
def bus(
    fastq_files: List[Path],
    output_dir: Path,
    index: Optional[Path] = None,
    txnames: Optional[Path] = None,
    ec_file: Optional[Path] = None,
    fragment_file: Optional[Path] = None,
    long: bool = False,
    platform: Optional[str] = None,
    fragment_length: Optional[float] = None,
    sd: Optional[float] = None,
    threads: int = 1,
    genemap: Optional[Path] = None,
    gtf: Optional[Path] = None,
    bootstrap_samples: int = 0,
    matrix_to_files: bool = False,
    matrix_to_directories: bool = False,
    seed: int = 42,
    plaintext: bool = False,
):
    """
    Generates BUS files for single-cell sequencing from FASTQ files.

    Parameters:
    - fastq_files: List of FASTQ files (plaintext or gzipped).
    - output_dir: Directory to write output to.
    - index: Filename for the kallisto index to be used.
    - txnames: File with names of transcripts (required if index not supplied).
    - ec_file: File containing equivalence classes (default: from index).
    - fragment_file: File containing fragment length distribution.
    - long: Use version of EM for long reads.
    - platform: Sequencing platform (e.g., PacBio or ONT).
    - fragment_length: Estimated average fragment length.
    - sd: Estimated standard deviation of fragment length.
    - threads: Number of threads to use (default: 1).
    - genemap: File for mapping transcripts to genes.
    - gtf: GTF file for transcriptome information.
    - bootstrap_samples: Number of bootstrap samples (default: 0).
    - matrix_to_files: Reorganize matrix output into abundance tsv files.
    - matrix_to_directories: Reorganize matrix output into abundance tsv files across multiple directories.
    - seed: Seed for bootstrap sampling (default: 42).
    - plaintext: Output plaintext only, not HDF5.
    """
    if not fastq_files or len(fastq_files) == 0:
        raise ValueError("At least one FASTQ file must be provided in fastq_files.")
    for f in fastq_files:
        if not f.exists():
            raise FileNotFoundError(f"FASTQ file not found: {f}")

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    if index is None and txnames is None:
        raise ValueError("Either index or txnames must be provided.")

    if index is not None and not index.exists():
        raise FileNotFoundError(f"Index file not found: {index}")

    if txnames is not None and not txnames.exists():
        raise FileNotFoundError(f"txnames file not found: {txnames}")

    if ec_file is not None and not ec_file.exists():
        raise FileNotFoundError(f"ec_file not found: {ec_file}")

    if fragment_file is not None and not fragment_file.exists():
        raise FileNotFoundError(f"fragment_file not found: {fragment_file}")

    if genemap is not None and not genemap.exists():
        raise FileNotFoundError(f"genemap file not found: {genemap}")

    if gtf is not None and not gtf.exists():
        raise FileNotFoundError(f"gtf file not found: {gtf}")

    if bootstrap_samples < 0:
        raise ValueError("bootstrap_samples must be >= 0.")

    if seed < 0:
        raise ValueError("seed must be >= 0.")

    if threads < 1:
        raise ValueError("threads must be >= 1.")

    cmd = ["kallisto", "bus", "-o", str(output_dir), "-t", str(threads)]

    if index is not None:
        cmd += ["-i", str(index)]
    if txnames is not None:
        cmd += ["-T", str(txnames)]
    if ec_file is not None:
        cmd += ["-e", str(ec_file)]
    if fragment_file is not None:
        cmd += ["-f", str(fragment_file)]
    if long:
        cmd.append("--long")
    if platform is not None:
        if platform not in ["PacBio", "ONT"]:
            raise ValueError("platform must be 'PacBio' or 'ONT' if specified.")
        cmd += ["-p", platform]
    if fragment_length is not None:
        if fragment_length <= 0:
            raise ValueError("fragment_length must be > 0 if specified.")
        cmd += ["-l", str(fragment_length)]
    if sd is not None:
        if sd <= 0:
            raise ValueError("sd must be > 0 if specified.")
        cmd += ["-s", str(sd)]
    if genemap is not None:
        cmd += ["-g", str(genemap)]
    if gtf is not None:
        cmd += ["-G", str(gtf)]
    if bootstrap_samples != 0:
        cmd += ["-b", str(bootstrap_samples)]
    if matrix_to_files:
        cmd.append("--matrix-to-files")
    if matrix_to_directories:
        cmd.append("--matrix-to-directories")
    if seed != 42:
        cmd += ["--seed", str(seed)]
    if plaintext:
        cmd.append("--plaintext")

    # Add fastq files at the end
    cmd += [str(f) for f in fastq_files]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Output files: output_dir contains output.bus, matrix.ec, transcripts.txt, etc.
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output_dir)],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto bus failed with exit code {e.returncode}",
        }


@mcp.tool()
def h5dump(
    abundance_h5: Path,
    output_dir: Path,
):
    """
    Converts HDF5-formatted results to plaintext.

    Parameters:
    - abundance_h5: Path to the abundance.h5 file.
    - output_dir: Directory to write output to.
    """
    if not abundance_h5.exists():
        raise FileNotFoundError(f"abundance.h5 file not found: {abundance_h5}")

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["kallisto", "h5dump", "-o", str(output_dir), str(abundance_h5)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Output files are plaintext abundance files in output_dir
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output_dir)],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto h5dump failed with exit code {e.returncode}",
        }


@mcp.tool()
def inspect(
    index_file: Path,
    threads: int = 1,
):
    """
    Inspects and gives information about a kallisto index.

    Parameters:
    - index_file: Path to the kallisto index file.
    - threads: Number of threads to use (default: 1).
    """
    if not index_file.exists():
        raise FileNotFoundError(f"Index file not found: {index_file}")

    if threads < 1:
        raise ValueError("threads must be >= 1.")

    cmd = ["kallisto", "inspect", str(index_file)]
    if threads != 1:
        cmd += ["-t", str(threads)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Output is printed to stdout
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto inspect failed with exit code {e.returncode}",
        }


@mcp.tool()
def version():
    """
    Prints kallisto version information.
    """
    cmd = ["kallisto", "version"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout.strip(),
            "stderr": result.stderr,
            "output_files": [],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto version failed with exit code {e.returncode}",
        }


@mcp.tool()
def cite():
    """
    Prints kallisto citation information.
    """
    cmd = ["kallisto", "cite"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout.strip(),
            "stderr": result.stderr,
            "output_files": [],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"kallisto cite failed with exit code {e.returncode}",
        }


if __name__ == '__main__':
    mcp.run()