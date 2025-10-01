from fastmcp import FastMCP
from pathlib import Path
import subprocess
from typing import Optional, List

mcp = FastMCP()

@mcp.tool()
def fastqc(
    input_files: List[Path],
    outdir: Path = Path("fastqc_output"),
    threads: int = 1,
    extract: bool = False,
    quiet: bool = False,
    noextract: bool = False,
    format: Optional[str] = None,
    contamination: Optional[Path] = None,
):
    """
    Run FastQC quality control checks on raw sequence data files.

    Parameters:
    - input_files: List of input sequence files (fastq or bam).
    - outdir: Output directory for FastQC results (default: fastqc_output).
    - threads: Number of threads to use (default: 1).
    - extract: If True, extract the zipped output files.
    - quiet: If True, suppress all FastQC output.
    - noextract: If True, do not extract zipped output files.
    - format: Force input file format (fastq or bam).
    - contamination: Path to contamination file for contamination check.

    Returns:
    - command_executed: The full command line executed.
    - stdout: Standard output from FastQC.
    - stderr: Standard error from FastQC.
    - output_files: List of generated output files/directories.
    """
    # Validate input files
    if not input_files:
        raise ValueError("At least one input file must be provided.")
    #for f in input_files:
    #    if not f.exists():
    #        raise FileNotFoundError(f"Input file does not exist: {f}")

    # Validate threads
    if threads < 1:
        raise ValueError("threads must be >= 1")

    # Validate format
    if format is not None and format not in ("fastq", "bam"):
        raise ValueError("format must be either 'fastq' or 'bam' if specified")

    # Validate contamination file
    if contamination is not None and not contamination.exists():
        raise FileNotFoundError(f"Contamination file does not exist: {contamination}")

    # Prepare output directory
    outdir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = ["fastqc"]
    cmd += ["--outdir", str(outdir)]
    cmd += ["--threads", str(threads)]
    if extract:
        cmd.append("--extract")
    if quiet:
        cmd.append("--quiet")
    if noextract:
        cmd.append("--noextract")
    if format:
        cmd += ["--format", format]
    if contamination:
        cmd += ["--contamination", str(contamination)]
    cmd += [str(f) for f in input_files]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"FastQC failed with return code {e.returncode}"
        }

    # Collect output files: FastQC creates .html and .zip files per input file in outdir
    output_files = []
    #for f in input_files:
    #    base = f.stem
    #    html_file = outdir / f"{base}_fastqc.html"
    #    zip_file = outdir / f"{base}_fastqc.zip"
    #    if html_file.exists():
    #        output_files.append(str(html_file))
    #    if zip_file.exists():
    #        output_files.append(str(zip_file))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }

if __name__ == '__main__':
    mcp.run()