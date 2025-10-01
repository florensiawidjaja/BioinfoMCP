from fastmcp import FastMCP
from pathlib import Path
import subprocess

mcp = FastMCP()

@mcp.tool()
def computeGCBias(
    bamfile: str,
    effectiveGenomeSize: int,
    genome: str,
    fragmentLength: int = 200,
    GCbiasFrequenciesFile: str = "",
):
    """
    Compute GC bias from a BAM file.

    Parameters:
    - bamfile: Path to input BAM file (required).
    - effectiveGenomeSize: Effective genome size (required).
    - genome: Genome file in 2bit format (required).
    - fragmentLength: Fragment length, default 200.
    - GCbiasFrequenciesFile: Output file for GC bias frequencies (required).

    Returns:
    - command_executed: The full command line executed.
    - stdout: Standard output from the tool.
    - stderr: Standard error from the tool.
    - output_files: List containing the GC bias frequencies output file.
    """
    # Validate input files
    bam_path = Path(bamfile)
    if not bam_path.is_file():
        raise FileNotFoundError(f"BAM file not found: {bamfile}")

    genome_path = Path(genome)
    if not genome_path.is_file():
        raise FileNotFoundError(f"Genome file not found: {genome}")

    # Validate output file path
    freq_path = Path(GCbiasFrequenciesFile)
    if freq_path.exists() and not freq_path.is_file():
        raise ValueError(f"GCbiasFrequenciesFile path exists and is not a file: {GCbiasFrequenciesFile}")

    if effectiveGenomeSize <= 0:
        raise ValueError("effectiveGenomeSize must be a positive integer")

    if fragmentLength <= 0:
        raise ValueError("fragmentLength must be a positive integer")

    cmd = [
        "computeGCBias",
        "-b", str(bam_path),
        "--effectiveGenomeSize", str(effectiveGenomeSize),
        "-g", str(genome_path),
        "-l", str(fragmentLength),
        "--GCbiasFrequenciesFile", str(freq_path)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(freq_path)]
    }


if __name__ == '__main__':
    mcp.run()