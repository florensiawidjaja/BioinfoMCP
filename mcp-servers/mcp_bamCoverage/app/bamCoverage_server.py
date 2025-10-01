from fastmcp import FastMCP
from pathlib import Path
import subprocess

mcp = FastMCP()

@mcp.tool()
def bamCoverage(
    bam: Path,
    outFileName: Path,
    binSize: int = 50,
    normalizeUsing: str = "None",
    ignoreDuplicates: bool = False,
    minMappingQuality: int = 0,
    numberOfProcessors: int = 1,
) -> dict:
    """
    Generate a coverage file (bigWig) from a BAM file using bamCoverage.

    Parameters:
    - bam: Input BAM file path (required).
    - outFileName: Output coverage file path (required).
    - binSize: Size of the bins (default 50).
    - normalizeUsing: Normalization method (default "None").
    - ignoreDuplicates: Whether to ignore duplicates (default False).
    - minMappingQuality: Minimum mapping quality to consider (default 0).
    - numberOfProcessors: Number of processors to use (default 1).
    """
    # Validate input files
    if not bam.exists() or not bam.is_file():
        raise FileNotFoundError(f"Input BAM file does not exist: {bam}")
    if binSize <= 0:
        raise ValueError("binSize must be a positive integer")
    if minMappingQuality < 0:
        raise ValueError("minMappingQuality must be >= 0")
    if numberOfProcessors <= 0:
        raise ValueError("numberOfProcessors must be a positive integer")

    # Build command
    cmd = [
        "bamCoverage",
        "-b", str(bam),
        "-o", str(outFileName),
        "--binSize", str(binSize),
        "--minMappingQuality", str(minMappingQuality),
        "--numberOfProcessors", str(numberOfProcessors),
    ]
    if normalizeUsing != "None":
        cmd.extend(["--normalizeUsing", normalizeUsing])
    if ignoreDuplicates:
        cmd.append("--ignoreDuplicates")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bamCoverage failed with exit code {e.returncode}"
        }

    output_files = []
    if outFileName.exists():
        output_files.append(str(outFileName))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()