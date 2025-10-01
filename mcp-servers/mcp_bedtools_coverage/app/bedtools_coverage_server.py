from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional
import subprocess

mcp = FastMCP()


@mcp.tool()
def bedtools_coverage(
    a: Path,
    b: List[Path],
    abam: bool = False,
    hist: bool = False,
    d: bool = False,
    counts: bool = False,
    f: float = 1e-9,
    F: float = 1e-9,
    r: bool = False,
    e: bool = False,
    s: bool = False,
    S: bool = False,
    split: bool = False,
    sorted: bool = False,
    g: Optional[Path] = None,
    header: bool = False,
    sortout: bool = False,
    nobuf: bool = False,
    iobuf: Optional[str] = None,
) -> dict:
    """
    Compute depth and breadth of coverage of features in file B on features in file A using bedtools coverage.

    Parameters:
    - a: Path to file A (BAM/BED/GFF/VCF). Features in A are compared to B.
    - b: List of one or more paths to file(s) B (BAM/BED/GFF/VCF). Can include multiple files.
    - abam: Treat file A as BAM input.
    - hist: Report histogram of coverage for each feature in A and summary histogram.
    - d: Report depth at each position in each A feature (one-based positions).
    - counts: Only report count of overlaps, no fraction computations.
    - f: Minimum overlap required as fraction of A (default 1e-9).
    - F: Minimum overlap required as fraction of B (default 1e-9).
    - r: Require reciprocal fraction overlap for A and B.
    - e: Require minimum fraction satisfied for A OR B (instead of both).
    - s: Force strandedness; only report hits overlapping on same strand.
    - S: Require different strandedness; only report hits overlapping on opposite strand.
    - split: Treat split BAM or BED12 entries as distinct intervals.
    - sorted: Use memory-efficient sweeping algorithm; requires position-sorted input.
    - g: Genome file defining chromosome order (used with -sorted).
    - header: Print header from A file prior to results.
    - sortout: When multiple databases (-b), sort output DB hits for each record.
    - nobuf: Disable buffered output; print lines as generated.
    - iobuf: Integer size of read buffer (e.g. 4K, 10M). No effect with compressed files.

    Returns:
    A dict with keys:
    - command_executed: The full command line executed.
    - stdout: Standard output from bedtools coverage.
    - stderr: Standard error from bedtools coverage.
    - output_files: Empty list (output is stdout).
    """
    # Validate input files
    if not a.exists():
        raise FileNotFoundError(f"Input file A not found: {a}")
    for bf in b:
        if not bf.exists():
            raise FileNotFoundError(f"Input file B not found: {bf}")

    # Validate float parameters
    if not (0.0 <= f <= 1.0):
        raise ValueError(f"Parameter f must be between 0.0 and 1.0, got {f}")
    if not (0.0 <= F <= 1.0):
        raise ValueError(f"Parameter F must be between 0.0 and 1.0, got {F}")

    # Validate iobuf if provided
    if iobuf is not None:
        valid_suffixes = ('K', 'M', 'G')
        if len(iobuf) < 2 or not iobuf[:-1].isdigit() or iobuf[-1].upper() not in valid_suffixes:
            raise ValueError(f"iobuf must be integer followed by K/M/G suffix, got {iobuf}")

    # Build command
    cmd = ["bedtools", "coverage"]

    # -a parameter
    if abam:
        cmd.append("-abam")
    else:
        cmd.append("-a")
    cmd.append(str(a))

    # -b parameter(s)
    # Join multiple files with comma
    b_joined = ",".join(str(x) for x in b)
    cmd.extend(["-b", b_joined])

    # Optional flags
    if hist:
        cmd.append("-hist")
    if d:
        cmd.append("-d")
    if counts:
        cmd.append("-counts")
    if r:
        cmd.append("-r")
    if e:
        cmd.append("-e")
    if s:
        cmd.append("-s")
    if S:
        cmd.append("-S")
    if split:
        cmd.append("-split")
    if sorted:
        cmd.append("-sorted")
    if header:
        cmd.append("-header")
    if sortout:
        cmd.append("-sortout")
    if nobuf:
        cmd.append("-nobuf")
    if g is not None:
        if not g.exists():
            raise FileNotFoundError(f"Genome file g not found: {g}")
        cmd.extend(["-g", str(g)])

    # Parameters with values
    cmd.extend(["-f", str(f)])
    cmd.extend(["-F", str(F)])

    if iobuf is not None:
        cmd.extend(["-iobuf", iobuf])

    # Run subprocess
    try:
        completed = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = completed.stdout
        stderr = completed.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"bedtools coverage failed with exit code {e.returncode}",
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [],
    }


if __name__ == '__main__':
    mcp.run()