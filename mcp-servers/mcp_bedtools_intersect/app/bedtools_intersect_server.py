from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional
import subprocess

mcp = FastMCP()

@mcp.tool()
def bedtools_intersect(
    a: Path,
    b: List[Path],
    abam: bool = False,
    ubam: bool = False,
    bed: bool = False,
    wa: bool = False,
    wb: bool = False,
    loj: bool = False,
    wo: bool = False,
    wao: bool = False,
    u: bool = False,
    c: bool = False,
    C: bool = False,
    v: bool = False,
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
    names: Optional[List[str]] = None,
    filenames: bool = False,
    sortout: bool = False,
    nobuf: bool = False,
    iobuf: Optional[str] = None,
) -> dict:
    """
    Run bedtools intersect to find overlaps between two sets of genomic features.

    Parameters:
    - a: Path to file A (BAM/BED/GFF/VCF). Required.
    - b: List of one or more files B (BAM/BED/GFF/VCF). Required.
    - abam: Treat A as BAM input (deprecated, autodetected).
    - ubam: Write uncompressed BAM output.
    - bed: Write output as BED when using BAM input.
    - wa: Write original entry in A for each overlap.
    - wb: Write original entry in B for each overlap.
    - loj: Left outer join; report all A features with or without overlaps.
    - wo: Write original A and B entries plus number of base pairs of overlap.
    - wao: Like -wo but also report A features without overlap with overlap=0.
    - u: Write original A entry once if any overlaps found in B.
    - c: For each A entry, report number of hits in B.
    - C: For each A entry, report number of overlaps with each B file separately.
    - v: Only report A entries with no overlap in B.
    - f: Minimum overlap fraction of A (default 1e-9).
    - F: Minimum overlap fraction of B (default 1e-9).
    - r: Require reciprocal overlap fraction for A and B.
    - e: Require minimum fraction satisfied for A OR B.
    - s: Force strandedness (overlaps on same strand only).
    - S: Require different strandedness (overlaps on opposite strand only).
    - split: Treat split BAM or BED12 entries as distinct intervals.
    - sorted: Use memory-efficient algorithm for sorted input.
    - g: Genome file defining chromosome order (used with -sorted).
    - header: Print header from A file before results.
    - names: Aliases for multiple -b files.
    - filenames: Show complete filename instead of fileId for multiple -b files.
    - sortout: Sort output DB hits for each record (multiple -b files).
    - nobuf: Disable buffered output.
    - iobuf: Integer size of read buffer (K/M/G suffix supported).

    Returns:
    - dict with keys: command_executed, stdout, stderr, output_files (empty list).
    """
    # Validate input files
    if not a.exists():
        raise FileNotFoundError(f"Input file -a does not exist: {a}")
    for bf in b:
        if not bf.exists():
            raise FileNotFoundError(f"Input file in -b does not exist: {bf}")

    # Validate float parameters
    if not (0.0 <= f <= 1.0):
        raise ValueError(f"Parameter -f must be between 0.0 and 1.0, got {f}")
    if not (0.0 <= F <= 1.0):
        raise ValueError(f"Parameter -F must be between 0.0 and 1.0, got {F}")

    # Validate iobuf if given
    if iobuf is not None:
        if not isinstance(iobuf, str) or not iobuf[:-1].isdigit() and not iobuf.isdigit():
            raise ValueError(f"Parameter -iobuf must be an integer optionally followed by K/M/G suffix, got {iobuf}")

    # Validate names length if given
    if names is not None:
        if len(names) != len(b):
            raise ValueError(f"Length of -names ({len(names)}) must match number of -b files ({len(b)})")

    # Validate genome file if given
    if g is not None and not g.exists():
        raise FileNotFoundError(f"Genome file -g does not exist: {g}")

    # Build command
    cmd = ["bedtools", "intersect"]

    # Input A
    if abam:
        cmd += ["-abam", str(a)]
    else:
        cmd += ["-a", str(a)]

    # Input B files
    for bf in b:
        cmd.append("-b")
        cmd.append(str(bf))

    # Output options
    if ubam:
        cmd.append("-ubam")
    if bed:
        cmd.append("-bed")
    if wa:
        cmd.append("-wa")
    if wb:
        cmd.append("-wb")
    if loj:
        cmd.append("-loj")
    if wo:
        cmd.append("-wo")
    if wao:
        cmd.append("-wao")
    if u:
        cmd.append("-u")
    if c:
        cmd.append("-c")
    if C:
        cmd.append("-C")
    if v:
        cmd.append("-v")

    # Overlap fractions and flags
    cmd += ["-f", str(f)]
    cmd += ["-F", str(F)]
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
    if g is not None:
        cmd += ["-g", str(g)]
    if header:
        cmd.append("-header")
    if names is not None:
        cmd += ["-names"] + names
    if filenames:
        cmd.append("-filenames")
    if sortout:
        cmd.append("-sortout")
    if nobuf:
        cmd.append("-nobuf")
    if iobuf is not None:
        cmd += ["-iobuf", iobuf]

    # Run subprocess
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else f"bedtools intersect failed with exit code {e.returncode}",
            "output_files": [],
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": [],
    }


if __name__ == '__main__':
    mcp.run()