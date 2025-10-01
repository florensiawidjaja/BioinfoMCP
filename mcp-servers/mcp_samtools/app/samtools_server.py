from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()

@mcp.tool()
def samtools_view(
    input_file: Path,
    regions: Optional[List[str]] = None,
    output_file: Optional[Path] = None,
    output_format: str = "sam",
    write_index: bool = False,
    threads: int = 1,
    filter_expr: Optional[str] = None,
    reference: Optional[Path] = None,
    input_fmt: Optional[str] = None,
    input_fmt_option: Optional[str] = None,
    output_fmt_option: Optional[str] = None,
    verbosity: int = 3,
) -> dict:
    """
    samtools view: Extract and convert alignments from SAM/BAM/CRAM files.
    Prints all alignments or those overlapping specified regions.
    Supports format conversion and indexing.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")
    if reference is not None and not reference.exists():
        raise FileNotFoundError(f"Reference file {reference} does not exist")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if output_format not in ("sam", "bam", "cram"):
        raise ValueError("output_format must be one of 'sam', 'bam', or 'cram'")

    cmd = ["samtools", "view"]
    if input_fmt:
        cmd.extend(["--input-fmt", input_fmt])
    if input_fmt_option:
        cmd.extend(["--input-fmt-option", input_fmt_option])
    cmd.extend(["-O", output_format])
    if output_fmt_option:
        cmd.extend(["--output-fmt-option", output_fmt_option])
    if write_index:
        cmd.append("--write-index")
    if verbosity != 3:
        cmd.extend(["--verbosity", str(verbosity)])
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    if filter_expr:
        cmd.extend(["--filter", filter_expr])
    if reference:
        cmd.extend(["--reference", str(reference)])

    if output_file:
        cmd.extend(["-o", str(output_file)])

    cmd.append(str(input_file))
    if regions:
        cmd.extend(regions)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools view failed with exit code {e.returncode}",
            "output_files": []
        }

    output_files = [str(output_file)] if output_file else []

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


@mcp.tool()
def samtools_sort(
    input_file: Path,
    output_file: Path,
    level: int = 0,
    max_mem: Optional[str] = None,
    output_format: str = "bam",
    sort_by_name: bool = False,
    tag: Optional[str] = None,
    tmp_prefix: Optional[Path] = None,
    threads: int = 1,
) -> dict:
    """
    samtools sort: Sort alignments by leftmost coordinates or by read name.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")
    if output_file.exists():
        raise FileExistsError(f"Output file {output_file} already exists")
    if not (0 <= level <= 9):
        raise ValueError("level must be between 0 and 9")
    if output_format not in ("sam", "bam", "cram"):
        raise ValueError("output_format must be one of 'sam', 'bam', or 'cram'")
    if threads < 1:
        raise ValueError("threads must be >= 1")

    cmd = ["samtools", "sort"]
    if level > 0:
        cmd.extend(["-l", str(level)])
    if max_mem:
        cmd.extend(["-m", max_mem])
    cmd.extend(["-o", str(output_file)])
    cmd.extend(["-O", output_format])
    if sort_by_name:
        cmd.append("-n")
    if tag:
        cmd.extend(["-t", tag])
    if tmp_prefix:
        cmd.extend(["-T", str(tmp_prefix)])
    if threads > 1:
        cmd.extend(["-@", str(threads)])

    cmd.append(str(input_file))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools sort failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(output_file)]
    }


@mcp.tool()
def samtools_index(
    input_file: Path,
    output_index: Optional[Path] = None,
    bai: bool = False,
    csi: bool = False,
    min_shift: int = 14,
    threads: int = 1,
) -> dict:
    """
    samtools index: Index a coordinate-sorted SAM, BAM or CRAM file for fast random access.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")
    if output_index and output_index.exists():
        raise FileExistsError(f"Output index file {output_index} already exists")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if bai and csi:
        raise ValueError("Cannot specify both bai and csi")

    cmd = ["samtools", "index"]
    if bai:
        cmd.append("-b")
    if csi:
        cmd.append("-c")
    if min_shift != 14:
        cmd.extend(["-m", str(min_shift)])
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(input_file))
    if output_index:
        cmd.append(str(output_index))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools index failed with exit code {e.returncode}",
            "output_files": []
        }

    output_files = [str(output_index)] if output_index else []
    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


@mcp.tool()
def samtools_flagstat(
    input_file: Path,
) -> dict:
    """
    samtools flagstat: Compute statistics from the input SAM/BAM/CRAM file.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")

    cmd = ["samtools", "flagstat", str(input_file)]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools flagstat failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": []
    }


@mcp.tool()
def samtools_merge(
    output_file: Path,
    input_files: List[Path],
    no_rg: bool = False,
    update_header: Optional[Path] = None,
    threads: int = 1,
) -> dict:
    """
    samtools merge: Merge multiple sorted alignment files into one sorted output file.
    """
    if not input_files:
        raise ValueError("At least one input file must be specified")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file {f} does not exist")
    if output_file.exists():
        raise FileExistsError(f"Output file {output_file} already exists")
    if threads < 1:
        raise ValueError("threads must be >= 1")

    cmd = ["samtools", "merge"]
    if no_rg:
        cmd.append("-n")
    if update_header:
        if not update_header.exists():
            raise FileNotFoundError(f"Header file {update_header} does not exist")
        cmd.extend(["-h", str(update_header)])
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(output_file))
    cmd.extend(str(f) for f in input_files)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools merge failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(output_file)]
    }


@mcp.tool()
def samtools_faidx(
    fasta_file: Path,
    regions: Optional[List[str]] = None,
) -> dict:
    """
    samtools faidx: Index a FASTA file or extract subsequences from indexed FASTA.
    """
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file {fasta_file} does not exist")

    cmd = ["samtools", "faidx", str(fasta_file)]
    if regions:
        cmd.extend(regions)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools faidx failed with exit code {e.returncode}",
            "output_files": []
        }

    output_files = []
    # If no regions specified, faidx creates fasta_file.fai index file
    if not regions:
        index_file = fasta_file.with_suffix(fasta_file.suffix + ".fai")
        if index_file.exists():
            output_files.append(str(index_file))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


@mcp.tool()
def samtools_fastq(
    input_file: Path,
    output_file: Optional[Path] = None,
    soft_clip: bool = False,
    threads: int = 1,
) -> dict:
    """
    samtools fastq: Convert BAM/CRAM to FASTQ format.
    Input must be collated by name.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")
    if output_file and output_file.exists():
        raise FileExistsError(f"Output file {output_file} already exists")
    if threads < 1:
        raise ValueError("threads must be >= 1")

    cmd = ["samtools", "fastq"]
    if soft_clip:
        cmd.append("--soft-clipped")
    if threads > 1:
        cmd.extend(["-@", str(threads)])
    cmd.append(str(input_file))
    if output_file:
        cmd.extend(["-o", str(output_file)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools fastq failed with exit code {e.returncode}",
            "output_files": []
        }

    output_files = [str(output_file)] if output_file else []

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


@mcp.tool()
def samtools_flag_convert(
    flags: str,
) -> dict:
    """
    samtools flags: Convert between textual and numeric flag representation.
    Input is a comma-separated string of flags or a numeric flag.
    """
    if not flags:
        raise ValueError("flags parameter must be provided")

    cmd = ["samtools", "flags", flags]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout.strip()
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools flags failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": []
    }


@mcp.tool()
def samtools_quickcheck(
    input_files: List[Path],
    verbose: bool = False,
) -> dict:
    """
    samtools quickcheck: Quickly check that input files appear intact.
    """
    if not input_files:
        raise ValueError("At least one input file must be specified")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file {f} does not exist")

    cmd = ["samtools", "quickcheck"]
    if verbose:
        cmd.append("-v")
    cmd.extend(str(f) for f in input_files)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        # quickcheck returns non-zero if files are corrupted
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools quickcheck failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": []
    }


@mcp.tool()
def samtools_stats(
    input_file: Path,
    regions: Optional[List[str]] = None,
) -> dict:
    """
    samtools stats: Collect statistics from BAM/CRAM files and output in text format.
    """
    if not input_file.exists():
        raise FileNotFoundError(f"Input file {input_file} does not exist")

    cmd = ["samtools", "stats", str(input_file)]
    if regions:
        cmd.extend(regions)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools stats failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": []
    }


@mcp.tool()
def samtools_depth(
    input_files: List[Path],
    regions: Optional[List[str]] = None,
) -> dict:
    """
    samtools depth: Compute read depth at each position or region.
    """
    if not input_files:
        raise ValueError("At least one input file must be specified")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file {f} does not exist")

    cmd = ["samtools", "depth"]
    for f in input_files:
        cmd.append(str(f))
    if regions:
        cmd.extend(regions)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"samtools depth failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": []
    }


if __name__ == '__main__':
    mcp.run()