from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List, Union
import subprocess

mcp = FastMCP()


@mcp.tool()
def minimap_index(
    target_fa: Path,
    output_index: Path,
    preset: Optional[str] = None,
    homopolymer_compressed: bool = False,
    kmer_length: int = 15,
    window_size: int = 10,
    syncmer_size: int = 10,
    max_target_bases: str = "8G",
    idx_no_seq: bool = False,
    alt_file: Optional[Path] = None,
    alt_drop_fraction: float = 0.15,
):
    """
    Create a minimizer index from target sequences.

    Parameters:
    - target_fa: Path to the target FASTA file.
    - output_index: Path to save the minimizer index (.mmi).
    - preset: Optional preset string to apply indexing presets.
    - homopolymer_compressed: Use homopolymer-compressed minimizers.
    - kmer_length: Minimizer k-mer length (default 15).
    - window_size: Minimizer window size (default 10).
    - syncmer_size: Syncmer submer size (default 10).
    - max_target_bases: Max target bases loaded into RAM for indexing (default "8G").
    - idx_no_seq: Do not store target sequences in the index.
    - alt_file: Optional path to ALT contigs list file.
    - alt_drop_fraction: Drop ALT hits by this fraction when ranking (default 0.15).
    """
    # Validate input files
    if not target_fa.is_file():
        raise FileNotFoundError(f"Target FASTA file not found: {target_fa}")
    if alt_file is not None and not alt_file.is_file():
        raise FileNotFoundError(f"ALT contigs file not found: {alt_file}")

    # Validate numeric parameters
    if kmer_length < 1:
        raise ValueError("kmer_length must be positive integer")
    if window_size < 1:
        raise ValueError("window_size must be positive integer")
    if syncmer_size < 1:
        raise ValueError("syncmer_size must be positive integer")
    if not (0.0 <= alt_drop_fraction <= 1.0):
        raise ValueError("alt_drop_fraction must be between 0 and 1")

    cmd = ["minimap2"]
    if preset:
        cmd.extend(["-x", preset])
    if homopolymer_compressed:
        cmd.append("-H")
    cmd.extend(["-k", str(kmer_length)])
    cmd.extend(["-w", str(window_size)])
    cmd.extend(["-j", str(syncmer_size)])
    cmd.extend(["-I", max_target_bases])
    if idx_no_seq:
        cmd.append("--idx-no-seq")
    cmd.extend(["-d", str(output_index)])
    if alt_file:
        cmd.extend(["--alt", str(alt_file)])
        cmd.extend(["--alt-drop", str(alt_drop_fraction)])
    cmd.append(str(target_fa))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Indexing failed with return code {e.returncode}"
        }

    output_files = []
    if output_index.is_file():
        output_files.append(str(output_index))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def minimap_map(
    target: Union[Path, str],
    query: Path,
    output: Optional[Path] = None,
    sam_output: bool = False,
    preset: Optional[str] = None,
    threads: int = 3,
    no_secondary: bool = False,
    max_query_length: Optional[int] = None,
    cs_tag: Optional[str] = None,  # None means no cs tag, "short" or "long"
    md_tag: bool = False,
    eqx_cigar: bool = False,
    soft_clip_supplementary: bool = False,
    secondary_seq: bool = False,
    seed: int = 11,
    io_threads_2: bool = False,
    max_bases_batch: str = "500M",
    paf_no_hit: bool = False,
    sam_hit_only: bool = False,
    read_group: Optional[str] = None,
    copy_comments: bool = False,
):
    """
    Map query sequences to target sequences or index.

    Parameters:
    - target: Path to target FASTA or minimap2 index (.mmi) file.
    - query: Path to query FASTA/FASTQ file.
    - output: Optional output file path. If None, output to stdout.
    - sam_output: Output SAM format with CIGAR (-a).
    - preset: Optional preset string to apply mapping presets.
    - threads: Number of threads to use (default 3).
    - no_secondary: Disable secondary alignments output.
    - max_query_length: Filter out query sequences longer than this length.
    - cs_tag: Output cs tag; None=no, "short" or "long".
    - md_tag: Output MD tag.
    - eqx_cigar: Output =/X CIGAR operators.
    - soft_clip_supplementary: Use soft clipping for supplementary alignments (-Y).
    - secondary_seq: Show query sequences for secondary alignments.
    - seed: Integer seed for randomizing equally best hits (default 11).
    - io_threads_2: Use two I/O threads during mapping (-2).
    - max_bases_batch: Number of bases loaded into memory per mini-batch (default "500M").
    - paf_no_hit: In PAF, output unmapped queries.
    - sam_hit_only: In SAM, do not output unmapped reads.
    - read_group: SAM read group line string (e.g. '@RG\tID:foo\tSM:bar').
    - copy_comments: Copy input FASTA/Q comments to output (-y).
    """
    # Validate input files
    if isinstance(target, Path):
        if not target.is_file():
            raise FileNotFoundError(f"Target file not found: {target}")
    else:
        # target is str, could be index or fasta, no file check here
        pass
    if not query.is_file():
        raise FileNotFoundError(f"Query file not found: {query}")

    if threads < 1:
        raise ValueError("threads must be positive integer")
    if max_query_length is not None and max_query_length < 1:
        raise ValueError("max_query_length must be positive integer if set")
    if seed < 0:
        raise ValueError("seed must be non-negative integer")

    cmd = ["minimap2"]
    if preset:
        cmd.extend(["-x", preset])
    if sam_output:
        cmd.append("-a")
    if no_secondary:
        cmd.append("--secondary=no")
    else:
        cmd.append("--secondary=yes")
    if max_query_length is not None:
        cmd.extend(["--max-qlen", str(max_query_length)])
    if cs_tag is not None:
        if cs_tag not in ("short", "long"):
            raise ValueError("cs_tag must be 'short', 'long', or None")
        if cs_tag == "short":
            cmd.append("--cs")
        else:
            cmd.append("--cs=long")
    if md_tag:
        cmd.append("--MD")
    if eqx_cigar:
        cmd.append("--eqx")
    if soft_clip_supplementary:
        cmd.append("-Y")
    if secondary_seq:
        cmd.append("--secondary-seq")
    cmd.extend(["-t", str(threads)])
    if io_threads_2:
        cmd.append("-2")
    cmd.extend(["-K", max_bases_batch])
    cmd.extend(["-s", str(seed)])
    if paf_no_hit:
        cmd.append("--paf-no-hit")
    if sam_hit_only:
        cmd.append("--sam-hit-only")
    if read_group:
        cmd.extend(["-R", read_group])
    if copy_comments:
        cmd.append("-y")

    # target can be index or fasta
    cmd.append(str(target))
    cmd.append(str(query))

    # Output handling
    if output is not None:
        output_path = Path(output)
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        stdout_target = open(output_path, "w")
    else:
        stdout_target = subprocess.PIPE

    try:
        result = subprocess.run(
            cmd,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        if output is None:
            stdout = result.stdout
        else:
            stdout = ""
    except subprocess.CalledProcessError as e:
        if output is not None:
            stdout_target.close()
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": [],
            "error": f"Mapping failed with return code {e.returncode}"
        }
    if output is not None:
        stdout_target.close()

    output_files = []
    if output is not None and Path(output).is_file():
        output_files.append(str(output))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": result.stderr if output is None else "",
        "output_files": output_files
    }


@mcp.tool()
def minimap_version():
    """
    Get minimap2 version string.
    """
    cmd = ["minimap2", "--version"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Failed to get version with return code {e.returncode}"
        }
    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout.strip(),
        "stderr": result.stderr,
        "output_files": []
    }


if __name__ == '__main__':
    mcp.run()