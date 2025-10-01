from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Union
import subprocess

mcp = FastMCP()


@mcp.tool()
def macs3_callpeak(
    treatment: List[Path],
    control: Optional[List[Path]] = None,
    name: str = "macs3_callpeak",
    format: str = "AUTO",
    outdir: Optional[Path] = None,
    bdg: bool = False,
    trackline: bool = False,
    gsize: str = "hs",
    tsize: int = 0,
    qvalue: float = 0.05,
    pvalue: float = 0.0,
    min_length: int = 0,
    max_gap: int = 0,
    nolambda: bool = False,
    slocal: int = 1000,
    llocal: int = 10000,
    nomodel: bool = False,
    extsize: int = 0,
    shift: int = 0,
    keep_dup: Union[str, int] = 1,
    broad: bool = False,
    broad_cutoff: float = 0.1,
    scale_to: str = "small",
    call_summits: bool = False,
    buffer_size: int = 100000,
    cutoff_analysis: bool = False,
    barcodes: Optional[Path] = None,
    max_count: Optional[int] = None,
):
    """
    Call significantly enriched regions (peaks) from alignment files using MACS3 callpeak.

    Parameters:
    - treatment (-t, --treatment): List of treatment alignment files (required).
    - control (-c, --control): List of control alignment files (optional).
    - name (-n, --name): Name string for experiment, used as prefix for output files.
    - format (-f, --format): Format of tag files. One of ELAND, BED, ELANDMULTI, ELANDEXPORT,
      SAM, BAM, BOWTIE, BAMPE, BEDPE, FRAG, or AUTO (default).
    - outdir (--outdir): Directory to save output files. Created if doesn't exist.
    - bdg (-B, --bdg): If True, output bedGraph files for fragment pileup and control lambda.
    - trackline (--trackline): If True, include UCSC genome browser trackline in output headers.
    - gsize (-g, --gsize): Effective genome size. Precompiled values: hs, mm, ce, dm or numeric string.
    - tsize (-s, --tsize): Size of sequencing tags. 0 means auto-detect.
    - qvalue (-q, --qvalue): q-value cutoff for significant peaks (default 0.05).
    - pvalue (-p, --pvalue): p-value cutoff. If >0, p-value cutoff used instead of q-value.
    - min_length (--min-length): Minimum length of called peak. Default 0 means use fragment size.
    - max_gap (--max-gap): Maximum gap between nearby regions to merge. Default 0 means use read length.
    - nolambda (--nolambda): If True, use background lambda as local lambda (no local bias correction).
    - slocal (--slocal): Small local region size in bp for local lambda calculation (default 1000).
    - llocal (--llocal): Large local region size in bp for local lambda calculation (default 10000).
    - nomodel (--nomodel): If True, bypass building shifting model.
    - extsize (--extsize): When nomodel is set, extend reads to this fixed fragment size.
    - shift (--shift): Shift cutting ends by this bp. Must be 0 if format is BAMPE or BEDPE.
    - keep_dup (--keep-dup): How to handle duplicate tags: 'auto', 'all', or integer number (default 1).
    - broad (--broad): If True, perform broad peak calling producing gappedPeak format.
    - broad_cutoff (--broad-cutoff): Cutoff for broad regions (default 0.1). Requires --broad.
    - scale_to (--scale-to): 'large' or 'small' to scale dataset depths (default 'small').
    - call_summits (--call-summits): If True, reanalyze signal profile to call subpeak summits.
    - buffer_size (--buffer-size): Buffer size for internal array (default 100000).
    - cutoff_analysis (--cutoff-analysis): If True, perform cutoff analysis and output report.
    - barcodes (--barcodes): Barcode list file, only valid if format is FRAG.
    - max_count (--max-count): Max count per fragment, only valid if format is FRAG.

    Returns:
    Dict with keys: command_executed, stdout, stderr, output_files.
    """
    # Validate input files
    if not treatment or len(treatment) == 0:
        raise ValueError("At least one treatment file must be specified in 'treatment' parameter.")
    for f in treatment:
        if not f.exists():
            raise FileNotFoundError(f"Treatment file not found: {f}")
    if control:
        for f in control:
            if not f.exists():
                raise FileNotFoundError(f"Control file not found: {f}")

    # Validate format
    valid_formats = {
        "ELAND", "BED", "ELANDMULTI", "ELANDEXPORT", "SAM", "BAM", "BOWTIE",
        "BAMPE", "BEDPE", "FRAG", "AUTO"
    }
    format_upper = format.upper()
    if format_upper not in valid_formats:
        raise ValueError(f"Invalid format '{format}'. Must be one of {valid_formats}.")

    # Validate keep_dup
    if isinstance(keep_dup, str):
        if keep_dup not in {"auto", "all"}:
            raise ValueError("keep_dup string value must be 'auto' or 'all'.")
    elif isinstance(keep_dup, int):
        if keep_dup < 0:
            raise ValueError("keep_dup integer value must be non-negative.")
    else:
        raise ValueError("keep_dup must be str ('auto','all') or non-negative int.")

    # Validate scale_to
    if scale_to not in {"large", "small"}:
        raise ValueError("scale_to must be 'large' or 'small'.")

    # Validate broad_cutoff only if broad is True
    if broad:
        if broad_cutoff <= 0 or broad_cutoff > 1:
            raise ValueError("broad_cutoff must be > 0 and <= 1 when broad is enabled.")
    else:
        if broad_cutoff != 0.1:
            raise ValueError("broad_cutoff option is only valid when broad is enabled.")

    # Validate shift for paired-end formats
    if format_upper in {"BAMPE", "BEDPE"} and shift != 0:
        raise ValueError("shift must be 0 when format is BAMPE or BEDPE.")

    # Validate tsize
    if tsize < 0:
        raise ValueError("tsize must be >= 0.")

    # Validate qvalue and pvalue
    if qvalue <= 0 or qvalue > 1:
        raise ValueError("qvalue must be > 0 and <= 1.")
    if pvalue < 0 or pvalue > 1:
        raise ValueError("pvalue must be >= 0 and <= 1.")

    # Validate min_length and max_gap
    if min_length < 0:
        raise ValueError("min_length must be >= 0.")
    if max_gap < 0:
        raise ValueError("max_gap must be >= 0.")

    # Validate slocal and llocal
    if slocal <= 0:
        raise ValueError("slocal must be > 0.")
    if llocal <= 0:
        raise ValueError("llocal must be > 0.")

    # Validate buffer_size
    if buffer_size <= 0:
        raise ValueError("buffer_size must be > 0.")

    # Validate max_count only if format is FRAG
    if max_count is not None:
        if format_upper != "FRAG":
            raise ValueError("--max-count is only valid when format is FRAG.")
        if max_count < 1:
            raise ValueError("max_count must be >= 1.")

    # Validate barcodes only if format is FRAG
    if barcodes is not None:
        if format_upper != "FRAG":
            raise ValueError("--barcodes option is only valid when format is FRAG.")
        if not barcodes.exists():
            raise FileNotFoundError(f"Barcode list file not found: {barcodes}")

    # Prepare output directory
    if outdir is not None:
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)
        outdir_str = str(outdir.resolve())
    else:
        outdir_str = None

    # Build command line
    cmd = ["macs3", "callpeak"]

    # Treatment files
    for f in treatment:
        cmd.extend(["-t", str(f.resolve())])

    # Control files
    if control:
        for f in control:
            cmd.extend(["-c", str(f.resolve())])

    # Name
    cmd.extend(["-n", name])

    # Format
    if format_upper != "AUTO":
        cmd.extend(["-f", format_upper])

    # Output directory
    if outdir_str:
        cmd.extend(["--outdir", outdir_str])

    # bdg
    if bdg:
        cmd.append("-B")

    # trackline
    if trackline:
        cmd.append("--trackline")

    # gsize
    if gsize:
        cmd.extend(["-g", gsize])

    # tsize
    if tsize > 0:
        cmd.extend(["-s", str(tsize)])

    # qvalue or pvalue
    if pvalue > 0:
        cmd.extend(["-p", str(pvalue)])
    else:
        cmd.extend(["-q", str(qvalue)])

    # min_length
    if min_length > 0:
        cmd.extend(["--min-length", str(min_length)])

    # max_gap
    if max_gap > 0:
        cmd.extend(["--max-gap", str(max_gap)])

    # nolambda
    if nolambda:
        cmd.append("--nolambda")

    # slocal and llocal
    cmd.extend(["--slocal", str(slocal)])
    cmd.extend(["--llocal", str(llocal)])

    # nomodel
    if nomodel:
        cmd.append("--nomodel")

    # extsize
    if extsize > 0:
        cmd.extend(["--extsize", str(extsize)])

    # shift
    if shift != 0:
        cmd.extend(["--shift", str(shift)])

    # keep_dup
    if isinstance(keep_dup, int):
        cmd.extend(["--keep-dup", str(keep_dup)])
    else:
        cmd.extend(["--keep-dup", keep_dup])

    # broad
    if broad:
        cmd.append("--broad")
        cmd.extend(["--broad-cutoff", str(broad_cutoff)])

    # scale_to
    if scale_to != "small":
        cmd.extend(["--scale-to", scale_to])

    # call_summits
    if call_summits:
        cmd.append("--call-summits")

    # buffer_size
    if buffer_size != 100000:
        cmd.extend(["--buffer-size", str(buffer_size)])

    # cutoff_analysis
    if cutoff_analysis:
        cmd.append("--cutoff-analysis")

    # barcodes
    if barcodes is not None:
        cmd.extend(["--barcodes", str(barcodes.resolve())])

    # max_count
    if max_count is not None:
        cmd.extend(["--max-count", str(max_count)])

    # Run command
    try:
        completed = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"MACS3 callpeak failed with return code {e.returncode}",
        }

    # Collect output files expected based on name and outdir
    output_files = []
    base_path = Path(outdir_str) if outdir_str else Path.cwd()
    # Required output files always generated:
    # NAME_peaks.xls, NAME_peaks.narrowPeak, NAME_summits.bed, NAME_model.r
    output_files.append(str(base_path / f"{name}_peaks.xls"))
    output_files.append(str(base_path / f"{name}_peaks.narrowPeak"))
    output_files.append(str(base_path / f"{name}_summits.bed"))
    output_files.append(str(base_path / f"{name}_model.r"))
    # Optional files
    if broad:
        output_files.append(str(base_path / f"{name}_peaks.broadPeak"))
        output_files.append(str(base_path / f"{name}_peaks.gappedPeak"))
    if bdg:
        output_files.append(str(base_path / f"{name}_treat_pileup.bdg"))
        output_files.append(str(base_path / f"{name}_control_lambda.bdg"))
    if cutoff_analysis:
        output_files.append(str(base_path / f"{name}_cutoff_analysis.txt"))

    return {
        "command_executed": " ".join(cmd),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "output_files": output_files,
    }


if __name__ == '__main__':
    mcp.run()