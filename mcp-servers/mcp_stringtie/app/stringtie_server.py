from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Union
import subprocess

mcp = FastMCP()


@mcp.tool()
def stringtie_assemble(
    input_bams: List[Path],
    guide_gtf: Optional[Path] = None,
    prefix: str = "STRG",
    output_gtf: Optional[Path] = None,
    cpus: int = 1,
    verbose: bool = False,
    min_anchor_len: int = 10,
    min_len: int = 200,
    min_anchor_cov: int = 1,
    min_iso: float = 0.01,
    min_bundle_cov: float = 1.0,
    max_gap: int = 50,
    no_trim: bool = False,
    min_multi_exon_cov: float = 1.0,
    min_single_exon_cov: float = 4.75,
    long_reads: bool = False,
    clean_only: bool = False,
    viral: bool = False,
    err_margin: int = 25,
    ptf_file: Optional[Path] = None,
    exclude_seqids: Optional[List[str]] = None,
    gene_abund_out: Optional[Path] = None,
    ballgown: bool = False,
    ballgown_dir: Optional[Path] = None,
    estimate_abund_only: bool = False,
    no_multimapping_correction: bool = False,
    mix: bool = False,
    conservative: bool = False,
    stranded_rf: bool = False,
    stranded_fr: bool = False,
    nascent: bool = False,
    nascent_output: bool = False,
    cram_ref: Optional[Path] = None,
) -> dict:
    """
    Run StringTie transcript assembly on one or more BAM/CRAM input files.

    Parameters:
    - input_bams: List of input BAM/CRAM files (at least one).
    - guide_gtf: Reference annotation GTF/GFF file to guide assembly.
    - prefix: Prefix for output transcripts (default: STRG).
    - output_gtf: Output GTF file path (default: stdout).
    - cpus: Number of threads to use (default: 1).
    - verbose: Enable verbose logging.
    - min_anchor_len: Minimum anchor length for junctions (default: 10).
    - min_len: Minimum assembled transcript length (default: 200).
    - min_anchor_cov: Minimum junction coverage (default: 1).
    - min_iso: Minimum isoform fraction (default: 0.01).
    - min_bundle_cov: Minimum reads per bp coverage for multi-exon transcripts (default: 1.0).
    - max_gap: Maximum gap allowed between read mappings (default: 50).
    - no_trim: Disable trimming of predicted transcripts based on coverage.
    - min_multi_exon_cov: Minimum coverage for multi-exon transcripts (default: 1.0).
    - min_single_exon_cov: Minimum coverage for single-exon transcripts (default: 4.75).
    - long_reads: Enable long reads processing.
    - clean_only: If long reads provided, clean and collapse reads but do not assemble.
    - viral: Enable viral mode for long reads.
    - err_margin: Window around erroneous splice sites (default: 25).
    - ptf_file: Load point-features from a 4-column feature file.
    - exclude_seqids: List of reference sequence IDs to exclude from assembly.
    - gene_abund_out: Output file for gene abundance estimation.
    - ballgown: Enable output of Ballgown table files in output GTF directory.
    - ballgown_dir: Directory path to output Ballgown table files.
    - estimate_abund_only: Only estimate abundance of given reference transcripts.
    - no_multimapping_correction: Disable multi-mapping correction.
    - mix: Both short and long read alignments provided (long reads must be 2nd BAM).
    - conservative: Conservative transcript assembly (same as -t -c 1.5 -f 0.05).
    - stranded_rf: Assume stranded library fr-firststrand.
    - stranded_fr: Assume stranded library fr-secondstrand.
    - nascent: Nascent aware assembly for rRNA-depleted RNAseq libraries.
    - nascent_output: Enables nascent and outputs assembled nascent transcripts.
    - cram_ref: Reference genome FASTA file for CRAM input.
    """
    # Validate inputs
    if len(input_bams) == 0:
        raise ValueError("At least one input BAM/CRAM file must be provided.")
    for bam in input_bams:
        if not bam.exists():
            raise FileNotFoundError(f"Input BAM/CRAM file not found: {bam}")
    if guide_gtf is not None and not guide_gtf.exists():
        raise FileNotFoundError(f"Guide GTF/GFF file not found: {guide_gtf}")
    if ptf_file is not None and not ptf_file.exists():
        raise FileNotFoundError(f"Point-feature file not found: {ptf_file}")
    if gene_abund_out is not None:
        gene_abund_out = Path(gene_abund_out)
    if output_gtf is not None:
        output_gtf = Path(output_gtf)
    if ballgown_dir is not None:
        ballgown_dir = Path(ballgown_dir)
        if not ballgown_dir.exists():
            raise FileNotFoundError(f"Ballgown directory does not exist: {ballgown_dir}")
    if cram_ref is not None and not cram_ref.exists():
        raise FileNotFoundError(f"CRAM reference FASTA file not found: {cram_ref}")
    if exclude_seqids is not None:
        if not all(isinstance(s, str) for s in exclude_seqids):
            raise ValueError("exclude_seqids must be a list of strings")

    # Validate numeric parameters
    if cpus < 1:
        raise ValueError("cpus must be >= 1")
    if min_anchor_len < 0:
        raise ValueError("min_anchor_len must be >= 0")
    if min_len < 0:
        raise ValueError("min_len must be >= 0")
    if min_anchor_cov < 0:
        raise ValueError("min_anchor_cov must be >= 0")
    if not (0.0 <= min_iso <= 1.0):
        raise ValueError("min_iso must be between 0 and 1")
    if min_bundle_cov < 0:
        raise ValueError("min_bundle_cov must be >= 0")
    if max_gap < 0:
        raise ValueError("max_gap must be >= 0")
    if min_multi_exon_cov < 0:
        raise ValueError("min_multi_exon_cov must be >= 0")
    if min_single_exon_cov < 0:
        raise ValueError("min_single_exon_cov must be >= 0")
    if err_margin < 0:
        raise ValueError("err_margin must be >= 0")

    # Build command
    cmd = ["stringtie"]
    # Input BAMs
    for bam in input_bams:
        cmd.append(str(bam))
    # Guide annotation
    if guide_gtf:
        cmd.extend(["-G", str(guide_gtf)])
    # Prefix
    if prefix:
        cmd.extend(["-l", prefix])
    # Output GTF
    if output_gtf:
        cmd.extend(["-o", str(output_gtf)])
    else:
        # If no output specified, stringtie writes to stdout by default
        pass
    # CPUs
    cmd.extend(["-p", str(cpus)])
    # Verbose
    if verbose:
        cmd.append("-v")
    # Min anchor length
    cmd.extend(["-a", str(min_anchor_len)])
    # Min transcript length
    cmd.extend(["-m", str(min_len)])
    # Min junction coverage
    cmd.extend(["-j", str(min_anchor_cov)])
    # Min isoform fraction
    cmd.extend(["-f", str(min_iso)])
    # Min bundle coverage (reads per bp coverage for multi-exon)
    cmd.extend(["-c", str(min_bundle_cov)])
    # Max gap
    cmd.extend(["-g", str(max_gap)])
    # No trimming
    if no_trim:
        cmd.append("-t")
    # Coverage thresholds for multi-exon and single-exon transcripts
    cmd.extend(["-c", str(min_multi_exon_cov)])  # -c is min reads per bp coverage multi-exon
    cmd.extend(["-s", str(min_single_exon_cov)])  # -s is min reads per bp coverage single-exon
    # Long reads processing
    if long_reads:
        cmd.append("-L")
    # Clean only (no assembly)
    if clean_only:
        cmd.append("-R")
    # Viral mode
    if viral:
        cmd.append("--viral")
    # Error margin
    cmd.extend(["-E", str(err_margin)])
    # Point features file
    if ptf_file:
        cmd.extend(["--ptf", str(ptf_file)])
    # Exclude seqids
    if exclude_seqids:
        cmd.extend(["-x", ",".join(exclude_seqids)])
    # Gene abundance output
    if gene_abund_out:
        cmd.extend(["-A", str(gene_abund_out)])
    # Ballgown output
    if ballgown:
        cmd.append("-B")
    if ballgown_dir:
        cmd.extend(["-b", str(ballgown_dir)])
    # Estimate abundance only
    if estimate_abund_only:
        cmd.append("-e")
    # No multi-mapping correction
    if no_multimapping_correction:
        cmd.append("-u")
    # Mix mode
    if mix:
        cmd.append("--mix")
    # Conservative mode
    if conservative:
        cmd.append("--conservative")
    # Strandedness
    if stranded_rf:
        cmd.append("--rf")
    if stranded_fr:
        cmd.append("--fr")
    # Nascent
    if nascent:
        cmd.append("-N")
    if nascent_output:
        cmd.append("--nasc")
    # CRAM reference
    if cram_ref:
        cmd.extend(["--cram-ref", str(cram_ref)])

    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"StringTie assembly failed with exit code {e.returncode}",
            "output_files": [],
        }

    output_files = []
    if output_gtf:
        output_files.append(str(output_gtf))
    if gene_abund_out:
        output_files.append(str(gene_abund_out))
    if ballgown_dir:
        # Ballgown files are created inside this directory
        output_files.append(str(ballgown_dir))
    elif ballgown and output_gtf:
        # Ballgown files created in output GTF directory
        output_files.append(str(output_gtf.parent))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files,
    }


@mcp.tool()
def stringtie_merge(
    input_gtfs: List[Path],
    guide_gtf: Optional[Path] = None,
    output_gtf: Optional[Path] = None,
    min_len: int = 50,
    min_cov: float = 0.0,
    min_fpkm: float = 1.0,
    min_tpm: float = 1.0,
    min_iso: float = 0.01,
    max_gap: int = 250,
    keep_retained_introns: bool = False,
    prefix: str = "MSTRG",
) -> dict:
    """
    Merge multiple StringTie GTF files into a unified non-redundant set of isoforms.

    Parameters:
    - input_gtfs: List of input GTF files to merge (at least one).
    - guide_gtf: Reference annotation GTF/GFF3 to include in the merging.
    - output_gtf: Output merged GTF file (default: stdout).
    - min_len: Minimum input transcript length to include (default: 50).
    - min_cov: Minimum input transcript coverage to include (default: 0).
    - min_fpkm: Minimum input transcript FPKM to include (default: 1.0).
    - min_tpm: Minimum input transcript TPM to include (default: 1.0).
    - min_iso: Minimum isoform fraction (default: 0.01).
    - max_gap: Gap between transcripts to merge together (default: 250).
    - keep_retained_introns: Keep merged transcripts with retained introns.
    - prefix: Name prefix for output transcripts (default: MSTRG).
    """
    # Validate inputs
    if len(input_gtfs) == 0:
        raise ValueError("At least one input GTF file must be provided.")
    for gtf in input_gtfs:
        if not gtf.exists():
            raise FileNotFoundError(f"Input GTF file not found: {gtf}")
    if guide_gtf is not None and not guide_gtf.exists():
        raise FileNotFoundError(f"Guide GTF/GFF3 file not found: {guide_gtf}")
    if output_gtf is not None:
        output_gtf = Path(output_gtf)

    # Validate numeric parameters
    if min_len < 0:
        raise ValueError("min_len must be >= 0")
    if min_cov < 0:
        raise ValueError("min_cov must be >= 0")
    if min_fpkm < 0:
        raise ValueError("min_fpkm must be >= 0")
    if min_tpm < 0:
        raise ValueError("min_tpm must be >= 0")
    if not (0.0 <= min_iso <= 1.0):
        raise ValueError("min_iso must be between 0 and 1")
    if max_gap < 0:
        raise ValueError("max_gap must be >= 0")

    # Build command
    cmd = ["stringtie", "--merge"]
    # Guide annotation
    if guide_gtf:
        cmd.extend(["-G", str(guide_gtf)])
    # Output GTF
    if output_gtf:
        cmd.extend(["-o", str(output_gtf)])
    # Min transcript length
    cmd.extend(["-m", str(min_len)])
    # Min coverage
    cmd.extend(["-c", str(min_cov)])
    # Min FPKM
    cmd.extend(["-F", str(min_fpkm)])
    # Min TPM
    cmd.extend(["-T", str(min_tpm)])
    # Min isoform fraction
    cmd.extend(["-f", str(min_iso)])
    # Max gap
    cmd.extend(["-g", str(max_gap)])
    # Keep retained introns
    if keep_retained_introns:
        cmd.append("-i")
    # Prefix
    if prefix:
        cmd.extend(["-l", prefix])
    # Input GTFs
    for gtf in input_gtfs:
        cmd.append(str(gtf))

    # Run command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"StringTie merge failed with exit code {e.returncode}",
            "output_files": [],
        }

    output_files = []
    if output_gtf:
        output_files.append(str(output_gtf))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files,
    }


@mcp.tool()
def stringtie_version() -> dict:
    """
    Print the StringTie version.
    """
    cmd = ["stringtie", "--version"]
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        stdout = result.stdout.strip()
        stderr = result.stderr.strip()
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"StringTie version command failed with exit code {e.returncode}",
            "output_files": [],
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [],
    }


if __name__ == '__main__':
    mcp.run()