from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


@mcp.tool()
def bamqc(
    bam: Path,
    paint_chromosome_limits: bool = False,
    cov_hist_lim: int = 50,
    dup_rate_lim: int = 2,
    genome_gc_distr: Optional[str] = None,
    feature_file: Optional[Path] = None,
    homopolymer_min_size: int = 3,
    collect_overlap_pairs: bool = False,
    nr: int = 1000,
    nt: int = 8,
    nw: int = 400,
    output_genome_coverage: Optional[Path] = None,
    outside_stats: bool = False,
    outdir: Optional[Path] = None,
    outfile: str = "report.pdf",
    outformat: str = "HTML",
    sequencing_protocol: str = "non-strand-specific",
    skip_duplicated: bool = False,
    skip_dup_mode: int = 0,
):
    """
    Perform BAM QC analysis on a BAM file.

    Parameters:
    - bam: Input BAM file path.
    - paint_chromosome_limits: Paint chromosome limits inside charts.
    - cov_hist_lim: Upstream limit for targeted per-bin coverage histogram (default 50).
    - dup_rate_lim: Upstream limit for duplication rate histogram (default 2).
    - genome_gc_distr: Species to compare with genome GC distribution: HUMAN or MOUSE.
    - feature_file: Feature file with regions of interest in GFF/GTF or BED format.
    - homopolymer_min_size: Minimum size for homopolymer in indel analysis (default 3).
    - collect_overlap_pairs: Collect statistics of overlapping paired-end reads.
    - nr: Number of reads analyzed in a chunk (default 1000).
    - nt: Number of threads (default 8).
    - nw: Number of windows (default 400).
    - output_genome_coverage: File to save per base non-zero coverage.
    - outside_stats: Report info for regions outside feature-file regions.
    - outdir: Output folder for HTML report and raw data.
    - outfile: Output file for PDF report (default "report.pdf").
    - outformat: Output report format PDF or HTML (default HTML).
    - sequencing_protocol: Library protocol: strand-specific-forward, strand-specific-reverse, or non-strand-specific (default).
    - skip_duplicated: Skip duplicate alignments from analysis.
    - skip_dup_mode: Type of duplicates to skip (0=flagged only, 1=estimated only, 2=both; default 0).
    """
    # Validate input file
    if not bam.exists() or not bam.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam}")

    # Validate feature_file if provided
    if feature_file is not None:
        if not feature_file.exists() or not feature_file.is_file():
            raise FileNotFoundError(f"Feature file not found: {feature_file}")

    # Validate outformat
    outformat_upper = outformat.upper()
    if outformat_upper not in ("PDF", "HTML"):
        raise ValueError("outformat must be 'PDF' or 'HTML'")

    # Validate sequencing_protocol
    valid_protocols = {"strand-specific-forward", "strand-specific-reverse", "non-strand-specific"}
    if sequencing_protocol not in valid_protocols:
        raise ValueError(f"sequencing_protocol must be one of {valid_protocols}")

    # Validate skip_dup_mode
    if skip_dup_mode not in (0, 1, 2):
        raise ValueError("skip_dup_mode must be 0, 1, or 2")

    # Prepare output directory
    if outdir is None:
        outdir = bam.parent / (bam.stem + "_qualimap")
    outdir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [
        "qualimap", "bamqc",
        "-bam", str(bam),
        "-cl", str(cov_hist_lim),
        "-dl", str(dup_rate_lim),
        "-hm", str(homopolymer_min_size),
        "-nr", str(nr),
        "-nt", str(nt),
        "-nw", str(nw),
        "-outdir", str(outdir),
        "-outfile", outfile,
        "-outformat", outformat_upper,
        "-p", sequencing_protocol,
        "-sdmode", str(skip_dup_mode),
    ]

    if paint_chromosome_limits:
        cmd.append("-c")
    if genome_gc_distr is not None:
        genome_gc_distr_upper = genome_gc_distr.upper()
        if genome_gc_distr_upper not in ("HUMAN", "MOUSE"):
            raise ValueError("genome_gc_distr must be 'HUMAN' or 'MOUSE'")
        cmd.extend(["-gd", genome_gc_distr_upper])
    if feature_file is not None:
        cmd.extend(["-gff", str(feature_file)])
    if collect_overlap_pairs:
        cmd.append("-ip")
    if output_genome_coverage is not None:
        cmd.extend(["-oc", str(output_genome_coverage)])
    if outside_stats:
        cmd.append("-os")
    if skip_duplicated:
        cmd.append("-sd")

    # Run command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap bamqc failed with exit code {e.returncode}"
        }

    # Collect output files: HTML report folder and PDF if generated
    output_files = []
    if outdir.exists():
        output_files.append(str(outdir.resolve()))
    pdf_path = outdir / outfile
    if pdf_path.exists():
        output_files.append(str(pdf_path.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def rnaseq(
    bam: Path,
    gtf: Path,
    algorithm: str = "uniquely-mapped-reads",
    num_pr_bases: int = 100,
    num_tr_bias: int = 1000,
    output_counts: Optional[Path] = None,
    outdir: Optional[Path] = None,
    outfile: str = "report.pdf",
    outformat: str = "HTML",
    sequencing_protocol: str = "non-strand-specific",
    paired: bool = False,
    sorted_flag: bool = False,
):
    """
    Perform RNA-seq QC analysis.

    Parameters:
    - bam: Input BAM file path.
    - gtf: Annotations file in Ensembl GTF format.
    - algorithm: Counting algorithm: uniquely-mapped-reads (default) or proportional.
    - num_pr_bases: Number of upstream/downstream bases to compute 5'-3' bias (default 100).
    - num_tr_bias: Number of top highly expressed transcripts to compute 5'-3' bias (default 1000).
    - output_counts: Path to output computed counts.
    - outdir: Output folder for HTML report and raw data.
    - outfile: Output file for PDF report (default "report.pdf").
    - outformat: Output report format PDF or HTML (default HTML).
    - sequencing_protocol: Library protocol: strand-specific-forward, strand-specific-reverse, or non-strand-specific (default).
    - paired: Flag for paired-end experiments (count fragments instead of reads).
    - sorted_flag: Flag indicating input BAM is sorted by name.
    """
    # Validate input files
    if not bam.exists() or not bam.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam}")
    if not gtf.exists() or not gtf.is_file():
        raise FileNotFoundError(f"GTF file not found: {gtf}")

    # Validate algorithm
    if algorithm not in ("uniquely-mapped-reads", "proportional"):
        raise ValueError("algorithm must be 'uniquely-mapped-reads' or 'proportional'")

    # Validate outformat
    outformat_upper = outformat.upper()
    if outformat_upper not in ("PDF", "HTML"):
        raise ValueError("outformat must be 'PDF' or 'HTML'")

    # Validate sequencing_protocol
    valid_protocols = {"strand-specific-forward", "strand-specific-reverse", "non-strand-specific"}
    if sequencing_protocol not in valid_protocols:
        raise ValueError(f"sequencing_protocol must be one of {valid_protocols}")

    # Prepare output directory
    if outdir is None:
        outdir = bam.parent / (bam.stem + "_rnaseq_qualimap")
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "qualimap", "rnaseq",
        "-bam", str(bam),
        "-gtf", str(gtf),
        "-a", algorithm,
        "-npb", str(num_pr_bases),
        "-ntb", str(num_tr_bias),
        "-outdir", str(outdir),
        "-outfile", outfile,
        "-outformat", outformat_upper,
        "-p", sequencing_protocol,
    ]

    if output_counts is not None:
        cmd.extend(["-oc", str(output_counts)])
    if paired:
        cmd.append("-pe")
    if sorted_flag:
        cmd.append("-s")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap rnaseq failed with exit code {e.returncode}"
        }

    output_files = []
    if outdir.exists():
        output_files.append(str(outdir.resolve()))
    pdf_path = outdir / outfile
    if pdf_path.exists():
        output_files.append(str(pdf_path.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def multi_bamqc(
    data: Path,
    paint_chromosome_limits: bool = False,
    feature_file: Optional[Path] = None,
    homopolymer_min_size: int = 3,
    nr: int = 1000,
    nw: int = 400,
    outdir: Optional[Path] = None,
    outfile: str = "report.pdf",
    outformat: str = "HTML",
    run_bamqc: bool = False,
):
    """
    Perform multi-sample BAM QC analysis.

    Parameters:
    - data: File describing input data (2- or 3-column tab-delimited).
    - paint_chromosome_limits: Paint chromosome limits inside charts (only for -r mode).
    - feature_file: Feature file with regions of interest in GFF/GTF or BED format (only for -r mode).
    - homopolymer_min_size: Minimum size for homopolymer in indel analysis (default 3, only for -r mode).
    - nr: Number of reads analyzed in a chunk (default 1000, only for -r mode).
    - nw: Number of windows (default 400, only for -r mode).
    - outdir: Output folder for HTML report and raw data.
    - outfile: Output file for PDF report (default "report.pdf").
    - outformat: Output report format PDF or HTML (default HTML).
    - run_bamqc: If True, run BAM QC first for each sample (-r mode).
    """
    if not data.exists() or not data.is_file():
        raise FileNotFoundError(f"Data file not found: {data}")

    outformat_upper = outformat.upper()
    if outformat_upper not in ("PDF", "HTML"):
        raise ValueError("outformat must be 'PDF' or 'HTML'")

    if outdir is None:
        outdir = data.parent / (data.stem + "_multi_bamqc_qualimap")
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "qualimap", "multi-bamqc",
        "-d", str(data),
        "-outdir", str(outdir),
        "-outfile", outfile,
        "-outformat", outformat_upper,
    ]

    if paint_chromosome_limits:
        cmd.append("-c")
    if feature_file is not None:
        cmd.extend(["-gff", str(feature_file)])
    if homopolymer_min_size != 3:
        cmd.extend(["-hm", str(homopolymer_min_size)])
    if nr != 1000:
        cmd.extend(["-nr", str(nr)])
    if nw != 400:
        cmd.extend(["-nw", str(nw)])
    if run_bamqc:
        cmd.append("-r")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap multi-bamqc failed with exit code {e.returncode}"
        }

    output_files = []
    if outdir.exists():
        output_files.append(str(outdir.resolve()))
    pdf_path = outdir / outfile
    if pdf_path.exists():
        output_files.append(str(pdf_path.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def counts(
    data: Path,
    compare: bool = False,
    info: Optional[Path] = None,
    threshold: Optional[int] = None,
    outdir: Optional[Path] = None,
    outfile: str = "report.pdf",
    outformat: str = "HTML",
    rscriptpath: Optional[Path] = None,
    species: Optional[str] = None,
):
    """
    Perform counts QC analysis.

    Parameters:
    - data: File describing input data (4-column tab-delimited).
    - compare: Perform comparison of conditions (max 2).
    - info: Path to info file with gene GC-content, length, and type.
    - threshold: Threshold for number of counts.
    - outdir: Output folder for HTML report and raw data.
    - outfile: Output file for PDF report (default "report.pdf").
    - outformat: Output report format PDF or HTML (default HTML).
    - rscriptpath: Path to Rscript executable (default assumes in system PATH).
    - species: Use built-in info file for species: HUMAN or MOUSE.
    """
    if not data.exists() or not data.is_file():
        raise FileNotFoundError(f"Data file not found: {data}")

    outformat_upper = outformat.upper()
    if outformat_upper not in ("PDF", "HTML"):
        raise ValueError("outformat must be 'PDF' or 'HTML'")

    if species is not None:
        species_upper = species.upper()
        if species_upper not in ("HUMAN", "MOUSE"):
            raise ValueError("species must be 'HUMAN' or 'MOUSE'")
    else:
        species_upper = None

    if outdir is None:
        outdir = data.parent / (data.stem + "_counts_qualimap")
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "qualimap", "counts",
        "-d", str(data),
        "-outdir", str(outdir),
        "-outfile", outfile,
        "-outformat", outformat_upper,
    ]

    if compare:
        cmd.append("-c")
    if info is not None:
        if not info.exists() or not info.is_file():
            raise FileNotFoundError(f"Info file not found: {info}")
        cmd.extend(["-i", str(info)])
    if threshold is not None:
        if threshold < 0:
            raise ValueError("threshold must be non-negative")
        cmd.extend(["-k", str(threshold)])
    if rscriptpath is not None:
        if not rscriptpath.exists() or not rscriptpath.is_file():
            raise FileNotFoundError(f"Rscript executable not found: {rscriptpath}")
        cmd.extend(["-R", str(rscriptpath)])
    if species_upper is not None:
        cmd.extend(["-s", species_upper])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap counts failed with exit code {e.returncode}"
        }

    output_files = []
    if outdir.exists():
        output_files.append(str(outdir.resolve()))
    pdf_path = outdir / outfile
    if pdf_path.exists():
        output_files.append(str(pdf_path.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def clustering(
    sample: List[Path],
    control: List[Path],
    regions: Path,
    bin_size: int = 100,
    clusters: str = "",
    expr: Optional[str] = None,
    fragment_length: Optional[int] = None,
    upstream_offset: int = 2000,
    downstream_offset: int = 500,
    names: Optional[List[str]] = None,
    outdir: Optional[Path] = None,
    outformat: str = "HTML",
    viz: Optional[str] = None,
):
    """
    Perform clustering of epigenomic signals.

    Parameters:
    - sample: List of sample BAM file paths (comma-separated).
    - control: List of control BAM file paths (comma-separated).
    - regions: Path to regions file.
    - bin_size: Size of the bin (default 100).
    - clusters: Comma-separated list of cluster sizes.
    - expr: Name of the experiment.
    - fragment_length: Smoothing length of a fragment.
    - upstream_offset: Upstream offset (default 2000).
    - downstream_offset: Downstream offset (default 500).
    - names: Comma-separated names of replicates.
    - outdir: Output folder.
    - outformat: Output report format PDF or HTML (default HTML).
    - viz: Visualization type: heatmap or line.
    """
    # Validate input files
    for f in sample:
        if not f.exists() or not f.is_file():
            raise FileNotFoundError(f"Sample BAM file not found: {f}")
    for f in control:
        if not f.exists() or not f.is_file():
            raise FileNotFoundError(f"Control BAM file not found: {f}")
    if not regions.exists() or not regions.is_file():
        raise FileNotFoundError(f"Regions file not found: {regions}")

    outformat_upper = outformat.upper()
    if outformat_upper not in ("PDF", "HTML"):
        raise ValueError("outformat must be 'PDF' or 'HTML'")

    if viz is not None and viz not in ("heatmap", "line"):
        raise ValueError("viz must be 'heatmap' or 'line'")

    if outdir is None:
        outdir = regions.parent / "clustering_qualimap"
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "qualimap", "clustering",
        "-sample", ",".join(str(p) for p in sample),
        "-control", ",".join(str(p) for p in control),
        "-regions", str(regions),
        "-b", str(bin_size),
        "-l", str(upstream_offset),
        "-r", str(downstream_offset),
        "-outdir", str(outdir),
        "-outformat", outformat_upper,
    ]

    if clusters:
        cmd.extend(["-c", clusters])
    if expr is not None:
        cmd.extend(["-expr", expr])
    if fragment_length is not None:
        cmd.extend(["-f", str(fragment_length)])
    if names is not None and len(names) > 0:
        cmd.extend(["-name", ",".join(names)])
    if viz is not None:
        cmd.extend(["-viz", viz])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap clustering failed with exit code {e.returncode}"
        }

    output_files = []
    if outdir.exists():
        output_files.append(str(outdir.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


@mcp.tool()
def comp_counts(
    bam: Path,
    gtf: Path,
    algorithm: str = "uniquely-mapped-reads",
    attribute_id: str = "gene_id",
    out: Optional[Path] = None,
    sequencing_protocol: str = "non-strand-specific",
    paired: bool = False,
    sorted_flag: Optional[str] = None,
    feature_type: str = "exon",
):
    """
    Compute counts from mapping data.

    Parameters:
    - bam: Mapping file in BAM format.
    - gtf: Region file in GTF, GFF or BED format.
    - algorithm: Counting algorithm: uniquely-mapped-reads (default) or proportional.
    - attribute_id: GTF attribute to be used as feature ID (default "gene_id").
    - out: Path to output file.
    - sequencing_protocol: Library protocol: strand-specific-forward, strand-specific-reverse, or non-strand-specific (default).
    - paired: Flag for paired-end experiments (count fragments instead of reads).
    - sorted_flag: Indicates if input file is sorted by name (only for paired-end).
    - feature_type: Value of third column of GTF considered for counting (default "exon").
    """
    if not bam.exists() or not bam.is_file():
        raise FileNotFoundError(f"BAM file not found: {bam}")
    if not gtf.exists() or not gtf.is_file():
        raise FileNotFoundError(f"GTF file not found: {gtf}")

    valid_algorithms = {"uniquely-mapped-reads", "proportional"}
    if algorithm not in valid_algorithms:
        raise ValueError(f"algorithm must be one of {valid_algorithms}")

    valid_protocols = {"strand-specific-forward", "strand-specific-reverse", "non-strand-specific"}
    if sequencing_protocol not in valid_protocols:
        raise ValueError(f"sequencing_protocol must be one of {valid_protocols}")

    if out is None:
        out = bam.parent / (bam.stem + ".counts")

    cmd = [
        "qualimap", "comp-counts",
        "-bam", str(bam),
        "-gtf", str(gtf),
        "-a", algorithm,
        "-id", attribute_id,
        "-out", str(out),
        "-p", sequencing_protocol,
        "-type", feature_type,
    ]

    if paired:
        cmd.append("-pe")
    if sorted_flag is not None:
        cmd.extend(["-s", sorted_flag])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Qualimap comp-counts failed with exit code {e.returncode}"
        }

    output_files = []
    if out.exists():
        output_files.append(str(out.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()