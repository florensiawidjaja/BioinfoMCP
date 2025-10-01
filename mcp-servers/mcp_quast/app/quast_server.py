from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


@mcp.tool()
def quast(
    contig_files: List[Path],
    output_dir: Optional[Path] = None,
    reference: Optional[Path] = None,
    features: Optional[str] = None,
    min_contig: int = 500,
    threads: int = 1,
    split_scaffolds: bool = False,
    labels: Optional[str] = None,
    labels_from_parent_dir: bool = False,
    eukaryote: bool = False,
    fungus: bool = False,
    large: bool = False,
    kmer_stats: bool = False,
    kmer_size: int = 101,
    circos: bool = False,
    gene_finding: bool = False,
    mgm: bool = False,
    glimmer: bool = False,
    gene_thresholds: Optional[str] = None,
    rna_finding: bool = False,
    conserved_genes_finding: bool = False,
    operons: Optional[Path] = None,
    est_ref_size: Optional[int] = None,
    contig_thresholds: Optional[str] = None,
    x_for_nx: int = 90,
    use_all_alignments: bool = False,
    min_alignment: int = 65,
    min_identity: float = 95.0,
    ambiguity_usage: str = "one",
    ambiguity_score: float = 0.99,
    strict_na: bool = False,
    extensive_mis_size: Optional[int] = None,
    local_mis_size: int = 200,
    scaffold_gap_max_size: int = 10000,
    unaligned_part_size: int = 500,
    skip_unaligned_mis_contigs: bool = False,
    fragmented: bool = False,
    fragmented_max_indent: int = 50,
    upper_bound_assembly: bool = False,
    upper_bound_min_con: Optional[int] = None,
    est_insert_size: int = 255,
    report_all_metrics: bool = False,
    plots_format: str = "pdf",
    memory_efficient: bool = False,
    space_efficient: bool = False,
    no_check: bool = False,
    no_plots: bool = False,
    no_krona: bool = False,
    no_html: bool = False,
    no_icarus: bool = False,
    no_snps: bool = False,
    no_gc: bool = False,
    no_sv: bool = False,
    no_read_stats: bool = False,
    fast: bool = False,
    silent: bool = False,
    test: bool = False,
    test_sv: bool = False,
    pe1: Optional[List[Path]] = None,
    pe2: Optional[List[Path]] = None,
    pe12: Optional[List[Path]] = None,
    mp1: Optional[List[Path]] = None,
    mp2: Optional[List[Path]] = None,
    mp12: Optional[List[Path]] = None,
    single: Optional[List[Path]] = None,
    pacbio: Optional[List[Path]] = None,
    nanopore: Optional[List[Path]] = None,
    ref_bam: Optional[Path] = None,
    ref_sam: Optional[Path] = None,
    bam: Optional[str] = None,
    sam: Optional[str] = None,
    sv_bedpe: Optional[Path] = None,
):
    """
    QUAST: Quality Assessment Tool for genome assemblies.
    Evaluates genome assemblies by computing various metrics and generating reports.
    Supports multiple assemblies, references, gene annotations, reads, and advanced options.

    Parameters:
    - contig_files: List of input assembly contig FASTA files (compressed allowed).
    - output_dir: Output directory path. Default is quast_results/results_<date_time>.
    - reference: Reference genome FASTA file (optional).
    - features: Genomic features file path or prefix with feature type (e.g. "CDS:features.gff").
    - min_contig: Minimum contig length to consider (default 500).
    - threads: Number of threads to use (default 1).
    - split_scaffolds: Treat assemblies as scaffolds and split by Ns.
    - labels: Comma-separated assembly labels.
    - labels_from_parent_dir: Use parent directory names as assembly labels.
    - eukaryote: Genome is eukaryotic (affects gene finding and other modules).
    - fungus: Genome is fungal (special mode for gene finding).
    - large: Optimize for large genomes (>100 Mbp).
    - kmer_stats: Compute k-mer based quality metrics.
    - kmer_size: k-mer size for kmer_stats (default 101).
    - circos: Generate Circos plot files.
    - gene_finding: Enable gene finding.
    - mgm: Force use of MetaGeneMark for gene finding.
    - glimmer: Use GlimmerHMM for gene finding.
    - gene_thresholds: Comma-separated gene length thresholds.
    - rna_finding: Enable ribosomal RNA gene finding.
    - conserved_genes_finding: Enable BUSCO conserved genes search.
    - operons: File with operon positions.
    - est_ref_size: Estimated reference genome size if no reference provided.
    - contig_thresholds: Comma-separated contig length thresholds.
    - x_for_nx: x value for Nx metrics (0-100, default 90).
    - use_all_alignments: Use all alignments instead of filtering ambiguous ones.
    - min_alignment: Minimum alignment length (default 65).
    - min_identity: Minimum alignment identity % (default 95.0).
    - ambiguity_usage: How to handle ambiguous alignments ('none','one','all').
    - ambiguity_score: Score threshold for ambiguous alignments (0.8-1.0).
    - strict_na: Break contigs at all misassemblies including local ones.
    - extensive_mis_size: Threshold for extensive misassemblies relocation size.
    - local_mis_size: Threshold for local misassemblies size (default 200).
    - scaffold_gap_max_size: Max scaffold gap size for misassembly detection (default 10000).
    - unaligned_part_size: Threshold for unaligned fragment size (default 500).
    - skip_unaligned_mis_contigs: Do not distinguish contigs with >50% unaligned bases.
    - fragmented: Reference genome is fragmented.
    - fragmented_max_indent: Max indent for fragmented translocation detection (default 50).
    - upper_bound_assembly: Simulate upper bound assembly using reads.
    - upper_bound_min_con: Minimal connecting reads for upper bound scaffolding.
    - est_insert_size: Estimated paired reads insert size (default 255).
    - report_all_metrics: Keep all quality metrics in report.
    - plots_format: Format for plots (emf, eps, pdf, png, ps, raw, rgba, svg, svgz).
    - memory_efficient: Use memory efficient mode (slower).
    - space_efficient: Use space efficient mode (less aux files).
    - no_check: Skip input FASTA checks and corrections.
    - no_plots: Do not generate static plots.
    - no_krona: Do not generate Krona pie charts.
    - no_html: Do not build HTML reports.
    - no_icarus: Do not build Icarus viewers.
    - no_snps: Do not report SNP statistics.
    - no_gc: Do not compute GC content and plots.
    - no_sv: Do not run structural variant calling.
    - no_read_stats: Do not align reads against assemblies.
    - fast: Shortcut for speedup options except no-check.
    - silent: Suppress detailed stdout output.
    - test: Run test on test_data folder.
    - test_sv: Run test with SV detection on test_data.
    - pe1, pe2, pe12: Lists of paired-end read files (forward, reverse, interlaced).
    - mp1, mp2, mp12: Lists of mate-pair read files (forward, reverse, interlaced).
    - single: List of unpaired read files.
    - pacbio: List of PacBio reads files.
    - nanopore: List of Oxford Nanopore reads files.
    - ref_bam: BAM file with reads aligned to reference.
    - ref_sam: SAM file with reads aligned to reference.
    - bam: Comma-separated BAM files aligned to assemblies.
    - sam: Comma-separated SAM files aligned to assemblies.
    - sv_bedpe: BEDPE file with structural variations.

    Returns:
    Dict with keys: command_executed, stdout, stderr, output_files (list).
    """
    import shlex

    # Validate contig files
    if not contig_files or len(contig_files) == 0:
        raise ValueError("At least one contig file must be specified.")
    for f in contig_files:
        if not f.exists():
            raise FileNotFoundError(f"Contig file not found: {f}")

    # Validate optional files
    if reference is not None and not reference.exists():
        raise FileNotFoundError(f"Reference genome file not found: {reference}")
    if features is not None:
        # features can be prefixed with feature type like "CDS:filename"
        if ":" in features:
            _, feat_path = features.split(":", 1)
            feat_path_obj = Path(feat_path)
            if not feat_path_obj.exists():
                raise FileNotFoundError(f"Features file not found: {feat_path_obj}")
        else:
            feat_path_obj = Path(features)
            if not feat_path_obj.exists():
                raise FileNotFoundError(f"Features file not found: {feat_path_obj}")
    if operons is not None and not operons.exists():
        raise FileNotFoundError(f"Operons file not found: {operons}")
    if ref_bam is not None and not ref_bam.exists():
        raise FileNotFoundError(f"Reference BAM file not found: {ref_bam}")
    if ref_sam is not None and not ref_sam.exists():
        raise FileNotFoundError(f"Reference SAM file not found: {ref_sam}")
    if sv_bedpe is not None and not sv_bedpe.exists():
        raise FileNotFoundError(f"SV BEDPE file not found: {sv_bedpe}")

    # Validate numeric parameters
    if min_contig < 0:
        raise ValueError("min_contig must be non-negative")
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if kmer_size < 1:
        raise ValueError("kmer_size must be positive")
    if not (0 <= x_for_nx <= 100):
        raise ValueError("x_for_nx must be between 0 and 100")
    if not (0.8 <= ambiguity_score <= 1.0):
        raise ValueError("ambiguity_score must be between 0.8 and 1.0")
    if local_mis_size < 0:
        raise ValueError("local_mis_size must be non-negative")
    if scaffold_gap_max_size < 0:
        raise ValueError("scaffold_gap_max_size must be non-negative")
    if unaligned_part_size < 0:
        raise ValueError("unaligned_part_size must be non-negative")
    if fragmented_max_indent < 0:
        raise ValueError("fragmented_max_indent must be non-negative")
    if est_insert_size < 0:
        raise ValueError("est_insert_size must be non-negative")
    if extensive_mis_size is not None and extensive_mis_size < 0:
        raise ValueError("extensive_mis_size must be non-negative if specified")
    if upper_bound_min_con is not None and upper_bound_min_con < 0:
        raise ValueError("upper_bound_min_con must be non-negative if specified")

    # Validate ambiguity_usage
    if ambiguity_usage not in ("none", "one", "all"):
        raise ValueError("ambiguity_usage must be one of: 'none', 'one', 'all'")

    # Validate plots_format
    valid_plot_formats = {"emf", "eps", "pdf", "png", "ps", "raw", "rgba", "svg", "svgz"}
    if plots_format.lower() not in valid_plot_formats:
        raise ValueError(f"plots_format must be one of {valid_plot_formats}")

    # Prepare command line
    cmd = ["quast.py"]

    # Add contig files
    for f in contig_files:
        cmd.append(str(f))

    # Output directory
    if output_dir is not None:
        output_dir = Path(output_dir)
        cmd.extend(["-o", str(output_dir)])

    # Reference genome
    if reference is not None:
        cmd.extend(["-r", str(reference)])

    # Features
    if features is not None:
        cmd.extend(["--features", features])

    # Min contig length
    cmd.extend(["-m", str(min_contig)])

    # Threads
    cmd.extend(["-t", str(threads)])

    # Split scaffolds
    if split_scaffolds:
        cmd.append("-s")

    # Labels
    if labels is not None:
        cmd.extend(["-l", labels])

    # Labels from parent directory
    if labels_from_parent_dir:
        cmd.append("-L")

    # Eukaryote
    if eukaryote:
        cmd.append("-e")

    # Fungus
    if fungus:
        cmd.append("--fungus")

    # Large genome mode
    if large:
        cmd.append("--large")

    # k-mer stats
    if kmer_stats:
        cmd.append("-k")

    # k-mer size
    if kmer_stats and kmer_size != 101:
        cmd.extend(["--k-mer-size", str(kmer_size)])

    # Circos plot
    if circos:
        cmd.append("--circos")

    # Gene finding
    if gene_finding:
        cmd.append("-f")

    # Force MetaGeneMark
    if mgm:
        cmd.append("--mgm")

    # GlimmerHMM gene finding
    if glimmer:
        cmd.append("--glimmer")

    # Gene thresholds
    if gene_thresholds is not None:
        cmd.extend(["--gene-thresholds", gene_thresholds])

    # RNA finding
    if rna_finding:
        cmd.append("--rna-finding")

    # Conserved genes finding (BUSCO)
    if conserved_genes_finding:
        cmd.append("-b")

    # Operons file
    if operons is not None:
        cmd.extend(["--operons", str(operons)])

    # Estimated reference size
    if est_ref_size is not None:
        cmd.extend(["--est-ref-size", str(est_ref_size)])

    # Contig thresholds
    if contig_thresholds is not None:
        cmd.extend(["--contig-thresholds", contig_thresholds])

    # x for Nx metrics
    if x_for_nx != 90:
        cmd.extend(["--x-for-Nx", str(x_for_nx)])

    # Use all alignments
    if use_all_alignments:
        cmd.append("-u")

    # Minimum alignment length
    if min_alignment != 65:
        cmd.extend(["-i", str(min_alignment)])

    # Minimum identity
    if min_identity != 95.0:
        cmd.extend(["--min-identity", str(min_identity)])

    # Ambiguity usage
    if ambiguity_usage != "one":
        cmd.extend(["-a", ambiguity_usage])

    # Ambiguity score
    if ambiguity_score != 0.99:
        cmd.extend(["--ambiguity-score", str(ambiguity_score)])

    # Strict NA
    if strict_na:
        cmd.append("--strict-NA")

    # Extensive mis-size
    if extensive_mis_size is not None:
        cmd.extend(["-x", str(extensive_mis_size)])

    # Local mis-size
    if local_mis_size != 200:
        cmd.extend(["--local-mis-size", str(local_mis_size)])

    # Scaffold gap max size
    if scaffold_gap_max_size != 10000:
        cmd.extend(["--scaffold-gap-max-size", str(scaffold_gap_max_size)])

    # Unaligned part size
    if unaligned_part_size != 500:
        cmd.extend(["--unaligned-part-size", str(unaligned_part_size)])

    # Skip unaligned mis contigs
    if skip_unaligned_mis_contigs:
        cmd.append("--skip-unaligned-mis-contigs")

    # Fragmented reference
    if fragmented:
        cmd.append("--fragmented")

    # Fragmented max indent
    if fragmented_max_indent != 50:
        cmd.extend(["--fragmented-max-indent", str(fragmented_max_indent)])

    # Upper bound assembly
    if upper_bound_assembly:
        cmd.append("--upper-bound-assembly")

    # Upper bound min connecting reads
    if upper_bound_min_con is not None:
        cmd.extend(["--upper-bound-min-con", str(upper_bound_min_con)])

    # Estimated insert size
    if est_insert_size != 255:
        cmd.extend(["--est-insert-size", str(est_insert_size)])

    # Report all metrics
    if report_all_metrics:
        cmd.append("--report-all-metrics")

    # Plots format
    if plots_format.lower() != "pdf":
        cmd.extend(["--plots-format", plots_format.lower()])

    # Memory efficient
    if memory_efficient:
        cmd.append("--memory-efficient")

    # Space efficient
    if space_efficient:
        cmd.append("--space-efficient")

    # No check
    if no_check:
        cmd.append("--no-check")

    # No plots
    if no_plots:
        cmd.append("--no-plots")

    # No krona
    if no_krona:
        cmd.append("--no-krona")

    # No html
    if no_html:
        cmd.append("--no-html")

    # No icarus
    if no_icarus:
        cmd.append("--no-icarus")

    # No snps
    if no_snps:
        cmd.append("--no-snps")

    # No gc
    if no_gc:
        cmd.append("--no-gc")

    # No sv
    if no_sv:
        cmd.append("--no-sv")

    # No read stats
    if no_read_stats:
        cmd.append("--no-read-stats")

    # Fast shortcut
    if fast:
        cmd.append("--fast")

    # Silent
    if silent:
        cmd.append("--silent")

    # Test
    if test:
        cmd.append("--test")

    # Test SV
    if test_sv:
        cmd.append("--test-sv")

    # Reads options: multiple allowed, add multiple times
    def add_multiple_files(option_name: str, files: Optional[List[Path]]):
        if files:
            for f in files:
                if not f.exists():
                    raise FileNotFoundError(f"Reads file not found: {f}")
                cmd.extend([option_name, str(f)])

    add_multiple_files("--pe1", pe1)
    add_multiple_files("--pe2", pe2)
    add_multiple_files("--pe12", pe12)
    add_multiple_files("--mp1", mp1)
    add_multiple_files("--mp2", mp2)
    add_multiple_files("--mp12", mp12)
    add_multiple_files("--single", single)
    add_multiple_files("--pacbio", pacbio)
    add_multiple_files("--nanopore", nanopore)

    # Reference BAM/SAM
    if ref_bam is not None:
        cmd.extend(["--ref-bam", str(ref_bam)])
    if ref_sam is not None:
        cmd.extend(["--ref-sam", str(ref_sam)])

    # BAM/SAM aligned to assemblies (comma-separated)
    if bam is not None:
        # Validate files exist
        bam_files = [Path(x) for x in bam.split(",")]
        for bf in bam_files:
            if not bf.exists():
                raise FileNotFoundError(f"BAM file not found: {bf}")
        cmd.extend(["--bam", bam])

    if sam is not None:
        sam_files = [Path(x) for x in sam.split(",")]
        for sf in sam_files:
            if not sf.exists():
                raise FileNotFoundError(f"SAM file not found: {sf}")
        cmd.extend(["--sam", sam])

    # SV BEDPE file
    if sv_bedpe is not None:
        cmd.extend(["--sv-bedpe", str(sv_bedpe)])

    # Execute command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(shlex.quote(c) for c in cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"QUAST execution failed with return code {e.returncode}",
        }

    # Determine output directory for reporting files
    if output_dir is None:
        # Default output dir pattern: quast_results/results_<date_time>
        # We cannot reliably guess exact directory name, so fallback to quast_results/latest if exists
        default_latest = Path("quast_results/latest")
        if default_latest.exists() and default_latest.is_dir():
            output_dir = default_latest
        else:
            output_dir = None

    # Collect main output files if output_dir is known
    output_files = []
    if output_dir is not None:
        # Common main report files
        main_reports = [
            "report.txt",
            "report.tsv",
            "report.tex",
            "icarus.html",
            "report.pdf",
            "report.html",
        ]
        for f in main_reports:
            p = output_dir / f
            if p.exists():
                output_files.append(str(p))

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


if __name__ == "__main__":
    mcp.run()