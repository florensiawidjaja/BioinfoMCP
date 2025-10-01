from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Union
import subprocess

mcp = FastMCP()


@mcp.tool()
def salmon_index(
    transcripts_fasta: Path,
    index_dir: Path,
    decoys_file: Optional[Path] = None,
    kmer_size: int = 31,
) -> dict:
    """
    Build a Salmon index for the transcriptome.

    Parameters:
    - transcripts_fasta: Path to the FASTA file containing reference transcripts.
    - index_dir: Directory path where the index will be created.
    - decoys_file: Optional path to a file listing decoy sequences.
    - kmer_size: k-mer size for the index (default 31, recommended for reads >=75bp).

    Returns:
    - dict with command executed, stdout, stderr, and output_files (index directory).
    """
    # Validate inputs
    if not transcripts_fasta.is_file():
        raise FileNotFoundError(f"Transcripts FASTA file not found: {transcripts_fasta}")
    if decoys_file is not None and not decoys_file.is_file():
        raise FileNotFoundError(f"Decoys file not found: {decoys_file}")
    if kmer_size <= 0:
        raise ValueError("kmer_size must be a positive integer")
    index_dir = Path(index_dir)
    # Prepare command
    cmd = ["salmon", "index", "-t", str(transcripts_fasta), "-i", str(index_dir), "-k", str(kmer_size)]
    if decoys_file:
        cmd.extend(["--decoys", str(decoys_file)])

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = [str(index_dir)]
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Salmon index failed with exit code {e.returncode}",
        }


@mcp.tool()
def salmon_quant(
    index_or_transcripts: Path,
    lib_type: str,
    output_dir: Path,
    reads_1: Optional[List[Path]] = None,
    reads_2: Optional[List[Path]] = None,
    single_reads: Optional[List[Path]] = None,
    alignments: Optional[List[Path]] = None,
    validate_mappings: bool = False,
    mimic_bt2: bool = False,
    mimic_strict_bt2: bool = False,
    meta: bool = False,
    recover_orphans: bool = False,
    hard_filter: bool = False,
    skip_quant: bool = False,
    allow_dovetail: bool = False,
    threads: int = 0,
    dump_eq: bool = False,
    incompat_prior: float = 0.01,
    fld_mean: Optional[float] = None,
    fld_sd: Optional[float] = None,
    min_score_fraction: Optional[float] = None,
    bandwidth: Optional[int] = None,
    max_mmpextension: Optional[int] = None,
    ma: Optional[int] = None,
    mp: Optional[int] = None,
    go: Optional[int] = None,
    ge: Optional[int] = None,
    range_factorization_bins: Optional[int] = None,
    use_em: bool = False,
    vb_prior: Optional[float] = None,
    per_transcript_prior: bool = False,
    num_bootstraps: int = 0,
    num_gibbs_samples: int = 0,
    seq_bias: bool = False,
    num_bias_samples: Optional[int] = None,
    gc_bias: bool = False,
    pos_bias: bool = False,
    bias_speed_samp: int = 5,
    write_unmapped_names: bool = False,
    write_mappings: Optional[Union[bool, Path]] = False,
) -> dict:
    """
    Quantify transcript abundances using Salmon in mapping-based or alignment-based mode.

    Parameters:
    - index_or_transcripts: Path to Salmon index directory (mapping-based mode) or transcripts FASTA (alignment-based mode).
    - lib_type: Library type string (e.g. IU, SF, OSR, or 'A' for automatic).
    - output_dir: Directory to write quantification results.
    - reads_1: List of paths to left reads files (paired-end).
    - reads_2: List of paths to right reads files (paired-end).
    - single_reads: List of paths to single-end reads files.
    - alignments: List of paths to SAM/BAM alignment files (alignment-based mode).
    - validate_mappings: Enable selective alignment (--validateMappings).
    - mimic_bt2: Mimic Bowtie2 mapping parameters.
    - mimic_strict_bt2: Mimic strict Bowtie2 mapping parameters.
    - meta: Enable metagenomic mode.
    - recover_orphans: Enable orphan rescue (with selective alignment).
    - hard_filter: Use hard filtering (with selective alignment).
    - skip_quant: Skip quantification step.
    - allow_dovetail: Allow dovetailing mappings.
    - threads: Number of threads to use (0 means auto-detect).
    - dump_eq: Dump equivalence classes.
    - incompat_prior: Prior probability for incompatible mappings (default 0.01).
    - fld_mean: Mean fragment length (single-end only).
    - fld_sd: Fragment length standard deviation (single-end only).
    - min_score_fraction: Minimum score fraction for valid mapping (with --validateMappings).
    - bandwidth: Bandwidth for ksw2 alignment (selective alignment).
    - max_mmpextension: Max extension length for selective alignment.
    - ma: Match score for alignment.
    - mp: Mismatch penalty for alignment.
    - go: Gap open penalty.
    - ge: Gap extension penalty.
    - range_factorization_bins: Fidelity parameter for range factorization.
    - use_em: Use EM algorithm instead of VBEM.
    - vb_prior: VBEM prior value.
    - per_transcript_prior: Use per-transcript prior instead of per-nucleotide.
    - num_bootstraps: Number of bootstrap samples.
    - num_gibbs_samples: Number of Gibbs samples (mutually exclusive with bootstraps).
    - seq_bias: Enable sequence-specific bias correction.
    - num_bias_samples: Number of reads to learn sequence bias from.
    - gc_bias: Enable fragment GC bias correction.
    - pos_bias: Enable positional bias correction.
    - bias_speed_samp: Sampling factor for bias speedup (default 5).
    - write_unmapped_names: Write unmapped read names.
    - write_mappings: Write mapping info; False=no, True=stdout, Path=filename.

    Returns:
    - dict with command executed, stdout, stderr, and output_files (output directory).
    """
    # Validate inputs
    index_or_transcripts = Path(index_or_transcripts)
    output_dir = Path(output_dir)
    if not index_or_transcripts.exists():
        raise FileNotFoundError(f"Index directory or transcripts file not found: {index_or_transcripts}")
    if reads_1 is None:
        reads_1 = []
    if reads_2 is None:
        reads_2 = []
    if single_reads is None:
        single_reads = []
    if alignments is None:
        alignments = []
    # Validate read files existence
    for f in reads_1 + reads_2 + single_reads + alignments:
        if not Path(f).exists():
            raise FileNotFoundError(f"Input file not found: {f}")
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if num_bootstraps > 0 and num_gibbs_samples > 0:
        raise ValueError("num_bootstraps and num_gibbs_samples are mutually exclusive")

    cmd = ["salmon", "quant"]

    # Determine mode: mapping-based (index) or alignment-based (transcripts + alignments)
    if index_or_transcripts.is_dir():
        # mapping-based mode
        cmd.extend(["-i", str(index_or_transcripts)])
    else:
        # alignment-based mode
        cmd.extend(["-t", str(index_or_transcripts)])

    cmd.extend(["-l", lib_type])
    cmd.extend(["-o", str(output_dir)])

    # Reads input
    if alignments:
        # alignment-based mode: provide -a with alignment files
        for aln in alignments:
            cmd.extend(["-a", str(aln)])
    else:
        # mapping-based mode: provide reads
        if single_reads:
            # single-end reads
            for r in single_reads:
                cmd.extend(["-r", str(r)])
        else:
            # paired-end reads
            if len(reads_1) == 0 or len(reads_2) == 0:
                raise ValueError("Paired-end reads require both reads_1 and reads_2 lists to be non-empty")
            if len(reads_1) != len(reads_2):
                raise ValueError("reads_1 and reads_2 must have the same number of files")
            for r1 in reads_1:
                cmd.append("-1")
                cmd.append(str(r1))
            for r2 in reads_2:
                cmd.append("-2")
                cmd.append(str(r2))

    # Flags and options
    if validate_mappings:
        cmd.append("--validateMappings")
    if mimic_bt2:
        cmd.append("--mimicBT2")
    if mimic_strict_bt2:
        cmd.append("--mimicStrictBT2")
    if meta:
        cmd.append("--meta")
    if recover_orphans:
        cmd.append("--recoverOrphans")
    if hard_filter:
        cmd.append("--hardFilter")
    if skip_quant:
        cmd.append("--skipQuant")
    if allow_dovetail:
        cmd.append("--allowDovetail")
    if threads > 0:
        cmd.extend(["-p", str(threads)])
    if dump_eq:
        cmd.append("--dumpEq")
    if incompat_prior != 0.01:
        if incompat_prior < 0.0 or incompat_prior > 1.0:
            raise ValueError("incompat_prior must be between 0 and 1")
        cmd.extend(["--incompatPrior", str(incompat_prior)])
    if fld_mean is not None:
        if fld_mean <= 0:
            raise ValueError("fld_mean must be positive")
        cmd.extend(["--fldMean", str(fld_mean)])
    if fld_sd is not None:
        if fld_sd <= 0:
            raise ValueError("fld_sd must be positive")
        cmd.extend(["--fldSD", str(fld_sd)])
    if min_score_fraction is not None:
        if not (0.0 <= min_score_fraction <= 1.0):
            raise ValueError("min_score_fraction must be between 0 and 1")
        cmd.extend(["--minScoreFraction", str(min_score_fraction)])
    if bandwidth is not None:
        if bandwidth <= 0:
            raise ValueError("bandwidth must be positive")
        cmd.extend(["--bandwidth", str(bandwidth)])
    if max_mmpextension is not None:
        if max_mmpextension <= 0:
            raise ValueError("max_mmpextension must be positive")
        cmd.extend(["--maxMMPExtension", str(max_mmpextension)])
    if ma is not None:
        if ma <= 0:
            raise ValueError("ma (match score) must be positive")
        cmd.extend(["--ma", str(ma)])
    if mp is not None:
        if mp >= 0:
            raise ValueError("mp (mismatch penalty) must be negative")
        cmd.extend(["--mp", str(mp)])
    if go is not None:
        if go <= 0:
            raise ValueError("go (gap open penalty) must be positive")
        cmd.extend(["--go", str(go)])
    if ge is not None:
        if ge <= 0:
            raise ValueError("ge (gap extension penalty) must be positive")
        cmd.extend(["--ge", str(ge)])
    if range_factorization_bins is not None:
        if range_factorization_bins <= 0:
            raise ValueError("range_factorization_bins must be positive")
        cmd.extend(["--rangeFactorizationBins", str(range_factorization_bins)])
    if use_em:
        cmd.append("--useEM")
    if vb_prior is not None:
        if vb_prior < 0:
            raise ValueError("vb_prior must be non-negative")
        cmd.extend(["--vbPrior", str(vb_prior)])
    if per_transcript_prior:
        cmd.append("--perTranscriptPrior")
    if num_bootstraps > 0:
        cmd.extend(["--numBootstraps", str(num_bootstraps)])
    if num_gibbs_samples > 0:
        cmd.extend(["--numGibbsSamples", str(num_gibbs_samples)])
    if seq_bias:
        cmd.append("--seqBias")
    if num_bias_samples is not None:
        if num_bias_samples <= 0:
            raise ValueError("num_bias_samples must be positive")
        cmd.extend(["--numBiasSamples", str(num_bias_samples)])
    if gc_bias:
        cmd.append("--gcBias")
    if pos_bias:
        cmd.append("--posBias")
    if bias_speed_samp <= 0:
        raise ValueError("bias_speed_samp must be positive")
    cmd.extend(["--biasSpeedSamp", str(bias_speed_samp)])
    if write_unmapped_names:
        cmd.append("--writeUnmappedNames")
    if write_mappings:
        if isinstance(write_mappings, bool):
            if write_mappings:
                # write to stdout
                cmd.append("--writeMappings")
        else:
            # write_mappings is a Path
            cmd.append(f"--writeMappings={str(write_mappings)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = [str(output_dir)]
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"Salmon quant failed with exit code {e.returncode}",
        }


if __name__ == '__main__':
    mcp.run()