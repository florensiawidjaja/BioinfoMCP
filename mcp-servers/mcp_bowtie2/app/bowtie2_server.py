from fastmcp import FastMCP
from typing import Optional, List
from pathlib import Path
import subprocess

mcp = FastMCP()


@mcp.tool()
def bowtie2_align(
    index_base: str,
    #mate1_files: Optional[List[str]] = None,
    #mate2_files: Optional[List[str]] = None,
    mate1_files: str = None,
    mate2_files: str = None,
    unpaired_files: Optional[List[str]] = None,
    interleaved: Optional[Path] = None,
    sra_accession: Optional[str] = None,
    bam_unaligned: Optional[Path] = None,
    sam_output: Optional[Path] = None,
    input_format_fastq: bool = True,
    tab5: bool = False,
    tab6: bool = False,
    qseq: bool = False,
    fasta: bool = False,
    one_seq_per_line: bool = False,
    kmer_fasta: Optional[Path] = None,
    kmer_int: Optional[int] = None,
    kmer_i: Optional[int] = None,
    reads_on_cmdline: Optional[List[str]] = None,
    skip_reads: int = 0,
    max_reads: Optional[int] = None,
    trim5: int = 0,
    trim3: int = 0,
    trim_to: Optional[str] = None,
    phred33: bool = False,
    phred64: bool = False,
    solexa_quals: bool = False,
    int_quals: bool = False,
    very_fast: bool = False,
    fast: bool = False,
    sensitive: bool = False,
    very_sensitive: bool = False,
    very_fast_local: bool = False,
    fast_local: bool = False,
    sensitive_local: bool = False,
    very_sensitive_local: bool = False,
    mismatches_seed: int = 0,
    seed_length: Optional[int] = None,
    seed_interval_func: Optional[str] = None,
    n_ceil_func: Optional[str] = None,
    dpad: int = 15,
    gbar: int = 4,
    ignore_quals: bool = False,
    nofw: bool = False,
    norc: bool = False,
    no_1mm_upfront: bool = False,
    end_to_end: bool = True,
    local: bool = False,
    match_bonus: int = 0,
    mp_max: int = 6,
    mp_min: int = 2,
    np_penalty: int = 1,
    rdg_open: int = 5,
    rdg_extend: int = 3,
    rfg_open: int = 5,
    rfg_extend: int = 3,
    score_min_func: Optional[str] = None,
    k: Optional[int] = None,
    a: bool = False,
    D: int = 15,
    R: int = 2,
    minins: int = 0,
    maxins: int = 500,
    fr: bool = True,
    rf: bool = False,
    ff: bool = False,
    no_mixed: bool = False,
    no_discordant: bool = False,
    dovetail: bool = False,
    no_contain: bool = False,
    no_overlap: bool = False,
    align_paired_reads: bool = False,
    preserve_tags: bool = False,
    quiet: bool = False,
    met_file: Optional[Path] = None,
    met_stderr: Optional[Path] = None,
    met_interval: int = 1,
    no_unal: bool = False,
    no_hd: bool = False,
    no_sq: bool = False,
    rg_id: Optional[str] = None,
    rg_fields: Optional[List[str]] = None,
    omit_sec_seq: bool = False,
    soft_clipped_unmapped_tlen: bool = False,
    sam_no_qname_trunc: bool = False,
    xeq: bool = False,
    sam_append_comment: bool = False,
    sam_opt_config: Optional[str] = None,
    offrate: Optional[int] = None,
    threads: int = 1,
    reorder: bool = False,
    mm: bool = False,
    qc_filter: bool = False,
    seed: int = 0,
    non_deterministic: bool = False,
) -> dict:
    """
    Bowtie2 aligner: aligns sequencing reads to a reference genome index and outputs SAM alignments.

    Parameters:
    - index_base: basename of the Bowtie2 index files.
    - mate1_files: A file containing mate 1 reads (comma-separated).
    - mate2_files: A file containing mate 2 reads (comma-separated).
    - unpaired_files: list of files containing unpaired reads (comma-separated).
    - interleaved: interleaved FASTQ file containing paired reads.
    - sra_accession: SRA accession to fetch reads from.
    - bam_unaligned: BAM file with unaligned reads.
    - sam_output: output SAM file path.
    - input_format_fastq: input reads are FASTQ (default True).
    - tab5, tab6, qseq, fasta, one_seq_per_line: input format flags.
    - kmer_fasta, kmer_int, kmer_i: k-mer extraction from fasta input.
    - reads_on_cmdline: reads given on command line.
    - skip_reads: skip first N reads.
    - max_reads: limit number of reads to align.
    - trim5, trim3: trim bases from 5' or 3' ends.
    - trim_to: trim reads exceeding length from 3' or 5'.
    - phred33, phred64, solexa_quals, int_quals: quality encoding options.
    - very_fast, fast, sensitive, very_sensitive: preset options for end-to-end mode.
    - very_fast_local, fast_local, sensitive_local, very_sensitive_local: preset options for local mode.
    - mismatches_seed: number of mismatches allowed in seed.
    - seed_length: seed substring length.
    - seed_interval_func: function governing seed interval.
    - n_ceil_func: function governing max ambiguous chars.
    - dpad, gbar: gap padding and disallow gap near ends.
    - ignore_quals: ignore quality values in mismatch penalty.
    - nofw, norc: disable forward or reverse strand alignment.
    - no_1mm_upfront: disable 1-mismatch end-to-end search upfront.
    - end_to_end, local: alignment mode flags.
    - match_bonus: match bonus in local mode.
    - mp_max, mp_min: mismatch penalties max and min.
    - np_penalty: penalty for ambiguous characters.
    - rdg_open, rdg_extend: read gap open and extend penalties.
    - rfg_open, rfg_extend: reference gap open and extend penalties.
    - score_min_func: minimum score function.
    - k: max number of distinct valid alignments to report.
    - a: report all valid alignments.
    - D, R: effort options controlling search.
    - minins, maxins: min and max fragment length for paired-end.
    - fr, rf, ff: mate orientation flags.
    - no_mixed, no_discordant: disable mixed or discordant alignments.
    - dovetail, no_contain, no_overlap: paired-end overlap behavior.
    - align_paired_reads: align paired BAM reads.
    - preserve_tags: preserve BAM tags.
    - quiet: suppress non-error output.
    - met_file, met_stderr, met_interval: metrics output options.
    - no_unal, no_hd, no_sq: suppress SAM output lines.
    - rg_id, rg_fields: read group header and fields.
    - omit_sec_seq: omit SEQ and QUAL in secondary alignments.
    - soft_clipped_unmapped_tlen: consider soft-clipped bases unmapped in TLEN.
    - sam_no_qname_trunc: disable truncation of read names.
    - xeq: use '='/'X' in CIGAR.
    - sam_append_comment: append FASTA/FASTQ comment to SAM.
    - sam_opt_config: configure SAM optional fields.
    - offrate: override index offrate.
    - threads: number of parallel threads.
    - reorder: guarantee output order matches input.
    - mm: use memory-mapped I/O for index.
    - qc_filter: filter reads failing QSEQ filter.
    - seed: seed for pseudo-random generator.
    - non_deterministic: use current time for random seed.

    Returns:
    dict with keys: command_executed, stdout, stderr, output_files (list).
    """
    import shlex

    # Validate mutually exclusive options
    if end_to_end and local:
        raise ValueError("Options --end-to-end and --local are mutually exclusive.")
    if k is not None and a:
        raise ValueError("Options -k and -a are mutually exclusive.")
    if trim_to is not None and (trim5 > 0 or trim3 > 0):
        raise ValueError("--trim-to and -3/-5 are mutually exclusive.")
    if phred33 and phred64:
        raise ValueError("--phred33 and --phred64 are mutually exclusive.")
    if mate1_files is not None and interleaved is not None:
        raise ValueError("Cannot specify both -1 and --interleaved.")
    if mate2_files is not None and interleaved is not None:
        raise ValueError("Cannot specify both -2 and --interleaved.")
    if (mate1_files is None) != (mate2_files is None):
        raise ValueError("Both -1 and -2 must be specified together for paired-end reads.")

    # Validate input files exist
    def check_files_exist(files: Optional[List[str]], param_name: str):
        if files:
            for f in files:
                if f != "-" and not Path(f).exists():
                    raise FileNotFoundError(f"Input file '{f}' specified in {param_name} does not exist.")

    #check_files_exist(mate1_files, "-1")
    #check_files_exist(mate2_files, "-2")
    check_files_exist(unpaired_files, "-U")
    if interleaved is not None and not interleaved.exists():
        raise FileNotFoundError(f"Interleaved file '{interleaved}' does not exist.")
    if bam_unaligned is not None and not bam_unaligned.exists():
        raise FileNotFoundError(f"BAM file '{bam_unaligned}' does not exist.")
    if kmer_fasta is not None and not kmer_fasta.exists():
        raise FileNotFoundError(f"K-mer fasta file '{kmer_fasta}' does not exist.")
    if sam_output is not None:
        sam_output = Path(sam_output)
        if sam_output.exists() and not sam_output.is_file():
            raise ValueError(f"Output SAM path '{sam_output}' exists and is not a file.")

    # Build command
    cmd = ["bowtie2"]

    # Index base (required)
    cmd.extend(["-x", index_base])

    # Input reads
    if mate1_files is not None and mate2_files is not None:
        cmd.extend(["-1", mate1_files])
        cmd.extend(["-2", mate2_files])
        #cmd.extend(["-1", ",".join(mate1_files)])
        #cmd.extend(["-2", ",".join(mate2_files)])
    elif unpaired_files is not None:
        cmd.extend(["-U", ",".join(unpaired_files)])
    elif interleaved is not None:
        cmd.extend(["--interleaved", str(interleaved)])
    elif sra_accession is not None:
        cmd.extend(["--sra-acc", sra_accession])
    elif bam_unaligned is not None:
        cmd.extend(["-b", str(bam_unaligned)])
    elif reads_on_cmdline is not None:
        # -c option: reads given on command line
        cmd.extend(["-c"])
        cmd.extend(reads_on_cmdline)
    elif kmer_fasta is not None and kmer_int is not None and kmer_i is not None:
        cmd.extend(["-F", f"{kmer_int},i:{kmer_i}"])
        cmd.append(str(kmer_fasta))
    else:
        raise ValueError("No input reads specified. Provide -1/-2, -U, --interleaved, --sra-acc, -b, -c, or -F options.")

    # Output SAM
    if sam_output is not None:
        cmd.extend(["-S", str(sam_output)])

    # Input format options
    if input_format_fastq:
        cmd.append("-q")
    if tab5:
        cmd.append("--tab5")
    if tab6:
        cmd.append("--tab6")
    if qseq:
        cmd.append("--qseq")
    if fasta:
        cmd.append("-f")
    if one_seq_per_line:
        cmd.append("-r")

    # Skip and limit reads
    if skip_reads > 0:
        cmd.extend(["-s", str(skip_reads)])
    if max_reads is not None:
        cmd.extend(["-u", str(max_reads)])

    # Trimming
    if trim5 > 0:
        cmd.extend(["-5", str(trim5)])
    if trim3 > 0:
        cmd.extend(["-3", str(trim3)])
    if trim_to is not None:
        # trim_to format: [3:|5:]<int>
        cmd.extend(["--trim-to", trim_to])

    # Quality encoding
    if phred33:
        cmd.append("--phred33")
    if phred64:
        cmd.append("--phred64")
    if solexa_quals:
        cmd.append("--solexa-quals")
    if int_quals:
        cmd.append("--int-quals")

    # Presets
    if very_fast:
        cmd.append("--very-fast")
    if fast:
        cmd.append("--fast")
    if sensitive:
        cmd.append("--sensitive")
    if very_sensitive:
        cmd.append("--very-sensitive")
    if very_fast_local:
        cmd.append("--very-fast-local")
    if fast_local:
        cmd.append("--fast-local")
    if sensitive_local:
        cmd.append("--sensitive-local")
    if very_sensitive_local:
        cmd.append("--very-sensitive-local")

    # Alignment options
    if mismatches_seed not in (0, 1):
        raise ValueError("-N must be 0 or 1")
    cmd.extend(["-N", str(mismatches_seed)])

    if seed_length is not None:
        cmd.extend(["-L", str(seed_length)])

    if seed_interval_func is not None:
        cmd.extend(["-i", seed_interval_func])

    if n_ceil_func is not None:
        cmd.extend(["--n-ceil", n_ceil_func])

    cmd.extend(["--dpad", str(dpad)])
    cmd.extend(["--gbar", str(gbar)])

    if ignore_quals:
        cmd.append("--ignore-quals")
    if nofw:
        cmd.append("--nofw")
    if norc:
        cmd.append("--norc")
    if no_1mm_upfront:
        cmd.append("--no-1mm-upfront")

    if end_to_end:
        cmd.append("--end-to-end")
    if local:
        cmd.append("--local")

    cmd.extend(["--ma", str(match_bonus)])
    cmd.extend(["--mp", f"{mp_max},{mp_min}"])
    cmd.extend(["--np", str(np_penalty)])
    cmd.extend(["--rdg", f"{rdg_open},{rdg_extend}"])
    cmd.extend(["--rfg", f"{rfg_open},{rfg_extend}"])

    if score_min_func is not None:
        cmd.extend(["--score-min", score_min_func])

    # Reporting options
    if k is not None:
        if k < 1:
            raise ValueError("-k must be >= 1")
        cmd.extend(["-k", str(k)])
    if a:
        cmd.append("-a")

    # Effort options
    cmd.extend(["-D", str(D)])
    cmd.extend(["-R", str(R)])

    # Paired-end options
    cmd.extend(["-I", str(minins)])
    cmd.extend(["-X", str(maxins)])

    if fr:
        cmd.append("--fr")
    if rf:
        cmd.append("--rf")
    if ff:
        cmd.append("--ff")

    if no_mixed:
        cmd.append("--no-mixed")
    if no_discordant:
        cmd.append("--no-discordant")
    if dovetail:
        cmd.append("--dovetail")
    if no_contain:
        cmd.append("--no-contain")
    if no_overlap:
        cmd.append("--no-overlap")

    # BAM options
    if align_paired_reads:
        cmd.append("--align-paired-reads")
    if preserve_tags:
        cmd.append("--preserve-tags")

    # Output options
    if quiet:
        cmd.append("--quiet")
    if met_file is not None:
        cmd.extend(["--met-file", str(met_file)])
    if met_stderr is not None:
        cmd.extend(["--met-stderr", str(met_stderr)])
    cmd.extend(["--met", str(met_interval)])

    if no_unal:
        cmd.append("--no-unal")
    if no_hd:
        cmd.append("--no-hd")
    if no_sq:
        cmd.append("--no-sq")

    if rg_id is not None:
        cmd.extend(["--rg-id", rg_id])
    if rg_fields is not None:
        for field in rg_fields:
            cmd.extend(["--rg", field])

    if omit_sec_seq:
        cmd.append("--omit-sec-seq")
    if soft_clipped_unmapped_tlen:
        cmd.append("--soft-clipped-unmapped-tlen")
    if sam_no_qname_trunc:
        cmd.append("--sam-no-qname-trunc")
    if xeq:
        cmd.append("--xeq")
    if sam_append_comment:
        cmd.append("--sam-append-comment")
    if sam_opt_config is not None:
        cmd.extend(["--sam-opt-config", sam_opt_config])

    if offrate is not None:
        cmd.extend(["-o", str(offrate)])

    cmd.extend(["-p", str(threads)])

    if reorder:
        cmd.append("--reorder")
    if mm:
        cmd.append("--mm")
    if qc_filter:
        cmd.append("--qc-filter")

    cmd.extend(["--seed", str(seed)])

    if non_deterministic:
        cmd.append("--non-deterministic")

    # Run command
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
            "error": f"Bowtie2 alignment failed with return code {e.returncode}",
            "output_files": [],
        }

    output_files = []
    if sam_output is not None:
        output_files.append(str(sam_output))

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


@mcp.tool()
def bowtie2_build(
    reference_in: List[str],
    index_base: str,
    fasta: bool = False,
    sequences_on_cmdline: bool = False,
    large_index: bool = False,
    noauto: bool = False,
    packed: bool = False,
    bmax: Optional[int] = None,
    bmaxdivn: Optional[int] = None,
    dcv: Optional[int] = None,
    nodc: bool = False,
    noref: bool = False,
    justref: bool = False,
    offrate: Optional[int] = None,
    ftabchars: Optional[int] = None,
    seed: Optional[int] = None,
    cutoff: Optional[int] = None,
    quiet: bool = False,
    threads: int = 1,
) -> dict:
    """
    Bowtie2-build indexer: builds a Bowtie2 index from reference sequences.

    Parameters:
    - reference_in: list of FASTA files or sequences (if -c).
    - index_base: basename for output index files.
    - fasta: input files are FASTA format.
    - sequences_on_cmdline: sequences given on command line (-c).
    - large_index: force building large index.
    - noauto: disable automatic parameter selection.
    - packed: use packed DNA representation.
    - bmax: max suffixes per block.
    - bmaxdivn: max suffixes per block as fraction of reference length.
    - dcv: period for difference-cover sample.
    - nodc: disable difference-cover sample.
    - noref: do not build bitpacked reference portions.
    - justref: build only bitpacked reference portions.
    - offrate: override offrate.
    - ftabchars: ftab lookup table size.
    - seed: seed for random number generator.
    - cutoff: index only first N bases.
    - quiet: suppress output except errors.
    - threads: number of threads.

    Returns:
    dict with keys: command_executed, stdout, stderr, output_files (list).
    """
    import shlex

    # Validate input files if not sequences on cmdline
    if not sequences_on_cmdline:
        for f in reference_in:
            if not Path(f).exists():
                raise FileNotFoundError(f"Reference input file '{f}' does not exist.")

    cmd = ["bowtie2-build"]

    if fasta:
        cmd.append("-f")
    if sequences_on_cmdline:
        cmd.append("-c")
    if large_index:
        cmd.append("--large-index")
    if noauto:
        cmd.append("-a")
    if packed:
        cmd.append("-p")
    if bmax is not None:
        cmd.extend(["--bmax", str(bmax)])
    if bmaxdivn is not None:
        cmd.extend(["--bmaxdivn", str(bmaxdivn)])
    if dcv is not None:
        cmd.extend(["--dcv", str(dcv)])
    if nodc:
        cmd.append("--nodc")
    if noref:
        cmd.append("-r")
    if justref:
        cmd.append("-3")
    if offrate is not None:
        cmd.extend(["-o", str(offrate)])
    if ftabchars is not None:
        cmd.extend(["-t", str(ftabchars)])
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    if cutoff is not None:
        cmd.extend(["--cutoff", str(cutoff)])
    if quiet:
        cmd.append("-q")
    cmd.extend(["--threads", str(threads)])

    # Add reference input and index base
    if sequences_on_cmdline:
        # reference_in are sequences separated by commas
        cmd.append(",".join(reference_in))
    else:
        # reference_in are files separated by commas
        cmd.append(",".join(reference_in))
    cmd.append(index_base)

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
            "error": f"bowtie2-build failed with return code {e.returncode}",
            "output_files": [],
        }

    # Output files: 6 files with suffixes .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2
    suffixes = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    if large_index:
        suffixes = [s.replace(".bt2", ".bt2l") for s in suffixes]

    output_files = [f"{index_base}{s}" for s in suffixes]

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


@mcp.tool()
def bowtie2_inspect(
    index_base: str,
    across: int = 60,
    names: bool = False,
    summary: bool = False,
    output: Optional[Path] = None,
    verbose: bool = False,
) -> dict:
    """
    Bowtie2-inspect index inspector: extracts information from a Bowtie2 index.

    Parameters:
    - index_base: basename of the index to inspect.
    - across: number of bases per line in FASTA output (default 60).
    - names: print reference sequence names only.
    - summary: print summary of index.
    - output: output file path (default stdout).
    - verbose: print verbose output.

    Returns:
    dict with keys: command_executed, stdout, stderr, output_files (list).
    """
    import shlex

    cmd = ["bowtie2-inspect"]

    cmd.extend(["-a", str(across)])

    if names:
        cmd.append("-n")
    if summary:
        cmd.append("-s")
    if output is not None:
        cmd.extend(["-o", str(output)])
    if verbose:
        cmd.append("-v")

    cmd.append(index_base)

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
            "error": f"bowtie2-inspect failed with return code {e.returncode}",
            "output_files": [],
        }

    output_files = []
    if output is not None:
        output_files.append(str(output))

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


if __name__ == '__main__':
    mcp.run()