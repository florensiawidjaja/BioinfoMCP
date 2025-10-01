from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


def _validate_func_option(func: str) -> None:
    """Validate function option format F,B,A where F in {C,L,S,G} and B,A are floats."""
    parts = func.split(",")
    if len(parts) != 3:
        raise ValueError(f"Function option must have 3 parts separated by commas: {func}")
    F, B, A = parts
    if F not in {"C", "L", "S", "G"}:
        raise ValueError(f"Function type must be one of C,L,S,G but got {F}")
    try:
        float(B)
        float(A)
    except ValueError:
        raise ValueError(f"Constant term and coefficient must be floats: {B}, {A}")


def _validate_int_pair(value: str, name: str) -> (int, int):
    """Validate a comma-separated pair of integers."""
    parts = value.split(",")
    if len(parts) != 2:
        raise ValueError(f"{name} must be two comma-separated integers")
    try:
        i1 = int(parts[0])
        i2 = int(parts[1])
    except ValueError:
        raise ValueError(f"{name} values must be integers")
    return i1, i2


@mcp.tool()
def hisat2_align(
    index_basename: str,
    mate1: Optional[str] = None,
    mate2: Optional[str] = None,
    unpaired: Optional[str] = None,
    sra_acc: Optional[str] = None,
    sam_output: Optional[str] = None,
    # Input format options
    fastq: bool = True,
    qseq: bool = False,
    fasta: bool = False,
    one_seq_per_line: bool = False,
    reads_on_cmdline: bool = False,
    skip: int = 0,
    upto: int = 0,
    trim5: int = 0,
    trim3: int = 0,
    phred33: bool = False,
    phred64: bool = False,
    solexa_quals: bool = False,
    int_quals: bool = False,
    # Alignment options
    n_ceil: str = "L,0,0.15",
    ignore_quals: bool = False,
    nofw: bool = False,
    norc: bool = False,
    # Scoring options
    mp: str = "6,2",
    sp: str = "2,1",
    no_softclip: bool = False,
    np: int = 1,
    rdg: str = "5,3",
    rfg: str = "5,3",
    score_min: str = "L,0,-0.2",
    # Spliced alignment options
    pen_cansplice: int = 0,
    pen_noncansplice: int = 12,
    pen_canintronlen: str = "G,-8,1",
    pen_noncanintronlen: str = "G,-8,1",
    min_intronlen: int = 20,
    max_intronlen: int = 500000,
    known_splicesite_infile: Optional[str] = None,
    novel_splicesite_outfile: Optional[str] = None,
    novel_splicesite_infile: Optional[str] = None,
    no_temp_splicesite: bool = False,
    no_spliced_alignment: bool = False,
    rna_strandness: Optional[str] = None,
    tmo: bool = False,
    dta: bool = False,
    dta_cufflinks: bool = False,
    avoid_pseudogene: bool = False,
    no_templatelen_adjustment: bool = False,
    # Reporting options
    k: int = 5,
    max_seeds: int = 10,
    all_alignments: bool = False,
    secondary: bool = False,
    # Paired-end options
    minins: int = 0,
    maxins: int = 500,
    fr: bool = True,
    rf: bool = False,
    ff: bool = False,
    no_mixed: bool = False,
    no_discordant: bool = False,
    # Output options
    time: bool = False,
    un: Optional[str] = None,
    un_gz: Optional[str] = None,
    un_bz2: Optional[str] = None,
    al: Optional[str] = None,
    al_gz: Optional[str] = None,
    al_bz2: Optional[str] = None,
    un_conc: Optional[str] = None,
    un_conc_gz: Optional[str] = None,
    un_conc_bz2: Optional[str] = None,
    al_conc: Optional[str] = None,
    al_conc_gz: Optional[str] = None,
    al_conc_bz2: Optional[str] = None,
    quiet: bool = False,
    summary_file: Optional[str] = None,
    new_summary: bool = False,
    met_file: Optional[str] = None,
    met_stderr: bool = False,
    met: int = 1,
    # SAM options
    no_unal: bool = False,
    no_hd: bool = False,
    no_sq: bool = False,
    rg_id: Optional[str] = None,
    rg: Optional[List[str]] = None,
    remove_chrname: bool = False,
    add_chrname: bool = False,
    omit_sec_seq: bool = False,
    # Performance options
    offrate: Optional[int] = None,
    threads: int = 1,
    reorder: bool = False,
    mm: bool = False,
    # Other options
    qc_filter: bool = False,
    seed: int = 0,
    non_deterministic: bool = False,
):
    """
    Run HISAT2 alignment with comprehensive options.

    Parameters:
    - index_basename: Basename of the HISAT2 index files.
    - mate1: Comma-separated list of mate 1 files.
    - mate2: Comma-separated list of mate 2 files.
    - unpaired: Comma-separated list of unpaired read files.
    - sra_acc: Comma-separated list of SRA accession numbers.
    - sam_output: Output SAM file path.
    - fastq, qseq, fasta, one_seq_per_line, reads_on_cmdline: Input format flags.
    - skip, upto, trim5, trim3: Read processing options.
    - phred33, phred64, solexa_quals, int_quals: Quality encoding options.
    - n_ceil: Function string for max ambiguous chars allowed.
    - ignore_quals, nofw, norc: Alignment behavior flags.
    - mp, sp, no_softclip, np, rdg, rfg, score_min: Scoring options.
    - pen_cansplice, pen_noncansplice, pen_canintronlen, pen_noncanintronlen: Splice penalties.
    - min_intronlen, max_intronlen: Intron length constraints.
    - known_splicesite_infile, novel_splicesite_outfile, novel_splicesite_infile: Splice site files.
    - no_temp_splicesite, no_spliced_alignment: Spliced alignment flags.
    - rna_strandness: Strand-specific info.
    - tmo, dta, dta_cufflinks, avoid_pseudogene, no_templatelen_adjustment: RNA-seq options.
    - k, max_seeds, all_alignments, secondary: Reporting and alignment count options.
    - minins, maxins, fr, rf, ff, no_mixed, no_discordant: Paired-end options.
    - time: Print wall-clock time.
    - un, un_gz, un_bz2, al, al_gz, al_bz2, un_conc, un_conc_gz, un_conc_bz2, al_conc, al_conc_gz, al_conc_bz2: Output read files.
    - quiet, summary_file, new_summary, met_file, met_stderr, met: Output and metrics options.
    - no_unal, no_hd, no_sq, rg_id, rg, remove_chrname, add_chrname, omit_sec_seq: SAM output options.
    - offrate, threads, reorder, mm: Performance options.
    - qc_filter, seed, non_deterministic: Other options.
    """
    # Validate index basename path (no extension)
    if not index_basename:
        raise ValueError("index_basename must be specified")
    # Validate input files if provided
    def _check_files_csv(csv: Optional[str], name: str):
        if csv:
            for f in csv.split(","):
                if f != "-" and not Path(f).exists():
                    raise FileNotFoundError(f"{name} file does not exist: {f}")

    _check_files_csv(mate1, "mate1")
    _check_files_csv(mate2, "mate2")
    _check_files_csv(unpaired, "unpaired")
    _check_files_csv(known_splicesite_infile, "known_splicesite_infile")
    _check_files_csv(novel_splicesite_infile, "novel_splicesite_infile")

    # Validate function options
    _validate_func_option(n_ceil)
    _validate_func_option(score_min)
    _validate_func_option(pen_canintronlen)
    _validate_func_option(pen_noncanintronlen)

    # Validate comma-separated integer pairs
    mp_mx, mp_mn = _validate_int_pair(mp, "mp")
    sp_mx, sp_mn = _validate_int_pair(sp, "sp")
    rdg_open, rdg_extend = _validate_int_pair(rdg, "rdg")
    rfg_open, rfg_extend = _validate_int_pair(rfg, "rfg")

    # Validate strandness
    if rna_strandness is not None:
        if rna_strandness not in {"F", "R", "FR", "RF"}:
            raise ValueError("rna_strandness must be one of F, R, FR, RF")

    # Validate paired-end orientation flags
    if sum([fr, rf, ff]) > 1:
        raise ValueError("Only one of --fr, --rf, --ff can be specified")

    # Validate threads
    if threads < 1:
        raise ValueError("threads must be >= 1")

    # Validate skip, upto, trim5, trim3
    if skip < 0:
        raise ValueError("skip must be >= 0")
    if upto < 0:
        raise ValueError("upto must be >= 0")
    if trim5 < 0:
        raise ValueError("trim5 must be >= 0")
    if trim3 < 0:
        raise ValueError("trim3 must be >= 0")

    # Validate min_intronlen and max_intronlen
    if min_intronlen < 0:
        raise ValueError("min_intronlen must be >= 0")
    if max_intronlen < min_intronlen:
        raise ValueError("max_intronlen must be >= min_intronlen")

    # Validate k and max_seeds
    if k < 1:
        raise ValueError("k must be >= 1")
    if max_seeds < 1:
        raise ValueError("max_seeds must be >= 1")

    # Validate offrate if specified
    if offrate is not None and offrate < 1:
        raise ValueError("offrate must be >= 1")

    # Validate seed
    if seed < 0:
        raise ValueError("seed must be >= 0")

    # Build command line
    cmd = ["hisat2"]

    # Index basename
    cmd += ["-x", index_basename]

    # Input reads
    if mate1 and mate2:
        cmd += ["-1", mate1, "-2", mate2]
    elif unpaired:
        cmd += ["-U", unpaired]
    elif sra_acc:
        cmd += ["--sra-acc", sra_acc]
    else:
        raise ValueError("Must specify either mate1 and mate2, or unpaired, or sra_acc")

    # Output SAM file
    if sam_output:
        cmd += ["-S", sam_output]

    # Input format options
    if fastq:
        cmd.append("-q")
    if qseq:
        cmd.append("--qseq")
    if fasta:
        cmd.append("-f")
    if one_seq_per_line:
        cmd.append("-r")
    if reads_on_cmdline:
        cmd.append("-c")

    # Read processing
    if skip > 0:
        cmd += ["-s", str(skip)]
    if upto > 0:
        cmd += ["-u", str(upto)]
    if trim5 > 0:
        cmd += ["-5", str(trim5)]
    if trim3 > 0:
        cmd += ["-3", str(trim3)]

    # Quality encoding
    if phred33:
        cmd.append("--phred33")
    if phred64:
        cmd.append("--phred64")
    if solexa_quals:
        cmd.append("--solexa-quals")
    if int_quals:
        cmd.append("--int-quals")

    # Alignment options
    if n_ceil != "L,0,0.15":
        cmd += ["--n-ceil", n_ceil]
    if ignore_quals:
        cmd.append("--ignore-quals")
    if nofw:
        cmd.append("--nofw")
    if norc:
        cmd.append("--norc")

    # Scoring options
    if mp != "6,2":
        cmd += ["--mp", mp]
    if sp != "2,1":
        cmd += ["--sp", sp]
    if no_softclip:
        cmd.append("--no-softclip")
    if np != 1:
        cmd += ["--np", str(np)]
    if rdg != "5,3":
        cmd += ["--rdg", rdg]
    if rfg != "5,3":
        cmd += ["--rfg", rfg]
    if score_min != "L,0,-0.2":
        cmd += ["--score-min", score_min]

    # Spliced alignment options
    if pen_cansplice != 0:
        cmd += ["--pen-cansplice", str(pen_cansplice)]
    if pen_noncansplice != 12:
        cmd += ["--pen-noncansplice", str(pen_noncansplice)]
    if pen_canintronlen != "G,-8,1":
        cmd += ["--pen-canintronlen", pen_canintronlen]
    if pen_noncanintronlen != "G,-8,1":
        cmd += ["--pen-noncanintronlen", pen_noncanintronlen]
    if min_intronlen != 20:
        cmd += ["--min-intronlen", str(min_intronlen)]
    if max_intronlen != 500000:
        cmd += ["--max-intronlen", str(max_intronlen)]
    if known_splicesite_infile:
        cmd += ["--known-splicesite-infile", known_splicesite_infile]
    if novel_splicesite_outfile:
        cmd += ["--novel-splicesite-outfile", novel_splicesite_outfile]
    if novel_splicesite_infile:
        cmd += ["--novel-splicesite-infile", novel_splicesite_infile]
    if no_temp_splicesite:
        cmd.append("--no-temp-splicesite")
    if no_spliced_alignment:
        cmd.append("--no-spliced-alignment")
    if rna_strandness:
        cmd += ["--rna-strandness", rna_strandness]
    if tmo:
        cmd.append("--tmo")
    if dta:
        cmd.append("--dta")
    if dta_cufflinks:
        cmd.append("--dta-cufflinks")
    if avoid_pseudogene:
        cmd.append("--avoid-pseudogene")
    if no_templatelen_adjustment:
        cmd.append("--no-templatelen-adjustment")

    # Reporting options
    if k != 5:
        cmd += ["-k", str(k)]
    if max_seeds != 10:
        cmd += ["--max-seeds", str(max_seeds)]
    if all_alignments:
        cmd.append("-a")
    if secondary:
        cmd.append("--secondary")

    # Paired-end options
    if minins != 0:
        cmd += ["-I", str(minins)]
    if maxins != 500:
        cmd += ["-X", str(maxins)]
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

    # Output options
    if time:
        cmd.append("-t")
    if un:
        cmd += ["--un", un]
    if un_gz:
        cmd += ["--un-gz", un_gz]
    if un_bz2:
        cmd += ["--un-bz2", un_bz2]
    if al:
        cmd += ["--al", al]
    if al_gz:
        cmd += ["--al-gz", al_gz]
    if al_bz2:
        cmd += ["--al-bz2", al_bz2]
    if un_conc:
        cmd += ["--un-conc", un_conc]
    if un_conc_gz:
        cmd += ["--un-conc-gz", un_conc_gz]
    if un_conc_bz2:
        cmd += ["--un-conc-bz2", un_conc_bz2]
    if al_conc:
        cmd += ["--al-conc", al_conc]
    if al_conc_gz:
        cmd += ["--al-conc-gz", al_conc_gz]
    if al_conc_bz2:
        cmd += ["--al-conc-bz2", al_conc_bz2]
    if quiet:
        cmd.append("--quiet")
    if summary_file:
        cmd += ["--summary-file", summary_file]
    if new_summary:
        cmd.append("--new-summary")
    if met_file:
        cmd += ["--met-file", met_file]
    if met_stderr:
        cmd.append("--met-stderr")
    if met != 1:
        cmd += ["--met", str(met)]

    # SAM options
    if no_unal:
        cmd.append("--no-unal")
    if no_hd:
        cmd.append("--no-hd")
    if no_sq:
        cmd.append("--no-sq")
    if rg_id:
        cmd += ["--rg-id", rg_id]
    if rg:
        for rg_field in rg:
            cmd += ["--rg", rg_field]
    if remove_chrname:
        cmd.append("--remove-chrname")
    if add_chrname:
        cmd.append("--add-chrname")
    if omit_sec_seq:
        cmd.append("--omit-sec-seq")

    # Performance options
    if offrate is not None:
        cmd += ["-o", str(offrate)]
    if threads != 1:
        cmd += ["-p", str(threads)]
    if reorder:
        cmd.append("--reorder")
    if mm:
        cmd.append("--mm")

    # Other options
    if qc_filter:
        cmd.append("--qc-filter")
    if seed != 0:
        cmd += ["--seed", str(seed)]
    if non_deterministic:
        cmd.append("--non-deterministic")

    # Run command
    try:
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
        stdout = completed.stdout
        stderr = completed.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"hisat2 failed with exit code {e.returncode}",
            "output_files": []
        }

    # Collect output files
    output_files = []
    if sam_output:
        output_files.append(str(Path(sam_output).resolve()))
    if un:
        output_files.append(str(Path(un).resolve()))
    if un_gz:
        output_files.append(str(Path(un_gz).resolve()))
    if un_bz2:
        output_files.append(str(Path(un_bz2).resolve()))
    if al:
        output_files.append(str(Path(al).resolve()))
    if al_gz:
        output_files.append(str(Path(al_gz).resolve()))
    if al_bz2:
        output_files.append(str(Path(al_bz2).resolve()))
    if un_conc:
        output_files.append(str(Path(un_conc).resolve()))
    if un_conc_gz:
        output_files.append(str(Path(un_conc_gz).resolve()))
    if un_conc_bz2:
        output_files.append(str(Path(un_conc_bz2).resolve()))
    if al_conc:
        output_files.append(str(Path(al_conc).resolve()))
    if al_conc_gz:
        output_files.append(str(Path(al_conc_gz).resolve()))
    if al_conc_bz2:
        output_files.append(str(Path(al_conc_bz2).resolve()))
    if summary_file:
        output_files.append(str(Path(summary_file).resolve()))
    if met_file:
        output_files.append(str(Path(met_file).resolve()))
    if known_splicesite_infile:
        output_files.append(str(Path(known_splicesite_infile).resolve()))
    if novel_splicesite_outfile:
        output_files.append(str(Path(novel_splicesite_outfile).resolve()))
    if novel_splicesite_infile:
        output_files.append(str(Path(novel_splicesite_infile).resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()