from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


@mcp.tool()
def trim_galore(
    input_files: List[Path],
    paired: bool = False,
    quality: int = 20,
    phred64: bool = False,
    fastqc: bool = False,
    fastqc_args: Optional[str] = None,
    adapter: Optional[str] = None,
    adapter2: Optional[str] = None,
    illumina: bool = False,
    nextera: bool = False,
    small_rna: bool = False,
    length: int = 20,
    max_length: Optional[int] = None,
    stringency: int = 1,
    error_rate: float = 0.1,
    gzip_output: bool = True,
    dont_gzip: bool = False,
    retain_unpaired: bool = False,
    max_n: Optional[int] = None,
    trim_n: bool = False,
    output_dir: Optional[Path] = None,
    no_report_file: bool = False,
    suppress_warn: bool = False,
    clip_R1: int = 0,
    clip_R2: int = 0,
    three_prime_clip_R1: int = 0,
    three_prime_clip_R2: int = 0,
    two_colour: bool = False,
    nextseq_trim: Optional[int] = None,
    path_to_cutadapt: Optional[Path] = None,
    basename: Optional[str] = None,
    cores: int = 1,
    demux: Optional[Path] = None,
    hardtrim5: Optional[int] = None,
    hardtrim3: Optional[int] = None,
    clock: bool = False,
    dual_index: bool = False,
    rrbs: bool = False,
    non_directional: bool = False,
    keep: bool = False,
    trim1: bool = False,
    length_1: Optional[int] = None,
    length_2: Optional[int] = None,
) -> dict:
    """
    Run trim_galore for adapter and quality trimming of sequencing reads.

    Parameters:
    - input_files: List of input FASTQ files (single or paired-end).
    - paired: If True, treat input as paired-end reads.
    - quality: Quality Phred score cutoff for trimming (default 20).
    - phred64: Use Phred+64 quality encoding instead of Phred+33.
    - fastqc: Run FastQC after trimming.
    - fastqc_args: Additional arguments to pass to FastQC.
    - adapter: Adapter sequence to trim (overrides auto-detection).
    - adapter2: Adapter sequence for second read in paired-end mode.
    - illumina: Use Illumina TruSeq adapter sequence.
    - nextera: Use Nextera adapter sequence.
    - small_rna: Use small RNA adapter sequence.
    - length: Minimum required sequence length after trimming (default 20).
    - max_length: Maximum length of reads to keep.
    - stringency: Minimum overlap length for adapter trimming (default 1).
    - error_rate: Maximum allowed error rate for adapter trimming (default 0.1).
    - gzip_output: Compress output files with gzip (default True).
    - dont_gzip: Do not gzip output files (overrides gzip_output).
    - retain_unpaired: Retain unpaired reads in paired-end mode.
    - max_n: Maximum number of Ns allowed in a read.
    - trim_n: Trim Ns from ends of reads.
    - output_dir: Directory to write output files.
    - no_report_file: Do not write report file.
    - suppress_warn: Suppress warnings.
    - clip_R1: Number of bases to clip from 5' end of read 1.
    - clip_R2: Number of bases to clip from 5' end of read 2.
    - three_prime_clip_R1: Number of bases to clip from 3' end of read 1.
    - three_prime_clip_R2: Number of bases to clip from 3' end of read 2.
    - two_colour: Use two-colour chemistry trimming (NextSeq).
    - nextseq_trim: Trim bases from 3' end for NextSeq reads.
    - path_to_cutadapt: Path to cutadapt executable.
    - basename: Preferred basename for output files.
    - cores: Number of cores to use (default 1).
    - demux: Barcode file for demultiplexing.
    - hardtrim5: Hard trim N bases from 5' end.
    - hardtrim3: Hard trim N bases from 3' end.
    - clock: Use clock UMI trimming mode.
    - dual_index: Use dual index trimming mode.
    - rrbs: Use RRBS mode.
    - non_directional: Use non-directional RRBS mode.
    - keep: Keep temporary files.
    - trim1: Trim first base from reads.
    - length_1: Minimum length for read 1.
    - length_2: Minimum length for read 2.

    Returns:
    Dictionary with keys: command_executed, stdout, stderr, output_files.
    """
    # Validate input files
    if not input_files:
        raise ValueError("At least one input file must be provided.")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file does not exist: {f}")

    # Validate cores
    if cores < 1:
        raise ValueError("cores must be >= 1")

    # Validate quality
    if quality < 0:
        raise ValueError("quality must be non-negative")

    # Validate stringency
    if stringency < 0:
        raise ValueError("stringency must be non-negative")

    # Validate error_rate
    if not (0.0 <= error_rate <= 1.0):
        raise ValueError("error_rate must be between 0 and 1")

    # Validate length
    if length < 0:
        raise ValueError("length must be non-negative")

    # Validate clip values
    for clip_val, name in [
        (clip_R1, "clip_R1"),
        (clip_R2, "clip_R2"),
        (three_prime_clip_R1, "three_prime_clip_R1"),
        (three_prime_clip_R2, "three_prime_clip_R2"),
    ]:
        if clip_val < 0:
            raise ValueError(f"{name} must be non-negative")

    # Validate hardtrim values
    if hardtrim5 is not None and hardtrim5 < 0:
        raise ValueError("hardtrim5 must be non-negative")
    if hardtrim3 is not None and hardtrim3 < 0:
        raise ValueError("hardtrim3 must be non-negative")

    # Validate length_1 and length_2
    if length_1 is not None and length_1 < 0:
        raise ValueError("length_1 must be non-negative")
    if length_2 is not None and length_2 < 0:
        raise ValueError("length_2 must be non-negative")

    # Validate output_dir
    if output_dir is not None:
        if not output_dir.exists():
            output_dir.mkdir(parents=True, exist_ok=True)
        elif not output_dir.is_dir():
            raise NotADirectoryError(f"output_dir is not a directory: {output_dir}")

    # Validate path_to_cutadapt
    cutadapt_cmd = "cutadapt"
    if path_to_cutadapt is not None:
        if not path_to_cutadapt.exists():
            raise FileNotFoundError(f"cutadapt executable not found at: {path_to_cutadapt}")
        cutadapt_cmd = str(path_to_cutadapt)

    # Build command line
    cmd = ["trim_galore"]

    # Paired-end
    if paired:
        cmd.append("--paired")

    # Quality cutoff
    cmd.extend(["-q", str(quality)])

    # Phred encoding
    if phred64:
        cmd.append("--phred64")
    else:
        cmd.append("--phred33")

    # FastQC
    if fastqc:
        cmd.append("--fastqc")
        if fastqc_args:
            cmd.extend(["--fastqc_args", fastqc_args])

    # Adapter options
    if adapter is not None:
        cmd.extend(["-a", adapter])
    else:
        if illumina:
            cmd.append("--illumina")
        elif nextera:
            cmd.append("--nextera")
        elif small_rna:
            cmd.append("--small_rna")

    if adapter2 is not None:
        cmd.extend(["-a2", adapter2])

    # Length options
    if length > 0:
        cmd.extend(["--length", str(length)])

    if max_length is not None:
        if max_length < length:
            raise ValueError("max_length must be >= length")
        cmd.extend(["--max_length", str(max_length)])

    # Stringency and error rate
    cmd.extend(["--stringency", str(stringency)])
    cmd.extend(["-e", str(error_rate)])

    # Gzip options
    if dont_gzip:
        cmd.append("--dont_gzip")
    elif gzip_output:
        cmd.append("--gzip")

    # Retain unpaired reads
    if retain_unpaired:
        cmd.append("--retain_unpaired")

    # max_n and trim_n
    if max_n is not None:
        if max_n < 0:
            raise ValueError("max_n must be non-negative")
        cmd.extend(["--max_n", str(max_n)])
    if trim_n:
        cmd.append("--trim-n")

    # Output directory
    if output_dir is not None:
        cmd.extend(["-o", str(output_dir)])

    # No report file
    if no_report_file:
        cmd.append("--no_report_file")

    # Suppress warnings
    if suppress_warn:
        cmd.append("--suppress_warn")

    # Clipping options
    if clip_R1 > 0:
        cmd.extend(["--clip_R1", str(clip_R1)])
    if clip_R2 > 0:
        cmd.extend(["--clip_R2", str(clip_R2)])
    if three_prime_clip_R1 > 0:
        cmd.extend(["--three_prime_clip_R1", str(three_prime_clip_R1)])
    if three_prime_clip_R2 > 0:
        cmd.extend(["--three_prime_clip_R2", str(three_prime_clip_R2)])

    # Two-colour / NextSeq trimming
    if two_colour:
        cmd.append("--2colour")
    if nextseq_trim is not None:
        if nextseq_trim < 0:
            raise ValueError("nextseq_trim must be non-negative")
        cmd.append(f"--nextseq-trim={nextseq_trim}")

    # Basename
    if basename is not None:
        cmd.extend(["--basename", basename])

    # Cores
    cmd.extend(["-j", str(cores)])

    # Demux barcode file
    if demux is not None:
        if not demux.exists():
            raise FileNotFoundError(f"Demux barcode file not found: {demux}")
        cmd.extend(["--demux", str(demux)])

    # Hard trim
    if hardtrim5 is not None:
        cmd.extend(["--hardtrim5", str(hardtrim5)])
    if hardtrim3 is not None:
        cmd.extend(["--hardtrim3", str(hardtrim3)])

    # Clock mode
    if clock:
        cmd.append("--clock")

    # Dual index mode
    if dual_index:
        cmd.append("--dual_index")

    # RRBS mode
    if rrbs:
        cmd.append("--rrbs")

    # Non-directional RRBS
    if non_directional:
        cmd.append("--non_directional")

    # Keep temporary files
    if keep:
        cmd.append("--keep")

    # Trim1
    if trim1:
        cmd.append("-t")

    # Length_1 and length_2
    if length_1 is not None:
        cmd.extend(["-r1", str(length_1)])
    if length_2 is not None:
        cmd.extend(["-r2", str(length_2)])

    # Append input files (convert Path to str)
    cmd.extend([str(f) for f in input_files])

    # Use specified cutadapt path if given
    if path_to_cutadapt is not None:
        # trim_galore uses cutadapt internally, set environment variable
        import os
        env = dict(**os.environ)
        env["PATH"] = f"{path_to_cutadapt.parent}:{env.get('PATH','')}"
    else:
        env = None

    try:
        completed = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"trim_galore failed with exit code {e.returncode}",
        }

    # Collect output files if output_dir is specified, else guess from input files and basename
    output_files = []
    if output_dir is not None:
        # Collect all files in output_dir that start with basename or input file stem
        base_names = []
        if basename is not None:
            base_names.append(basename)
        else:
            for f in input_files:
                base_names.append(f.stem)
        for base in base_names:
            # Possible suffixes from trim_galore output
            suffixes = [
                "_trimmed.fq.gz",
                "_trimmed.fq",
                "_val_1.fq.gz",
                "_val_2.fq.gz",
                "_val_1.fq",
                "_val_2.fq",
                ".fq.gz",
                ".fq",
            ]
            for suf in suffixes:
                candidate = output_dir / (base + suf)
                if candidate.exists():
                    output_files.append(str(candidate))
    else:
        # If no output_dir, try to guess output files in same dir as inputs
        base_names = []
        if basename is not None:
            base_names.append(basename)
        else:
            for f in input_files:
                base_names.append(f.stem)
        for base in base_names:
            suffixes = [
                "_trimmed.fq.gz",
                "_trimmed.fq",
                "_val_1.fq.gz",
                "_val_2.fq.gz",
                "_val_1.fq",
                "_val_2.fq",
                ".fq.gz",
                ".fq",
            ]
            for suf in suffixes:
                candidate = Path(input_files[0].parent) / (base + suf)
                if candidate.exists():
                    output_files.append(str(candidate))

    return {
        "command_executed": " ".join(cmd),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "output_files": output_files,
    }


if __name__ == '__main__':
    mcp.run()