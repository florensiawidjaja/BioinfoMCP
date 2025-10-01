from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


def validate_path_exists(path: Path, param_name: str):
    if not path.exists():
        raise FileNotFoundError(f"Input file for parameter '{param_name}' does not exist: {path}")


def validate_positive_int(value: int, param_name: str):
    if value < 1:
        raise ValueError(f"Parameter '{param_name}' must be a positive integer, got {value}")


def validate_quality_encoding(encoding: Optional[str]) -> Optional[str]:
    if encoding is None:
        return None
    if encoding not in ("phred33", "phred64"):
        raise ValueError("Quality encoding must be either 'phred33' or 'phred64'")
    return encoding


def build_steps_command(steps: List[str]) -> List[str]:
    # Each step is a string like "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10"
    # We just pass them as is to the command line
    return steps


@mcp.tool()
def trimmomatic_se(
    input_fastq: str,
    output_fastq: str,
    steps: List[str],
    threads: int = 1,
    phred33: bool = False,
    phred64: bool = False,
    trimlog: Optional[str] = None,
    trimmomatic_jar: str = "trimmomatic.jar",
) -> dict:
    """
    Run Trimmomatic in Single-End mode to trim and crop Illumina FASTQ data.

    Parameters:
    - input_fastq: Path to input FASTQ file (can be gzipped or bz2).
    - output_fastq: Path to output FASTQ file (can be gzipped or bz2).
    - steps: List of trimming steps to apply, e.g. ["ILLUMINACLIP:TruSeq3-SE.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"].
    - threads: Number of threads to use (default 1).
    - phred33: Use Phred+33 quality encoding (default False).
    - phred64: Use Phred+64 quality encoding (default False).
    - trimlog: Optional path to trim log file.
    - trimmomatic_jar: Path to the Trimmomatic jar file.

    Returns:
    A dictionary with command executed, stdout, stderr, and output files.
    """
    # Validate input/output paths
    input_path = Path(input_fastq)
    output_path = Path(output_fastq)
    validate_path_exists(input_path, "input_fastq")
    if output_path.exists():
        # Overwrite allowed but warn could be added if needed
        pass

    # Validate threads
    validate_positive_int(threads, "threads")

    # Validate quality encoding flags
    if phred33 and phred64:
        raise ValueError("Only one of phred33 or phred64 can be True")
    quality_flag = None
    if phred33:
        quality_flag = "-phred33"
    elif phred64:
        quality_flag = "-phred64"

    # Validate steps list
    if not steps or not isinstance(steps, list):
        raise ValueError("Parameter 'steps' must be a non-empty list of trimming steps")

    # Build command
    cmd = [
        "java",
        "-jar",
        trimmomatic_jar,
        "SE",
        "-threads",
        str(threads),
    ]
    if quality_flag:
        cmd.append(quality_flag)
    if trimlog:
        trimlog_path = Path(trimlog)
        cmd.extend(["-trimlog", str(trimlog_path)])
    cmd.extend([str(input_path), str(output_path)])
    cmd.extend(build_steps_command(steps))

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"Trimmomatic SE failed with exit code {e.returncode}",
            "output_files": [],
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(output_path)] + ([str(trimlog_path)] if trimlog else []),
    }


@mcp.tool()
def trimmomatic_pe(
    input_forward: Optional[str] = None,
    input_reverse: Optional[str] = None,
    output_paired_forward: Optional[str] = None,
    output_unpaired_forward: Optional[str] = None,
    output_paired_reverse: Optional[str] = None,
    output_unpaired_reverse: Optional[str] = None,
    basein: Optional[str] = None,
    baseout: Optional[str] = None,
    steps: List[str] = [],
    threads: int = 1,
    phred33: bool = False,
    phred64: bool = False,
    trimlog: Optional[str] = None,
    trimmomatic_jar: str = "trimmomatic.jar",
) -> dict:
    """
    Run Trimmomatic in Paired-End mode to trim and crop Illumina paired-end FASTQ data.

    Parameters:
    - input_forward: Path to forward input FASTQ file (gzipped or bz2 supported).
    - input_reverse: Path to reverse input FASTQ file.
    - output_paired_forward: Path to output forward paired reads file.
    - output_unpaired_forward: Path to output forward unpaired reads file.
    - output_paired_reverse: Path to output reverse paired reads file.
    - output_unpaired_reverse: Path to output reverse unpaired reads file.
    - basein: Optional base name for input files; if provided, input_forward and input_reverse are ignored.
    - baseout: Optional base name for output files; if provided, output files are derived automatically.
    - steps: List of trimming steps to apply, e.g. ["ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"].
    - threads: Number of threads to use (default 1).
    - phred33: Use Phred+33 quality encoding (default False).
    - phred64: Use Phred+64 quality encoding (default False).
    - trimlog: Optional path to trim log file.
    - trimmomatic_jar: Path to the Trimmomatic jar file.

    Returns:
    A dictionary with command executed, stdout, stderr, and output files.
    """
    # Validate threads
    validate_positive_int(threads, "threads")

    # Validate quality encoding flags
    if phred33 and phred64:
        raise ValueError("Only one of phred33 or phred64 can be True")
    quality_flag = None
    if phred33:
        quality_flag = "-phred33"
    elif phred64:
        quality_flag = "-phred64"

    # Validate steps list
    if not steps or not isinstance(steps, list):
        raise ValueError("Parameter 'steps' must be a non-empty list of trimming steps")

    # Determine input files
    if basein:
        basein_path = Path(basein)
        # We do not check existence of basein itself, but of derived files
        # Try to guess reverse file from forward file name patterns:
        # Common patterns: _R1_ -> _R2_, .1. -> .2., .f. -> .r.
        # We will try to find the forward file and reverse file by these patterns
        # We check existence of both files
        # Try common patterns in order:
        candidates = []
        basein_str = str(basein_path)
        # Pattern 1: _R1_ -> _R2_
        if "_R1_" in basein_str:
            forward_candidate = basein_path
            reverse_candidate = Path(basein_str.replace("_R1_", "_R2_"))
            candidates.append((forward_candidate, reverse_candidate))
        # Pattern 2: .1. -> .2.
        if ".1." in basein_str:
            forward_candidate = basein_path
            reverse_candidate = Path(basein_str.replace(".1.", ".2."))
            candidates.append((forward_candidate, reverse_candidate))
        # Pattern 3: .f. -> .r.
        if ".f." in basein_str:
            forward_candidate = basein_path
            reverse_candidate = Path(basein_str.replace(".f.", ".r."))
            candidates.append((forward_candidate, reverse_candidate))
        # If no pattern matched, fallback to basein + "_1.fastq" and basein + "_2.fastq"
        if not candidates:
            forward_candidate = basein_path.with_name(basein_path.stem + "_1" + basein_path.suffix)
            reverse_candidate = basein_path.with_name(basein_path.stem + "_2" + basein_path.suffix)
            candidates.append((forward_candidate, reverse_candidate))

        # Find first pair where both files exist
        found_pair = None
        for fwd, rev in candidates:
            if fwd.exists() and rev.exists():
                found_pair = (fwd, rev)
                break
        if not found_pair:
            raise FileNotFoundError(f"Could not find paired input files derived from basein '{basein}'")
        input_forward_path, input_reverse_path = found_pair
    else:
        if not input_forward or not input_reverse:
            raise ValueError("Either basein or both input_forward and input_reverse must be specified")
        input_forward_path = Path(input_forward)
        input_reverse_path = Path(input_reverse)
        validate_path_exists(input_forward_path, "input_forward")
        validate_path_exists(input_reverse_path, "input_reverse")

    # Determine output files
    if baseout:
        baseout_path = Path(baseout)
        # Derive 4 output files:
        # mySampleFiltered_1P.fq.gz - paired forward
        # mySampleFiltered_1U.fq.gz - unpaired forward
        # mySampleFiltered_2P.fq.gz - paired reverse
        # mySampleFiltered_2U.fq.gz - unpaired reverse
        suffix = baseout_path.suffix  # includes dot, e.g. ".fq.gz"
        stem = baseout_path.stem
        # If suffix is double extension like .fq.gz, stem will be "mySampleFiltered.fq"
        # So we handle double extensions by checking suffixes
        # We handle .fq.gz or .fastq.gz or .fq.bz2 etc.
        # To do this robustly, we check suffixes until no more known compression suffixes
        # For simplicity, we just use the full baseout string and append suffixes
        # So output files are baseout + "_1P" + suffix etc.
        # But to avoid double suffix, we remove compression suffixes from suffix and add them after
        # Known compression suffixes:
        compression_suffixes = [".gz", ".bz2"]
        comp_suffix = ""
        base_suffix = suffix
        for csuf in compression_suffixes:
            if suffix.endswith(csuf):
                comp_suffix = csuf
                base_suffix = suffix[:-len(csuf)]
                break
        # Compose output files
        paired_forward = str(baseout_path.with_name(baseout_path.stem + "_1P" + base_suffix + comp_suffix))
        unpaired_forward = str(baseout_path.with_name(baseout_path.stem + "_1U" + base_suffix + comp_suffix))
        paired_reverse = str(baseout_path.with_name(baseout_path.stem + "_2P" + base_suffix + comp_suffix))
        unpaired_reverse = str(baseout_path.with_name(baseout_path.stem + "_2U" + base_suffix + comp_suffix))
    else:
        # All 4 output files must be specified explicitly
        if (
            not output_paired_forward
            or not output_unpaired_forward
            or not output_paired_reverse
            or not output_unpaired_reverse
        ):
            raise ValueError(
                "Either baseout or all four output files (output_paired_forward, output_unpaired_forward, output_paired_reverse, output_unpaired_reverse) must be specified"
            )
        paired_forward = output_paired_forward
        unpaired_forward = output_unpaired_forward
        paired_reverse = output_paired_reverse
        unpaired_reverse = output_unpaired_reverse

    # Build command
    cmd = [
        "java",
        "-jar",
        trimmomatic_jar,
        "PE",
        "-threads",
        str(threads),
    ]
    if quality_flag:
        cmd.append(quality_flag)
    if trimlog:
        trimlog_path = Path(trimlog)
        cmd.extend(["-trimlog", str(trimlog_path)])

    # Input files or basein
    if basein:
        cmd.extend(["-basein", str(basein_path)])
    else:
        cmd.extend([str(input_forward_path), str(input_reverse_path)])

    # Output files or baseout
    if baseout:
        cmd.extend(["-baseout", str(baseout_path)])
    else:
        cmd.extend([paired_forward, unpaired_forward, paired_reverse, unpaired_reverse])

    # Add steps
    cmd.extend(build_steps_command(steps))

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"Trimmomatic PE failed with exit code {e.returncode}",
            "output_files": [],
        }

    output_files = []
    if baseout:
        output_files = [paired_forward, unpaired_forward, paired_reverse, unpaired_reverse]
    else:
        output_files = [paired_forward, unpaired_forward, paired_reverse, unpaired_reverse]
    if trimlog:
        output_files.append(str(trimlog_path))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files,
    }


if __name__ == "__main__":
    mcp.run()