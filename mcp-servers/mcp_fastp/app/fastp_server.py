from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()

@mcp.tool()
def fastp(
    in1: Path,
    out1: Optional[Path] = None,
    in2: Optional[Path] = None,
    out2: Optional[Path] = None,
    unpaired1: Optional[Path] = None,
    unpaired2: Optional[Path] = None,
    failed_out: Optional[Path] = None,
    merge: bool = False,
    merged_out: Optional[str] = None,
    include_unmerged: bool = False,
    phred64: bool = False,
    compression: int = 4,
    stdin: bool = False,
    stdout: bool = False,
    interleaved_in: bool = False,
    reads_to_process: int = 0,
    dont_overwrite: bool = False,
    fix_mgi_id: bool = False,
    disable_adapter_trimming: bool = False,
    adapter_sequence: Optional[str] = None,
    adapter_sequence_r2: Optional[str] = None,
    adapter_fasta: Optional[Path] = None,
    detect_adapter_for_pe: bool = False,
    trim_front1: int = 0,
    trim_tail1: int = 0,
    max_len1: int = 0,
    trim_front2: Optional[int] = None,
    trim_tail2: Optional[int] = None,
    max_len2: Optional[int] = None,
    dedup: bool = False,
    dup_calc_accuracy: int = 0,
    dont_eval_duplication: bool = False,
    trim_poly_g: bool = False,
    poly_g_min_len: int = 10,
    disable_trim_poly_g: bool = False,
    trim_poly_x: bool = False,
    poly_x_min_len: int = 10,
    cut_front: bool = False,
    cut_tail: bool = False,
    cut_right: bool = False,
    cut_window_size: int = 4,
    cut_mean_quality: int = 20,
    cut_front_window_size: Optional[int] = None,
    cut_front_mean_quality: Optional[int] = None,
    cut_tail_window_size: Optional[int] = None,
    cut_tail_mean_quality: Optional[int] = None,
    cut_right_window_size: Optional[int] = None,
    cut_right_mean_quality: Optional[int] = None,
    disable_quality_filtering: bool = False,
    qualified_quality_phred: int = 15,
    unqualified_percent_limit: int = 40,
    n_base_limit: int = 5,
    average_qual: int = 0,
    disable_length_filtering: bool = False,
    length_required: int = 15,
    length_limit: int = 0,
    low_complexity_filter: bool = False,
    complexity_threshold: int = 30,
    filter_by_index1: Optional[Path] = None,
    filter_by_index2: Optional[Path] = None,
    filter_by_index_threshold: int = 0,
    correction: bool = False,
    overlap_len_require: int = 30,
    overlap_diff_limit: int = 5,
    overlap_diff_percent_limit: int = 20,
    umi: bool = False,
    umi_loc: Optional[str] = None,
    umi_len: int = 0,
    umi_prefix: Optional[str] = None,
    umi_skip: int = 0,
    overrepresentation_analysis: bool = False,
    overrepresentation_sampling: int = 20,
    json: str = "fastp.json",
    html: str = "fastp.html",
    report_title: str = "fastp report",
    thread: int = 2,
    split: int = 0,
    split_by_lines: int = 0,
    split_prefix_digits: int = 4,
    verbose: bool = False,
) -> dict:
    """
    fastp: ultra-fast all-in-one FASTQ preprocessor for quality control and filtering.
    Supports single-end and paired-end data with extensive options for trimming,
    filtering, deduplication, UMI processing, and reporting.
    """
    # Validate input files
    if not stdin:
        if not in1.exists():
            raise FileNotFoundError(f"Input file read1 does not exist: {in1}")
        if in2 is not None and not in2.exists():
            raise FileNotFoundError(f"Input file read2 does not exist: {in2}")
    # Validate adapter fasta file if provided
    if adapter_fasta is not None and not adapter_fasta.exists():
        raise FileNotFoundError(f"Adapter fasta file does not exist: {adapter_fasta}")
    # Validate compression level
    if not (1 <= compression <= 9):
        raise ValueError("compression must be between 1 and 9")
    # Validate dup_calc_accuracy
    if not (0 <= dup_calc_accuracy <= 6):
        raise ValueError("dup_calc_accuracy must be between 0 and 6")
    # Validate quality cut parameters ranges
    if not (1 <= cut_window_size <= 1000):
        raise ValueError("cut_window_size must be between 1 and 1000")
    if not (1 <= cut_mean_quality <= 36):
        raise ValueError("cut_mean_quality must be between 1 and 36")
    # Validate unqualified_percent_limit
    if not (0 <= unqualified_percent_limit <= 100):
        raise ValueError("unqualified_percent_limit must be between 0 and 100")
    # Validate complexity_threshold
    if not (0 <= complexity_threshold <= 100):
        raise ValueError("complexity_threshold must be between 0 and 100")
    # Validate filter_by_index_threshold
    if filter_by_index_threshold < 0:
        raise ValueError("filter_by_index_threshold must be >= 0")
    # Validate thread count
    if thread < 1:
        raise ValueError("thread must be >= 1")
    # Validate split options
    if split != 0 and split_by_lines != 0:
        raise ValueError("Cannot enable both split and split_by_lines simultaneously")
    if split != 0 and not (2 <= split <= 999):
        raise ValueError("split must be between 2 and 999")
    if split_prefix_digits < 0 or split_prefix_digits > 10:
        raise ValueError("split_prefix_digits must be between 0 and 10")
    # Set defaults for trim_front2, trim_tail2, max_len2 if not specified
    if trim_front2 is None:
        trim_front2 = trim_front1
    if trim_tail2 is None:
        trim_tail2 = trim_tail1
    if max_len2 is None:
        max_len2 = max_len1
    # Check output files existence if dont_overwrite is True
    output_files_to_check = []
    if out1 is not None:
        output_files_to_check.append(out1)
    if out2 is not None:
        output_files_to_check.append(out2)
    if unpaired1 is not None:
        output_files_to_check.append(unpaired1)
    if unpaired2 is not None:
        output_files_to_check.append(unpaired2)
    if failed_out is not None:
        output_files_to_check.append(failed_out)
    if merged_out is not None and merged_out != "--stdout":
        output_files_to_check.append(Path(merged_out))
    output_files_to_check.append(Path(json))
    output_files_to_check.append(Path(html))
    if dont_overwrite:
        for f in output_files_to_check:
            if f.exists():
                raise FileExistsError(f"Output file exists and --dont_overwrite enabled: {f}")
    # Build command line
    cmd = ["fastp"]
    # Input/output
    if stdin:
        cmd.append("--stdin")
    else:
        cmd.extend(["-i", str(in1)])
        if out1 is not None:
            cmd.extend(["-o", str(out1)])
        if in2 is not None:
            cmd.extend(["-I", str(in2)])
            if out2 is not None:
                cmd.extend(["-O", str(out2)])
    if unpaired1 is not None:
        cmd.extend(["--unpaired1", str(unpaired1)])
    if unpaired2 is not None:
        cmd.extend(["--unpaired2", str(unpaired2)])
    if failed_out is not None:
        cmd.extend(["--failed_out", str(failed_out)])
    if merge:
        cmd.append("-m")
        if merged_out is not None:
            if merged_out == "--stdout":
                cmd.append("--merged_out")
                cmd.append("--stdout")
            else:
                cmd.extend(["--merged_out", merged_out])
        else:
            # merged_out must be specified or stdout enabled in merge mode
            raise ValueError("In merge mode, --merged_out or --stdout must be specified")
    if include_unmerged:
        cmd.append("--include_unmerged")
    if phred64:
        cmd.append("-6")
    cmd.extend(["-z", str(compression)])
    if stdout:
        cmd.append("--stdout")
    if interleaved_in:
        cmd.append("--interleaved_in")
    if reads_to_process > 0:
        cmd.extend(["--reads_to_process", str(reads_to_process)])
    if dont_overwrite:
        cmd.append("--dont_overwrite")
    if fix_mgi_id:
        cmd.append("--fix_mgi_id")
    if disable_adapter_trimming:
        cmd.append("-A")
    if adapter_sequence is not None:
        cmd.extend(["-a", adapter_sequence])
    if adapter_sequence_r2 is not None:
        cmd.extend(["--adapter_sequence_r2", adapter_sequence_r2])
    if adapter_fasta is not None:
        cmd.extend(["--adapter_fasta", str(adapter_fasta)])
    if detect_adapter_for_pe:
        cmd.append("--detect_adapter_for_pe")
    # Global trimming
    cmd.extend(["-f", str(trim_front1)])
    cmd.extend(["-t", str(trim_tail1)])
    cmd.extend(["-b", str(max_len1)])
    cmd.extend(["-F", str(trim_front2)])
    cmd.extend(["-T", str(trim_tail2)])
    cmd.extend(["-B", str(max_len2)])
    # Deduplication
    if dedup:
        cmd.append("-D")
    cmd.extend(["--dup_calc_accuracy", str(dup_calc_accuracy)])
    if dont_eval_duplication:
        cmd.append("--dont_eval_duplication")
    # PolyG trimming
    if trim_poly_g:
        cmd.append("-g")
    if disable_trim_poly_g:
        cmd.append("-G")
    cmd.extend(["--poly_g_min_len", str(poly_g_min_len)])
    # PolyX trimming
    if trim_poly_x:
        cmd.append("-x")
    cmd.extend(["--poly_x_min_len", str(poly_x_min_len)])
    # Per read cutting by quality
    if cut_front:
        cmd.append("-5")
    if cut_tail:
        cmd.append("-3")
    if cut_right:
        cmd.append("-r")
    cmd.extend(["-W", str(cut_window_size)])
    cmd.extend(["-M", str(cut_mean_quality)])
    if cut_front_window_size is not None:
        cmd.extend(["--cut_front_window_size", str(cut_front_window_size)])
    if cut_front_mean_quality is not None:
        cmd.extend(["--cut_front_mean_quality", str(cut_front_mean_quality)])
    if cut_tail_window_size is not None:
        cmd.extend(["--cut_tail_window_size", str(cut_tail_window_size)])
    if cut_tail_mean_quality is not None:
        cmd.extend(["--cut_tail_mean_quality", str(cut_tail_mean_quality)])
    if cut_right_window_size is not None:
        cmd.extend(["--cut_right_window_size", str(cut_right_window_size)])
    if cut_right_mean_quality is not None:
        cmd.extend(["--cut_right_mean_quality", str(cut_right_mean_quality)])
    # Quality filtering
    if disable_quality_filtering:
        cmd.append("-Q")
    cmd.extend(["-q", str(qualified_quality_phred)])
    cmd.extend(["-u", str(unqualified_percent_limit)])
    cmd.extend(["-n", str(n_base_limit)])
    cmd.extend(["-e", str(average_qual)])
    # Length filtering
    if disable_length_filtering:
        cmd.append("-L")
    cmd.extend(["-l", str(length_required)])
    cmd.extend(["--length_limit", str(length_limit)])
    # Low complexity filtering
    if low_complexity_filter:
        cmd.append("-y")
    cmd.extend(["-Y", str(complexity_threshold)])
    # Filter by index
    if filter_by_index1 is not None:
        if not filter_by_index1.exists():
            raise FileNotFoundError(f"filter_by_index1 file does not exist: {filter_by_index1}")
        cmd.extend(["--filter_by_index1", str(filter_by_index1)])
    if filter_by_index2 is not None:
        if not filter_by_index2.exists():
            raise FileNotFoundError(f"filter_by_index2 file does not exist: {filter_by_index2}")
        cmd.extend(["--filter_by_index2", str(filter_by_index2)])
    cmd.extend(["--filter_by_index_threshold", str(filter_by_index_threshold)])
    # Base correction by overlap analysis
    if correction:
        cmd.append("-c")
    cmd.extend(["--overlap_len_require", str(overlap_len_require)])
    cmd.extend(["--overlap_diff_limit", str(overlap_diff_limit)])
    cmd.extend(["--overlap_diff_percent_limit", str(overlap_diff_percent_limit)])
    # UMI processing
    if umi:
        cmd.append("-U")
        if umi_loc is not None:
            if umi_loc not in ("index1", "index2", "read1", "read2", "per_index", "per_read"):
                raise ValueError("umi_loc must be one of: index1, index2, read1, read2, per_index, per_read")
            cmd.extend(["--umi_loc", umi_loc])
        cmd.extend(["--umi_len", str(umi_len)])
        if umi_prefix is not None:
            cmd.extend(["--umi_prefix", umi_prefix])
        cmd.extend(["--umi_skip", str(umi_skip)])
    # Overrepresented sequence analysis
    if overrepresentation_analysis:
        cmd.append("-p")
    cmd.extend(["-P", str(overrepresentation_sampling)])
    # Reporting options
    cmd.extend(["-j", json])
    cmd.extend(["-h", html])
    cmd.extend(["-R", report_title])
    # Threading
    cmd.extend(["-w", str(thread)])
    # Output splitting
    if split != 0:
        cmd.extend(["-s", str(split)])
    if split_by_lines != 0:
        cmd.extend(["-S", str(split_by_lines)])
    cmd.extend(["-d", str(split_prefix_digits)])
    # Verbose
    if verbose:
        cmd.append("-V")

    # Run command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"fastp failed with return code {e.returncode}",
            "output_files": []
        }

    # Collect output files
    output_files = []
    if out1 is not None:
        output_files.append(str(out1))
    if out2 is not None:
        output_files.append(str(out2))
    if unpaired1 is not None:
        output_files.append(str(unpaired1))
    if unpaired2 is not None:
        output_files.append(str(unpaired2))
    if failed_out is not None:
        output_files.append(str(failed_out))
    if merged_out is not None and merged_out != "--stdout":
        output_files.append(str(merged_out))
    output_files.append(json)
    output_files.append(html)

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()