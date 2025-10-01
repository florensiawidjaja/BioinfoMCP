from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional
import subprocess

mcp = FastMCP()

@mcp.tool()
def gatk_ApplyBQSR(
    bqsr_recal_file: Path,
    input: List[Path],
    output: Path,
    reference: Optional[Path] = None,
    intervals: Optional[List[str]] = None,
    exclude_intervals: Optional[List[str]] = None,
    interval_merging_rule: str = "ALL",
    interval_set_rule: str = "UNION",
    interval_padding: int = 0,
    interval_exclusion_padding: int = 0,
    preserve_qscores_less_than: int = 6,
    quantize_quals: int = 0,
    static_quantized_quals: Optional[List[int]] = None,
    round_down_quantized: bool = False,
    emit_original_quals: bool = False,
    use_original_qualities: bool = False,
    create_output_bam_index: bool = True,
    create_output_bam_md5: bool = False,
    disable_bam_index_caching: bool = False,
    disable_sequence_dictionary_validation: bool = False,
    disable_tool_default_read_filters: bool = False,
    disable_read_filter: Optional[List[str]] = None,
    read_filter: Optional[List[str]] = None,
    inverted_read_filter: Optional[List[str]] = None,
    read_index: Optional[List[Path]] = None,
    read_validation_stringency: str = "SILENT",
    add_output_sam_program_record: bool = True,
    add_output_vcf_command_line: bool = True,
    gatk_config_file: Optional[Path] = None,
    cloud_index_prefetch_buffer: int = -1,
    cloud_prefetch_buffer: int = 40,
    gcs_max_retries: int = 20,
    gcs_project_for_requester_pays: str = "",
    lenient: bool = False,
    max_variants_per_shard: int = 0,
    quiet: bool = False,
    seconds_between_progress_updates: float = 10.0,
    sequence_dictionary: Optional[Path] = None,
    sites_only_vcf_output: bool = False,
    tmp_dir: Optional[Path] = None,
    use_jdk_deflater: bool = False,
    use_jdk_inflater: bool = False,
    verbosity: str = "INFO",
    version: bool = False,
    help: bool = False,
    arguments_file: Optional[List[Path]] = None,
    show_hidden: bool = False,
):
    """
    Apply base quality score recalibration (BQSR) to BAM/CRAM files using a recalibration table.
    This tool recalibrates base qualities of input reads based on the recalibration table produced by BaseRecalibrator.
    
    Parameters:
    - bqsr_recal_file: Input recalibration table file produced by BaseRecalibrator (required).
    - input: List of BAM/SAM/CRAM input files containing reads (required).
    - output: Output BAM/CRAM file with recalibrated reads (required).
    - reference: Reference sequence FASTA file.
    - intervals: Genomic intervals over which to operate.
    - exclude_intervals: Genomic intervals to exclude from processing.
    - interval_merging_rule: Rule for merging abutting intervals ("ALL" or "OVERLAPPING_ONLY").
    - interval_set_rule: Set merging approach for combining interval inputs ("UNION" or "INTERSECTION").
    - interval_padding: Padding (bp) to add to each included interval.
    - interval_exclusion_padding: Padding (bp) to add to each excluded interval.
    - preserve_qscores_less_than: Do not recalibrate bases with quality scores less than this threshold.
    - quantize_quals: Quantize quality scores to given number of levels (0 disables).
    - static_quantized_quals: Use static quantized quality scores to given levels (cannot be used with quantize_quals).
    - round_down_quantized: Round quals down to nearest quantized qual (only with static_quantized_quals).
    - emit_original_quals: Emit original base qualities under the OQ tag.
    - use_original_qualities: Use base quality scores from the OQ tag if present.
    - create_output_bam_index: Create BAM/CRAM index for output.
    - create_output_bam_md5: Create MD5 digest for output BAM/CRAM.
    - disable_bam_index_caching: Disable caching of BAM indexes.
    - disable_sequence_dictionary_validation: Disable sequence dictionary compatibility check.
    - disable_tool_default_read_filters: Disable all tool default read filters.
    - disable_read_filter: List of read filters to disable.
    - read_filter: List of read filters to apply.
    - inverted_read_filter: List of inverted read filters to apply.
    - read_index: List of indices for read inputs.
    - read_validation_stringency: Validation stringency ("STRICT", "LENIENT", "SILENT").
    - add_output_sam_program_record: Add PG tag to output SAM/BAM/CRAM files.
    - add_output_vcf_command_line: Add command line header line to output VCF files.
    - gatk_config_file: Configuration file for GATK.
    - cloud_index_prefetch_buffer: Size of cloud-only prefetch buffer (MB).
    - cloud_prefetch_buffer: Size of cloud-only prefetch buffer (MB).
    - gcs_max_retries: Max retries for GCS bucket errors.
    - gcs_project_for_requester_pays: Project to bill for requester pays buckets.
    - lenient: Lenient processing of VCF files.
    - max_variants_per_shard: Partition VCF output into shards with max records.
    - quiet: Suppress job-summary info on stderr.
    - seconds_between_progress_updates: Interval for output traversal statistics.
    - sequence_dictionary: Sequence dictionary file (.dict).
    - sites_only_vcf_output: Do not emit genotype fields in VCF output.
    - tmp_dir: Temporary directory to use.
    - use_jdk_deflater: Use JdkDeflater instead of IntelDeflater.
    - use_jdk_inflater: Use JdkInflater instead of IntelInflater.
    - verbosity: Logging verbosity ("ERROR", "WARNING", "INFO", "DEBUG").
    - version: Display version number.
    - help: Display help message.
    - arguments_file: List of argument files to add to command line.
    - show_hidden: Display hidden arguments.
    """
    # Validate required files
    if not bqsr_recal_file.exists():
        raise FileNotFoundError(f"Recalibration file not found: {bqsr_recal_file}")
    if not output.parent.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output.parent}")
    for f in input:
        if not f.exists():
            raise FileNotFoundError(f"Input file not found: {f}")
    if reference is not None and not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if gatk_config_file is not None and not gatk_config_file.exists():
        raise FileNotFoundError(f"GATK config file not found: {gatk_config_file}")
    if sequence_dictionary is not None and not sequence_dictionary.exists():
        raise FileNotFoundError(f"Sequence dictionary file not found: {sequence_dictionary}")
    if read_index is not None:
        for idx in read_index:
            if not idx.exists():
                raise FileNotFoundError(f"Read index file not found: {idx}")
    if arguments_file is not None:
        for argf in arguments_file:
            if not argf.exists():
                raise FileNotFoundError(f"Arguments file not found: {argf}")
    if tmp_dir is not None and not tmp_dir.exists():
        raise FileNotFoundError(f"Temporary directory not found: {tmp_dir}")

    # Validate parameter constraints
    if interval_merging_rule not in ("ALL", "OVERLAPPING_ONLY"):
        raise ValueError("interval_merging_rule must be 'ALL' or 'OVERLAPPING_ONLY'")
    if interval_set_rule not in ("UNION", "INTERSECTION"):
        raise ValueError("interval_set_rule must be 'UNION' or 'INTERSECTION'")
    if read_validation_stringency not in ("STRICT", "LENIENT", "SILENT"):
        raise ValueError("read_validation_stringency must be one of 'STRICT', 'LENIENT', 'SILENT'")
    if verbosity not in ("ERROR", "WARNING", "INFO", "DEBUG"):
        raise ValueError("verbosity must be one of 'ERROR', 'WARNING', 'INFO', 'DEBUG'")
    if preserve_qscores_less_than < 0:
        raise ValueError("preserve_qscores_less_than must be >= 0")
    if quantize_quals < 0:
        raise ValueError("quantize_quals must be >= 0")
    if max_variants_per_shard < 0:
        raise ValueError("max_variants_per_shard must be >= 0")
    if seconds_between_progress_updates < 0:
        raise ValueError("seconds_between_progress_updates must be >= 0")
    if interval_padding < 0:
        raise ValueError("interval_padding must be >= 0")
    if interval_exclusion_padding < 0:
        raise ValueError("interval_exclusion_padding must be >= 0")
    if static_quantized_quals is not None and quantize_quals != 0:
        raise ValueError("Cannot use static_quantized_quals and quantize_quals simultaneously")
    if round_down_quantized and quantize_quals != 0:
        raise ValueError("round_down_quantized cannot be used with quantize_quals")

    # Build command line
    cmd = ["gatk", "ApplyBQSR"]

    # Required arguments
    cmd += ["--bqsr-recal-file", str(bqsr_recal_file)]
    for in_file in input:
        cmd += ["-I", str(in_file)]
    cmd += ["-O", str(output)]

    # Optional arguments
    if reference:
        cmd += ["-R", str(reference)]
    if intervals:
        for interval in intervals:
            cmd += ["-L", interval]
    if exclude_intervals:
        for excl in exclude_intervals:
            cmd += ["-XL", excl]
    cmd += ["-imr", interval_merging_rule]
    cmd += ["-isr", interval_set_rule]
    if interval_padding > 0:
        cmd += ["-ip", str(interval_padding)]
    if interval_exclusion_padding > 0:
        cmd += ["-ixp", str(interval_exclusion_padding)]
    cmd += ["--preserve-qscores-less-than", str(preserve_qscores_less_than)]
    if quantize_quals != 0:
        cmd += ["--quantize-quals", str(quantize_quals)]
    if static_quantized_quals:
        for val in static_quantized_quals:
            cmd += ["--static-quantized-quals", str(val)]
    if round_down_quantized:
        cmd += ["--round-down-quantized"]
    if emit_original_quals:
        cmd += ["--emit-original-quals"]
    if use_original_qualities:
        cmd += ["--use-original-qualities"]
    if create_output_bam_index:
        cmd += ["-OBI"]
    else:
        cmd += ["--no-OBI"]
    if create_output_bam_md5:
        cmd += ["-OBM"]
    if disable_bam_index_caching:
        cmd += ["-DBIC"]
    if disable_sequence_dictionary_validation:
        cmd += ["--disable-sequence-dictionary-validation"]
    if disable_tool_default_read_filters:
        cmd += ["--disable-tool-default-read-filters"]
    if disable_read_filter:
        for filt in disable_read_filter:
            cmd += ["-DF", filt]
    if read_filter:
        for filt in read_filter:
            cmd += ["-RF", filt]
    if inverted_read_filter:
        for filt in inverted_read_filter:
            cmd += ["-XRF", filt]
    if read_index:
        for idx in read_index:
            cmd += ["--read-index", str(idx)]
    cmd += ["-VS", read_validation_stringency]
    if add_output_sam_program_record:
        cmd += ["--add-output-sam-program-record"]
    else:
        cmd += ["--no-add-output-sam-program-record"]
    if add_output_vcf_command_line:
        cmd += ["--add-output-vcf-command-line"]
    else:
        cmd += ["--no-add-output-vcf-command-line"]
    if gatk_config_file:
        cmd += ["--gatk-config-file", str(gatk_config_file)]
    if cloud_index_prefetch_buffer != -1:
        cmd += ["-CIPB", str(cloud_index_prefetch_buffer)]
    if cloud_prefetch_buffer != 40:
        cmd += ["-CPB", str(cloud_prefetch_buffer)]
    if gcs_max_retries != 20:
        cmd += ["-gcs-retries", str(gcs_max_retries)]
    if gcs_project_for_requester_pays:
        cmd += ["--gcs-project-for-requester-pays", gcs_project_for_requester_pays]
    if lenient:
        cmd += ["-LE"]
    if max_variants_per_shard != 0:
        cmd += ["--max-variants-per-shard", str(max_variants_per_shard)]
    if quiet:
        cmd += ["--QUIET"]
    if seconds_between_progress_updates != 10.0:
        cmd += ["-seconds-between-progress-updates", str(seconds_between_progress_updates)]
    if sequence_dictionary:
        cmd += ["-sequence-dictionary", str(sequence_dictionary)]
    if sites_only_vcf_output:
        cmd += ["--sites-only-vcf-output"]
    if tmp_dir:
        cmd += ["--tmp-dir", str(tmp_dir)]
    if use_jdk_deflater:
        cmd += ["-jdk-deflater"]
    if use_jdk_inflater:
        cmd += ["-jdk-inflater"]
    if verbosity != "INFO":
        cmd += ["-verbosity", verbosity]
    if version:
        cmd += ["--version"]
    if help:
        cmd += ["--help"]
    if arguments_file:
        for argf in arguments_file:
            cmd += ["--arguments_file", str(argf)]
    if show_hidden:
        cmd += ["--showHidden"]

    # Run command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "error": f"Command failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(output)]
    }


if __name__ == '__main__':
    mcp.run()