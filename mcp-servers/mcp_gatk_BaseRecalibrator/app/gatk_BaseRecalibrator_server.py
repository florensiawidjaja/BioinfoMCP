from fastmcp import FastMCP
from typing import List, Optional
from pathlib import Path
import subprocess

mcp = FastMCP()

@mcp.tool()
def gatk_BaseRecalibrator(
    input_reads: List[Path],
    reference: Path,
    known_sites: List[Path],
    output: Path,
    arguments_file: Optional[List[Path]] = None,
    binary_tag_name: Optional[str] = None,
    bqsr_baq_gap_open_penalty: float = 40.0,
    cloud_index_prefetch_buffer: int = -1,
    cloud_prefetch_buffer: int = 40,
    default_base_qualities: int = -1,
    deletions_default_quality: int = 45,
    disable_bam_index_caching: bool = False,
    disable_sequence_dictionary_validation: bool = False,
    gcs_max_retries: int = 20,
    gcs_project_for_requester_pays: str = "",
    indels_context_size: int = 3,
    insertions_default_quality: int = 45,
    interval_merging_rule: str = "ALL",
    intervals: Optional[List[str]] = None,
    low_quality_tail: int = 2,
    maximum_cycle_value: int = 500,
    mismatches_context_size: int = 2,
    mismatches_default_quality: int = -1,
    preserve_qscores_less_than: int = 6,
    quantizing_levels: int = 16,
    sites_only_vcf_output: bool = False,
    use_original_qualities: bool = False,
    add_output_sam_program_record: bool = True,
    add_output_vcf_command_line: bool = True,
    create_output_bam_index: bool = True,
    create_output_bam_md5: bool = False,
    create_output_variant_index: bool = True,
    create_output_variant_md5: bool = False,
    disable_read_filter: Optional[List[str]] = None,
    disable_tool_default_read_filters: bool = False,
    exclude_intervals: Optional[List[str]] = None,
    gatk_config_file: Optional[str] = None,
    gcs_retries: int = 20,
    help_flag: bool = False,
    interval_exclusion_padding: int = 0,
    interval_padding: int = 0,
    interval_set_rule: str = "UNION",
    inverted_read_filter: Optional[List[str]] = None,
    lenient: bool = False,
    max_variants_per_shard: int = 0,
    quiet: bool = False,
    read_filter: Optional[List[str]] = None,
    read_index: Optional[List[Path]] = None,
    read_validation_stringency: str = "SILENT",
    seconds_between_progress_updates: float = 10.0,
    sequence_dictionary: Optional[Path] = None,
    show_hidden: bool = False,
    tmp_dir: Optional[Path] = None,
    use_jdk_deflater: bool = False,
    use_jdk_inflater: bool = False,
    verbosity: str = "INFO",
    version: bool = False,
) -> dict:
    """
    Generates a recalibration table for Base Quality Score Recalibration (BQSR).
    This tool analyzes input reads against known polymorphic sites to produce a recalibration table.
    
    Parameters:
    - input_reads: List of BAM/SAM/CRAM files containing reads to be recalibrated.
    - reference: Reference genome fasta file.
    - known_sites: List of known polymorphic sites files (VCF/BCF/BED) to exclude from analysis.
    - output: Output recalibration table file path.
    - arguments_file: Optional list of argument files to add to the command line.
    - binary_tag_name: Optional binary tag covariate name.
    - bqsr_baq_gap_open_penalty: BQSR BAQ gap open penalty (Phred scaled).
    - cloud_index_prefetch_buffer: Size of cloud-only prefetch buffer (MB), -1 to default.
    - cloud_prefetch_buffer: Size of cloud-only prefetch buffer (MB).
    - default_base_qualities: Default base quality to assign if missing (-1 disables).
    - deletions_default_quality: Default quality for base deletions covariate.
    - disable_bam_index_caching: Disable BAM index caching to reduce memory.
    - disable_sequence_dictionary_validation: Disable sequence dictionary compatibility check.
    - gcs_max_retries: Max retries for GCS bucket errors.
    - gcs_project_for_requester_pays: Project to bill for requester pays buckets.
    - indels_context_size: k-mer context size for insertions/deletions covariate (1-13).
    - insertions_default_quality: Default quality for base insertions covariate.
    - interval_merging_rule: Interval merging rule for abutting intervals ("ALL" or "OVERLAPPING_ONLY").
    - intervals: List of genomic intervals to operate on.
    - low_quality_tail: Minimum quality for bases in read tails to be considered.
    - maximum_cycle_value: Max cycle value permitted for Cycle covariate.
    - mismatches_context_size: k-mer context size for base mismatches covariate (1-13).
    - mismatches_default_quality: Default quality for base mismatches covariate (-1 disables).
    - preserve_qscores_less_than: Bases with quality less than this are not recalibrated.
    - quantizing_levels: Number of distinct quality scores in quantized output.
    - sites_only_vcf_output: If true, omit genotype fields in VCF output.
    - use_original_qualities: Use base qualities from OQ tag if present.
    - add_output_sam_program_record: Add PG tag to output SAM/BAM/CRAM files.
    - add_output_vcf_command_line: Add command line header to VCF files.
    - create_output_bam_index: Create BAM/CRAM index for coordinate sorted output.
    - create_output_bam_md5: Create MD5 digest for BAM/SAM/CRAM output.
    - create_output_variant_index: Create VCF index for coordinate sorted output.
    - create_output_variant_md5: Create MD5 digest for VCF output.
    - disable_read_filter: List of read filters to disable.
    - disable_tool_default_read_filters: Disable all tool default read filters.
    - exclude_intervals: List of genomic intervals to exclude from processing.
    - gatk_config_file: Configuration file for GATK.
    - gcs_retries: Alias for gcs_max_retries.
    - help_flag: Display help message.
    - interval_exclusion_padding: Padding (bp) to add to excluded intervals.
    - interval_padding: Padding (bp) to add to included intervals.
    - interval_set_rule: Set merging approach for intervals ("UNION" or "INTERSECTION").
    - inverted_read_filter: List of inverted read filters.
    - lenient: Lenient processing of VCF files.
    - max_variants_per_shard: Partition VCF output into shards with max records.
    - quiet: Suppress job-summary info on stderr.
    - read_filter: List of read filters to apply.
    - read_index: List of indices for read inputs.
    - read_validation_stringency: Validation stringency ("STRICT", "LENIENT", "SILENT").
    - seconds_between_progress_updates: Interval in seconds for progress updates.
    - sequence_dictionary: Sequence dictionary file (.dict).
    - show_hidden: Display hidden arguments.
    - tmp_dir: Temporary directory path.
    - use_jdk_deflater: Use JdkDeflater compression.
    - use_jdk_inflater: Use JdkInflater decompression.
    - verbosity: Logging verbosity level ("ERROR", "WARNING", "INFO", "DEBUG").
    - version: Display version number.
    
    Returns:
    A dictionary with keys: command_executed, stdout, stderr, output_files.
    """
    # Validate input files
    if not input_reads:
        raise ValueError("At least one input read file must be specified.")
    for f in input_reads:
        if not f.exists():
            raise FileNotFoundError(f"Input read file not found: {f}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if not known_sites:
        raise ValueError("At least one known sites file must be specified.")
    for f in known_sites:
        if not f.exists():
            raise FileNotFoundError(f"Known sites file not found: {f}")
    if output.exists():
        # Allow overwrite but warn? Here we just proceed.
        pass
    if arguments_file:
        for f in arguments_file:
            if not f.exists():
                raise FileNotFoundError(f"Arguments file not found: {f}")
    if read_index:
        if len(read_index) != len(input_reads):
            raise ValueError("Number of read indices must match number of input reads.")
        for f in read_index:
            if not f.exists():
                raise FileNotFoundError(f"Read index file not found: {f}")
    if sequence_dictionary and not sequence_dictionary.exists():
        raise FileNotFoundError(f"Sequence dictionary file not found: {sequence_dictionary}")
    if tmp_dir and not tmp_dir.exists():
        raise FileNotFoundError(f"Temporary directory not found: {tmp_dir}")
    if interval_merging_rule not in ("ALL", "OVERLAPPING_ONLY"):
        raise ValueError("interval_merging_rule must be 'ALL' or 'OVERLAPPING_ONLY'")
    if interval_set_rule not in ("UNION", "INTERSECTION"):
        raise ValueError("interval_set_rule must be 'UNION' or 'INTERSECTION'")
    if read_validation_stringency not in ("STRICT", "LENIENT", "SILENT"):
        raise ValueError("read_validation_stringency must be one of 'STRICT', 'LENIENT', 'SILENT'")
    if verbosity not in ("ERROR", "WARNING", "INFO", "DEBUG"):
        raise ValueError("verbosity must be one of 'ERROR', 'WARNING', 'INFO', 'DEBUG'")
    if indels_context_size < 1 or indels_context_size > 13:
        raise ValueError("indels_context_size must be between 1 and 13 inclusive")
    if mismatches_context_size < 1 or mismatches_context_size > 13:
        raise ValueError("mismatches_context_size must be between 1 and 13 inclusive")
    if max_variants_per_shard < 0:
        raise ValueError("max_variants_per_shard must be zero or positive")
    if low_quality_tail < 0:
        raise ValueError("low_quality_tail must be zero or positive")
    if preserve_qscores_less_than < 0:
        raise ValueError("preserve_qscores_less_than must be zero or positive")
    if quantizing_levels < 1:
        raise ValueError("quantizing_levels must be positive")
    if maximum_cycle_value < 0:
        raise ValueError("maximum_cycle_value must be zero or positive")
    if cloud_index_prefetch_buffer < -1:
        raise ValueError("cloud_index_prefetch_buffer must be >= -1")
    if cloud_prefetch_buffer < 0:
        raise ValueError("cloud_prefetch_buffer must be zero or positive")
    if gcs_max_retries < 0:
        raise ValueError("gcs_max_retries must be zero or positive")
    if gcs_retries < 0:
        raise ValueError("gcs_retries must be zero or positive")
    if interval_exclusion_padding < 0:
        raise ValueError("interval_exclusion_padding must be zero or positive")
    if interval_padding < 0:
        raise ValueError("interval_padding must be zero or positive")
    if bqsr_baq_gap_open_penalty < 0:
        raise ValueError("bqsr_baq_gap_open_penalty must be non-negative")

    # Build command line
    cmd = ["gatk", "BaseRecalibrator"]

    # Input reads and indices
    for i, read_file in enumerate(input_reads):
        cmd.extend(["-I", str(read_file)])
        if read_index and read_index[i]:
            cmd.extend(["--read-index", str(read_index[i])])

    # Reference
    cmd.extend(["-R", str(reference)])

    # Known sites (multiple allowed)
    for ks in known_sites:
        cmd.extend(["--known-sites", str(ks)])

    # Output
    cmd.extend(["-O", str(output)])

    # Optional arguments
    if arguments_file:
        for argf in arguments_file:
            cmd.extend(["--arguments_file", str(argf)])
    if binary_tag_name:
        cmd.extend(["--binary-tag-name", binary_tag_name])
    cmd.extend(["--bqsr-baq-gap-open-penalty", str(bqsr_baq_gap_open_penalty)])
    cmd.extend(["--cloud-index-prefetch-buffer", str(cloud_index_prefetch_buffer)])
    cmd.extend(["--cloud-prefetch-buffer", str(cloud_prefetch_buffer)])
    cmd.extend(["--default-base-qualities", str(default_base_qualities)])
    cmd.extend(["--deletions-default-quality", str(deletions_default_quality)])
    if disable_bam_index_caching:
        cmd.append("--disable-bam-index-caching")
    if disable_sequence_dictionary_validation:
        cmd.append("--disable-sequence-dictionary-validation")
    cmd.extend(["--gcs-max-retries", str(gcs_max_retries)])
    if gcs_project_for_requester_pays:
        cmd.extend(["--gcs-project-for-requester-pays", gcs_project_for_requester_pays])
    cmd.extend(["--indels-context-size", str(indels_context_size)])
    cmd.extend(["--insertions-default-quality", str(insertions_default_quality)])
    cmd.extend(["--interval-merging-rule", interval_merging_rule])
    if intervals:
        for interval in intervals:
            cmd.extend(["-L", interval])
    cmd.extend(["--low-quality-tail", str(low_quality_tail)])
    cmd.extend(["--maximum-cycle-value", str(maximum_cycle_value)])
    cmd.extend(["--mismatches-context-size", str(mismatches_context_size)])
    cmd.extend(["--mismatches-default-quality", str(mismatches_default_quality)])
    cmd.extend(["--preserve-qscores-less-than", str(preserve_qscores_less_than)])
    cmd.extend(["--quantizing-levels", str(quantizing_levels)])
    if sites_only_vcf_output:
        cmd.append("--sites-only-vcf-output")
    if use_original_qualities:
        cmd.append("--use-original-qualities")
    if add_output_sam_program_record:
        cmd.append("--add-output-sam-program-record")
    if add_output_vcf_command_line:
        cmd.append("--add-output-vcf-command-line")
    if create_output_bam_index:
        cmd.append("--create-output-bam-index")
    if create_output_bam_md5:
        cmd.append("--create-output-bam-md5")
    if create_output_variant_index:
        cmd.append("--create-output-variant-index")
    if create_output_variant_md5:
        cmd.append("--create-output-variant-md5")
    if disable_read_filter:
        for filt in disable_read_filter:
            cmd.extend(["--disable-read-filter", filt])
    if disable_tool_default_read_filters:
        cmd.append("--disable-tool-default-read-filters")
    if exclude_intervals:
        for excl in exclude_intervals:
            cmd.extend(["-XL", excl])
    if gatk_config_file:
        cmd.extend(["--gatk-config-file", gatk_config_file])
    # gcs_retries alias
    cmd.extend(["--gcs-retries", str(gcs_retries)])
    if help_flag:
        cmd.append("--help")
    cmd.extend(["--interval-exclusion-padding", str(interval_exclusion_padding)])
    cmd.extend(["--interval-padding", str(interval_padding)])
    cmd.extend(["--interval-set-rule", interval_set_rule])
    if inverted_read_filter:
        for invf in inverted_read_filter:
            cmd.extend(["--inverted-read-filter", invf])
    if lenient:
        cmd.append("--lenient")
    if max_variants_per_shard > 0:
        cmd.extend(["--max-variants-per-shard", str(max_variants_per_shard)])
    if quiet:
        cmd.append("--QUIET")
    if read_filter:
        for rf in read_filter:
            cmd.extend(["--read-filter", rf])
    if read_validation_stringency:
        cmd.extend(["--read-validation-stringency", read_validation_stringency])
    cmd.extend(["--seconds-between-progress-updates", str(seconds_between_progress_updates)])
    if sequence_dictionary:
        cmd.extend(["--sequence-dictionary", str(sequence_dictionary)])
    if show_hidden:
        cmd.append("--showHidden")
    if tmp_dir:
        cmd.extend(["--tmp-dir", str(tmp_dir)])
    if use_jdk_deflater:
        cmd.append("--use-jdk-deflater")
    if use_jdk_inflater:
        cmd.append("--use-jdk-inflater")
    if verbosity:
        cmd.extend(["--verbosity", verbosity])
    if version:
        cmd.append("--version")

    # Run command
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
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