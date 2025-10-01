from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Set
import subprocess

mcp = FastMCP()

@mcp.tool()
def gatk_SelectVariants(
    reference: Path,
    variant: Path,
    output: Path,
    apply_jexl_filters_first: bool = False,
    arguments_file: Optional[List[Path]] = None,
    call_genotypes: bool = False,
    cloud_index_prefetch_buffer: int = -1,
    cloud_prefetch_buffer: int = 40,
    concordance: Optional[Path] = None,
    disable_bam_index_caching: bool = False,
    disable_sequence_dictionary_validation: bool = False,
    discordance: Optional[Path] = None,
    drop_genotype_annotation: Optional[List[str]] = None,
    drop_info_annotation: Optional[List[str]] = None,
    exclude_filtered: bool = False,
    exclude_ids: Optional[Set[str]] = None,
    exclude_non_variants: bool = False,
    exclude_sample_expressions: Optional[Set[str]] = None,
    exclude_sample_name: Optional[Set[str]] = None,
    gcs_max_retries: int = 20,
    gcs_project_for_requester_pays: str = "",
    genomicsdb_max_alternate_alleles: int = 50,
    genomicsdb_shared_posixfs_optimizations: bool = False,
    genomicsdb_use_bcf_codec: bool = False,
    genomicsdb_use_gcs_hdfs_connector: bool = False,
    help_flag: bool = False,
    ignore_non_ref_in_types: bool = False,
    input_reads: Optional[List[Path]] = None,
    interval_exclusion_padding: int = 0,
    interval_merging_rule: str = "ALL",
    interval_padding: int = 0,
    interval_set_rule: str = "UNION",
    intervals: Optional[List[str]] = None,
    invert_mendelian_violation: bool = False,
    inverted_read_filter: Optional[List[str]] = None,
    invert_select: bool = False,
    keep_ids: Optional[Set[str]] = None,
    keep_original_ac: bool = False,
    keep_original_dp: bool = False,
    lenient: bool = False,
    max_filtered_genotypes: int = 2147483647,
    max_fraction_filtered_genotypes: float = 1.0,
    max_indel_size: int = 2147483647,
    max_nocall_fraction: float = 1.0,
    max_nocall_number: int = 2147483647,
    max_variants_per_shard: int = 0,
    mendelian_violation: bool = False,
    mendelian_violation_qual_threshold: float = 0.0,
    min_filtered_genotypes: int = 0,
    min_fraction_filtered_genotypes: float = 0.0,
    min_indel_size: int = 0,
    pedigree: Optional[Path] = None,
    preserve_alleles: bool = False,
    quiet: bool = False,
    read_filter: Optional[List[str]] = None,
    read_index: Optional[List[Path]] = None,
    read_validation_stringency: str = "SILENT",
    remove_fraction_genotypes: float = 0.0,
    remove_unused_alternates: bool = False,
    restrict_alleles_to: str = "ALL",
    sample_expressions: Optional[Set[str]] = None,
    sample_name: Optional[Set[str]] = None,
    select: Optional[List[str]] = None,
    select_genotype_expressions: Optional[List[str]] = None,
    select_random_fraction: float = 1.0,
    select_type_to_exclude: Optional[List[str]] = None,
    select_type_to_include: Optional[List[str]] = None,
    set_filtered_gt_to_nocall: bool = False,
    sites_only_vcf_output: bool = False,
    version: bool = False,
    add_output_sam_program_record: bool = True,
    add_output_vcf_command_line: bool = True,
    create_output_bam_index: bool = True,
    create_output_bam_md5: bool = False,
    create_output_variant_index: bool = True,
    create_output_variant_md5: bool = False,
    disable_read_filter: Optional[List[str]] = None,
    disable_tool_default_read_filters: bool = False,
    exclude_intervals: Optional[List[str]] = None,
    gatk_config_file: Optional[Path] = None,
    seconds_between_progress_updates: float = 10.0,
    sequence_dictionary: Optional[Path] = None,
    tmp_dir: Optional[Path] = None,
    use_jdk_deflater: bool = False,
    use_jdk_inflater: bool = False,
    verbosity: str = "INFO",
    variant_output_filtering: Optional[str] = None,
    show_hidden: bool = False,
) -> dict:
    """
    Select a subset of variants from a VCF file using GATK SelectVariants.

    Parameters:
    - reference: Reference genome fasta file (required).
    - variant: Input VCF file containing variants (required).
    - output: Output VCF file path (required).
    - apply_jexl_filters_first: Apply JEXL-based filtering before subsetting samples.
    - arguments_file: List of argument files to add to the command line.
    - call_genotypes: Output called genotypes in final VCF.
    - cloud_index_prefetch_buffer: Size of cloud-only prefetch buffer in MB (-1 to default).
    - cloud_prefetch_buffer: Size of cloud-only prefetch buffer in MB.
    - concordance: VCF file for concordance comparison.
    - disable_bam_index_caching: Disable BAM index caching.
    - disable_sequence_dictionary_validation: Disable sequence dictionary validation.
    - discordance: VCF file for discordance comparison.
    - drop_genotype_annotation: List of genotype annotation keys to drop.
    - drop_info_annotation: List of info annotation keys to drop.
    - exclude_filtered: Exclude filtered sites.
    - exclude_ids: Set of variant rsIDs to exclude.
    - exclude_non_variants: Exclude non-variant sites.
    - exclude_sample_expressions: Set of sample regex expressions to exclude.
    - exclude_sample_name: Set of sample names to exclude.
    - gcs_max_retries: Max retries for GCS bucket errors.
    - gcs_project_for_requester_pays: GCS project for requester pays buckets.
    - genomicsdb_max_alternate_alleles: Max alternate alleles for GenomicsDB.
    - genomicsdb_shared_posixfs_optimizations: Enable GenomicsDB POSIX FS optimizations.
    - genomicsdb_use_bcf_codec: Use BCF codec for GenomicsDB.
    - genomicsdb_use_gcs_hdfs_connector: Use GCS HDFS connector.
    - help_flag: Display help message.
    - ignore_non_ref_in_types: Ignore NON_REF alleles in type determination.
    - input_reads: List of BAM/SAM/CRAM files containing reads.
    - interval_exclusion_padding: Padding bp for excluded intervals.
    - interval_merging_rule: Interval merging rule ("ALL" or "OVERLAPPING_ONLY").
    - interval_padding: Padding bp for included intervals.
    - interval_set_rule: Interval set merging rule ("UNION" or "INTERSECTION").
    - intervals: List of genomic intervals to operate on.
    - invert_mendelian_violation: Output non-mendelian violation sites only.
    - inverted_read_filter: List of inverted read filters.
    - invert_select: Invert selection criteria.
    - keep_ids: Set of variant rsIDs to keep.
    - keep_original_ac: Keep original AC, AF, AN annotations.
    - keep_original_dp: Keep original DP annotation.
    - lenient: Lenient VCF processing.
    - max_filtered_genotypes: Max number of filtered genotypes allowed.
    - max_fraction_filtered_genotypes: Max fraction of filtered genotypes allowed.
    - max_indel_size: Max indel size to include.
    - max_nocall_fraction: Max fraction of no-call genotypes allowed.
    - max_nocall_number: Max number of no-call genotypes allowed.
    - max_variants_per_shard: Max variants per shard (0 disables sharding).
    - mendelian_violation: Output mendelian violation sites only.
    - mendelian_violation_qual_threshold: Minimum GQ for mendelian violation.
    - min_filtered_genotypes: Min number of filtered genotypes required.
    - min_fraction_filtered_genotypes: Min fraction of filtered genotypes required.
    - min_indel_size: Min indel size to include.
    - pedigree: Pedigree file path.
    - preserve_alleles: Preserve original alleles, do not trim.
    - quiet: Suppress job-summary info.
    - read_filter: List of read filters to apply.
    - read_index: List of indices for read inputs.
    - read_validation_stringency: Validation stringency ("STRICT", "LENIENT", "SILENT").
    - remove_fraction_genotypes: Fraction of genotypes to randomly set to no-call.
    - remove_unused_alternates: Remove alternate alleles not present in genotypes.
    - restrict_alleles_to: Allelicity restriction ("ALL", "BIALLELIC", "MULTIALLELIC").
    - sample_expressions: Set of regex expressions to select samples.
    - sample_name: Set of sample names to include.
    - select: List of JEXL expressions to filter variants.
    - select_genotype_expressions: List of JEXL expressions to filter genotypes.
    - select_random_fraction: Fraction of variants to select randomly.
    - select_type_to_exclude: List of variant types to exclude.
    - select_type_to_include: List of variant types to include.
    - set_filtered_gt_to_nocall: Set filtered genotypes to no-call.
    - sites_only_vcf_output: Output sites-only VCF (no genotype fields).
    - version: Display version number.
    - add_output_sam_program_record: Add PG tag to created SAM/BAM/CRAM files.
    - add_output_vcf_command_line: Add command line header to VCF files.
    - create_output_bam_index: Create BAM/CRAM index.
    - create_output_bam_md5: Create MD5 digest for BAM/SAM/CRAM files.
    - create_output_variant_index: Create VCF index.
    - create_output_variant_md5: Create MD5 digest for VCF files.
    - disable_read_filter: List of read filters to disable.
    - disable_tool_default_read_filters: Disable all tool default read filters.
    - exclude_intervals: List of genomic intervals to exclude.
    - gatk_config_file: GATK configuration file.
    - seconds_between_progress_updates: Seconds between progress updates.
    - sequence_dictionary: Sequence dictionary file (.dict).
    - tmp_dir: Temporary directory.
    - use_jdk_deflater: Use JdkDeflater compression.
    - use_jdk_inflater: Use JdkInflater decompression.
    - verbosity: Logging verbosity level ("ERROR", "WARNING", "INFO", "DEBUG").
    - variant_output_filtering: Variant output filtering mode.
    - show_hidden: Show hidden arguments.

    Returns:
    dict with keys: command_executed, stdout, stderr, output_files
    """
    import shlex

    # Validate required files
    if not reference.is_file():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if not variant.is_file():
        raise FileNotFoundError(f"Variant file not found: {variant}")
    if pedigree is not None and not pedigree.is_file():
        raise FileNotFoundError(f"Pedigree file not found: {pedigree}")
    if gatk_config_file is not None and not gatk_config_file.is_file():
        raise FileNotFoundError(f"GATK config file not found: {gatk_config_file}")
    if sequence_dictionary is not None and not sequence_dictionary.is_file():
        raise FileNotFoundError(f"Sequence dictionary file not found: {sequence_dictionary}")
    if tmp_dir is not None and not tmp_dir.is_dir():
        raise NotADirectoryError(f"Temporary directory not found: {tmp_dir}")

    # Validate numeric ranges
    if not (0.0 <= remove_fraction_genotypes <= 1.0):
        raise ValueError("remove_fraction_genotypes must be between 0 and 1")
    if not (0.0 <= max_fraction_filtered_genotypes <= 1.0):
        raise ValueError("max_fraction_filtered_genotypes must be between 0 and 1")
    if not (0.0 <= min_fraction_filtered_genotypes <= 1.0):
        raise ValueError("min_fraction_filtered_genotypes must be between 0 and 1")
    if not (0.0 <= max_nocall_fraction <= 1.0):
        raise ValueError("max_nocall_fraction must be between 0 and 1")
    if select_random_fraction <= 0.0 or select_random_fraction > 1.0:
        raise ValueError("select_random_fraction must be > 0 and <= 1")

    # Validate enums
    valid_interval_merging_rules = {"ALL", "OVERLAPPING_ONLY"}
    if interval_merging_rule not in valid_interval_merging_rules:
        raise ValueError(f"interval_merging_rule must be one of {valid_interval_merging_rules}")

    valid_interval_set_rules = {"UNION", "INTERSECTION"}
    if interval_set_rule not in valid_interval_set_rules:
        raise ValueError(f"interval_set_rule must be one of {valid_interval_set_rules}")

    valid_read_validation_stringency = {"STRICT", "LENIENT", "SILENT"}
    if read_validation_stringency not in valid_read_validation_stringency:
        raise ValueError(f"read_validation_stringency must be one of {valid_read_validation_stringency}")

    valid_restrict_alleles_to = {"ALL", "BIALLELIC", "MULTIALLELIC"}
    if restrict_alleles_to not in valid_restrict_alleles_to:
        raise ValueError(f"restrict_alleles_to must be one of {valid_restrict_alleles_to}")

    valid_verbosity = {"ERROR", "WARNING", "INFO", "DEBUG"}
    if verbosity not in valid_verbosity:
        raise ValueError(f"verbosity must be one of {valid_verbosity}")

    valid_variant_output_filtering = {None, "STARTS_IN", "ENDS_IN", "OVERLAPS", "CONTAINED", "ANYWHERE"}
    if variant_output_filtering not in valid_variant_output_filtering:
        raise ValueError(f"variant_output_filtering must be one of {valid_variant_output_filtering}")

    valid_select_types = {"NO_VARIATION", "SNP", "MNP", "INDEL", "SYMBOLIC", "MIXED"}

    if select_type_to_include:
        for t in select_type_to_include:
            if t not in valid_select_types:
                raise ValueError(f"select_type_to_include contains invalid type: {t}")

    if select_type_to_exclude:
        for t in select_type_to_exclude:
            if t not in valid_select_types:
                raise ValueError(f"select_type_to_exclude contains invalid type: {t}")

    # Prepare command line
    cmd = ["gatk", "SelectVariants"]
    cmd += ["-R", str(reference)]
    cmd += ["-V", str(variant)]
    cmd += ["-O", str(output)]

    if apply_jexl_filters_first:
        cmd.append("--apply-jexl-filters-first")

    if arguments_file:
        for af in arguments_file:
            if not af.is_file():
                raise FileNotFoundError(f"Arguments file not found: {af}")
            cmd += ["--arguments_file", str(af)]

    if call_genotypes:
        cmd.append("--call-genotypes")

    if cloud_index_prefetch_buffer != -1:
        cmd += ["--cloud-index-prefetch-buffer", str(cloud_index_prefetch_buffer)]

    if cloud_prefetch_buffer != 40:
        cmd += ["--cloud-prefetch-buffer", str(cloud_prefetch_buffer)]

    if concordance:
        if not Path(concordance).exists():
            raise FileNotFoundError(f"Concordance file not found: {concordance}")
        cmd += ["--concordance", str(concordance)]

    if disable_bam_index_caching:
        cmd.append("--disable-bam-index-caching")

    if disable_sequence_dictionary_validation:
        cmd.append("--disable-sequence-dictionary-validation")

    if discordance:
        if not Path(discordance).exists():
            raise FileNotFoundError(f"Discordance file not found: {discordance}")
        cmd += ["--discordance", str(discordance)]

    if drop_genotype_annotation:
        for dga in drop_genotype_annotation:
            cmd += ["--drop-genotype-annotation", dga]

    if drop_info_annotation:
        for dia in drop_info_annotation:
            cmd += ["--drop-info-annotation", dia]

    if exclude_filtered:
        cmd.append("--exclude-filtered")

    if exclude_ids:
        for eid in exclude_ids:
            if eid.endswith(".list"):
                path = Path(eid)
                if not path.is_file():
                    raise FileNotFoundError(f"Exclude IDs list file not found: {eid}")
                cmd += ["--exclude-ids", str(path)]
            else:
                cmd += ["--exclude-ids", eid]

    if exclude_non_variants:
        cmd.append("--exclude-non-variants")

    if exclude_sample_expressions:
        for ese in exclude_sample_expressions:
            cmd += ["--exclude-sample-expressions", ese]

    if exclude_sample_name:
        for esn in exclude_sample_name:
            if esn.endswith(".args"):
                path = Path(esn)
                if not path.is_file():
                    raise FileNotFoundError(f"Exclude sample name args file not found: {esn}")
                cmd += ["--exclude-sample-name", str(path)]
            else:
                cmd += ["--exclude-sample-name", esn]

    if gcs_max_retries != 20:
        cmd += ["--gcs-max-retries", str(gcs_max_retries)]

    if gcs_project_for_requester_pays:
        cmd += ["--gcs-project-for-requester-pays", gcs_project_for_requester_pays]

    if genomicsdb_max_alternate_alleles != 50:
        cmd += ["--genomicsdb-max-alternate-alleles", str(genomicsdb_max_alternate_alleles)]

    if genomicsdb_shared_posixfs_optimizations:
        cmd.append("--genomicsdb-shared-posixfs-optimizations")

    if genomicsdb_use_bcf_codec:
        cmd.append("--genomicsdb-use-bcf-codec")

    if genomicsdb_use_gcs_hdfs_connector:
        cmd.append("--genomicsdb-use-gcs-hdfs-connector")

    if help_flag:
        cmd.append("--help")

    if ignore_non_ref_in_types:
        cmd.append("--ignore-non-ref-in-types")

    if input_reads:
        for ir in input_reads:
            if not Path(ir).is_file():
                raise FileNotFoundError(f"Input read file not found: {ir}")
            cmd += ["--input", str(ir)]

    if interval_exclusion_padding != 0:
        cmd += ["--interval-exclusion-padding", str(interval_exclusion_padding)]

    if interval_merging_rule != "ALL":
        cmd += ["--interval-merging-rule", interval_merging_rule]

    if interval_padding != 0:
        cmd += ["--interval-padding", str(interval_padding)]

    if interval_set_rule != "UNION":
        cmd += ["--interval-set-rule", interval_set_rule]

    if intervals:
        for interval in intervals:
            cmd += ["--intervals", interval]

    if invert_mendelian_violation:
        cmd.append("--invert-mendelian-violation")

    if inverted_read_filter:
        for irf in inverted_read_filter:
            cmd += ["--inverted-read-filter", irf]

    if invert_select:
        cmd.append("--invertSelect")

    if keep_ids:
        for kid in keep_ids:
            if kid.endswith(".list"):
                path = Path(kid)
                if not path.is_file():
                    raise FileNotFoundError(f"Keep IDs list file not found: {kid}")
                cmd += ["--keep-ids", str(path)]
            else:
                cmd += ["--keep-ids", kid]

    if keep_original_ac:
        cmd.append("--keep-original-ac")

    if keep_original_dp:
        cmd.append("--keep-original-dp")

    if lenient:
        cmd.append("--lenient")

    if max_filtered_genotypes != 2147483647:
        cmd += ["--max-filtered-genotypes", str(max_filtered_genotypes)]

    if max_fraction_filtered_genotypes != 1.0:
        cmd += ["--max-fraction-filtered-genotypes", str(max_fraction_filtered_genotypes)]

    if max_indel_size != 2147483647:
        cmd += ["--max-indel-size", str(max_indel_size)]

    if max_nocall_fraction != 1.0:
        cmd += ["--max-nocall-fraction", str(max_nocall_fraction)]

    if max_nocall_number != 2147483647:
        cmd += ["--max-nocall-number", str(max_nocall_number)]

    if max_variants_per_shard != 0:
        cmd += ["--max-variants-per-shard", str(max_variants_per_shard)]

    if mendelian_violation:
        cmd.append("--mendelian-violation")

    if mendelian_violation_qual_threshold != 0.0:
        cmd += ["--mendelian-violation-qual-threshold", str(mendelian_violation_qual_threshold)]

    if min_filtered_genotypes != 0:
        cmd += ["--min-filtered-genotypes", str(min_filtered_genotypes)]

    if min_fraction_filtered_genotypes != 0.0:
        cmd += ["--min-fraction-filtered-genotypes", str(min_fraction_filtered_genotypes)]

    if min_indel_size != 0:
        cmd += ["--min-indel-size", str(min_indel_size)]

    if pedigree:
        cmd += ["--pedigree", str(pedigree)]

    if preserve_alleles:
        cmd.append("--preserve-alleles")

    if quiet:
        cmd.append("--QUIET")

    if read_filter:
        for rf in read_filter:
            cmd += ["--read-filter", rf]

    if read_index:
        for ri in read_index:
            if not Path(ri).is_file():
                raise FileNotFoundError(f"Read index file not found: {ri}")
            cmd += ["--read-index", str(ri)]

    if read_validation_stringency != "SILENT":
        cmd += ["--read-validation-stringency", read_validation_stringency]

    if remove_fraction_genotypes != 0.0:
        cmd += ["--remove-fraction-genotypes", str(remove_fraction_genotypes)]

    if remove_unused_alternates:
        cmd.append("--remove-unused-alternates")

    if restrict_alleles_to != "ALL":
        cmd += ["--restrict-alleles-to", restrict_alleles_to]

    if sample_expressions:
        for se in sample_expressions:
            cmd += ["--sample-expressions", se]

    if sample_name:
        for sn in sample_name:
            if sn.endswith(".args"):
                path = Path(sn)
                if not path.is_file():
                    raise FileNotFoundError(f"Sample name args file not found: {sn}")
                cmd += ["--sample-name", str(path)]
            else:
                cmd += ["--sample-name", sn]

    if select:
        for s in select:
            cmd += ["--select", s]

    if select_genotype_expressions:
        for sge in select_genotype_expressions:
            cmd += ["--select-genotype", sge]

    if select_random_fraction != 1.0:
        cmd += ["--select-random-fraction", str(select_random_fraction)]

    if select_type_to_exclude:
        for ste in select_type_to_exclude:
            cmd += ["--select-type-to-exclude", ste]

    if select_type_to_include:
        for sti in select_type_to_include:
            cmd += ["--select-type-to-include", sti]

    if set_filtered_gt_to_nocall:
        cmd.append("--set-filtered-gt-to-nocall")

    if sites_only_vcf_output:
        cmd.append("--sites-only-vcf-output")

    if version:
        cmd.append("--version")

    if add_output_sam_program_record is False:
        cmd.append("--add-output-sam-program-record=false")
    else:
        cmd.append("--add-output-sam-program-record=true")

    if add_output_vcf_command_line is False:
        cmd.append("--add-output-vcf-command-line=false")
    else:
        cmd.append("--add-output-vcf-command-line=true")

    if create_output_bam_index is False:
        cmd.append("--create-output-bam-index=false")
    else:
        cmd.append("--create-output-bam-index=true")

    if create_output_bam_md5:
        cmd.append("--create-output-bam-md5")

    if create_output_variant_index is False:
        cmd.append("--create-output-variant-index=false")
    else:
        cmd.append("--create-output-variant-index=true")

    if create_output_variant_md5:
        cmd.append("--create-output-variant-md5")

    if disable_read_filter:
        for drf in disable_read_filter:
            cmd += ["--disable-read-filter", drf]

    if disable_tool_default_read_filters:
        cmd.append("--disable-tool-default-read-filters")

    if exclude_intervals:
        for exi in exclude_intervals:
            cmd += ["--exclude-intervals", exi]

    if gatk_config_file:
        cmd += ["--gatk-config-file", str(gatk_config_file)]

    if seconds_between_progress_updates != 10.0:
        cmd += ["--seconds-between-progress-updates", str(seconds_between_progress_updates)]

    if sequence_dictionary:
        cmd += ["--sequence-dictionary", str(sequence_dictionary)]

    if tmp_dir:
        cmd += ["--tmp-dir", str(tmp_dir)]

    if use_jdk_deflater:
        cmd.append("--use-jdk-deflater")

    if use_jdk_inflater:
        cmd.append("--use-jdk-inflater")

    if verbosity != "INFO":
        cmd += ["--verbosity", verbosity]

    if variant_output_filtering:
        cmd += ["--variant-output-filtering", variant_output_filtering]

    if show_hidden:
        cmd.append("--showHidden")

    # Run command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(shlex.quote(c) for c in cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": []
        }

    output_files = [str(output)]
    # If create_output_variant_index is True, expect .tbi or .idx index file
    if create_output_variant_index:
        idx_path = output.with_suffix(output.suffix + ".tbi")
        if idx_path.exists():
            output_files.append(str(idx_path))
        else:
            # Some VCF index files use .idx extension
            idx_alt = output.with_suffix(".idx")
            if idx_alt.exists():
                output_files.append(str(idx_alt))

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()