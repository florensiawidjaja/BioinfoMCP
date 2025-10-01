from fastmcp import FastMCP
from typing import List, Optional
from pathlib import Path
import subprocess

mcp = FastMCP()

@mcp.tool()
def gatk_HaplotypeCaller(
    input_bams: List[Path],
    output: Path,
    reference: Path,
    alleles: Optional[Path] = None,
    annotate_with_num_discovered_alleles: bool = False,
    annotation: Optional[List[str]] = None,
    annotation_group: Optional[List[str]] = None,
    annotations_to_exclude: Optional[List[str]] = None,
    arguments_file: Optional[List[Path]] = None,
    assembly_region_out: Optional[Path] = None,
    assembly_region_padding: int = 100,
    base_quality_score_threshold: int = 18,
    cloud_index_prefetch_buffer: int = -1,
    cloud_prefetch_buffer: int = 40,
    contamination_fraction_to_filter: float = 0.0,
    dbsnp: Optional[Path] = None,
    disable_bam_index_caching: bool = False,
    disable_sequence_dictionary_validation: bool = False,
    dont_use_dragstr_pair_hmm_scores: bool = False,
    dont_use_soft_clipped_bases: bool = False,
    dragen_mode: bool = False,
    dragstr_het_hom_ratio: int = 2,
    dragstr_params_path: Optional[Path] = None,
    enable_dynamic_read_disqualification_for_genotyping: bool = False,
    flow_order_for_annotations: Optional[List[str]] = None,
    founder_id: Optional[List[str]] = None,
    gcs_max_retries: int = 20,
    gcs_project_for_requester_pays: str = "",
    genotype_assignment_method: str = "USE_PLS_TO_ASSIGN",
    graph_output: Optional[Path] = None,
    help: bool = False,
    heterozygosity: float = 0.001,
    heterozygosity_stdev: float = 0.01,
    indel_heterozygosity: float = 0.000125,
    interval_merging_rule: str = "ALL",
    intervals: Optional[List[str]] = None,
    interval_exclusion_padding: int = 0,
    interval_padding: int = 0,
    interval_set_rule: str = "UNION",
    inverted_read_filter: Optional[List[str]] = None,
    keep_boundary_flows: bool = False,
    kmer_size: Optional[List[int]] = None,
    lenient: bool = False,
    likelihood_calculation_engine: str = "PairHMM",
    linked_de_bruijn_graph: bool = False,
    mapping_quality_threshold_for_genotyping: int = 20,
    max_alternate_alleles: int = 6,
    max_assembly_region_size: int = 300,
    max_effective_depth_adjustment_for_frd: int = 0,
    max_genotype_count: int = 1024,
    max_mnp_distance: int = 0,
    max_num_haplotypes_in_population: int = 128,
    max_prob_propagation_distance: int = 50,
    max_reads_per_alignment_start: int = 50,
    max_unpruned_variants: int = 100,
    max_variants_per_shard: int = 0,
    min_assembly_region_size: int = 50,
    min_base_quality_score: int = 10,
    min_dangling_branch_length: int = 4,
    min_pruning: int = 2,
    native_pair_hmm_threads: int = 4,
    native_pair_hmm_use_double_precision: bool = False,
    num_pruning_samples: int = 1,
    num_reference_samples_if_no_call: int = 0,
    output_mode: str = "EMIT_VARIANTS_ONLY",
    pedigree: Optional[Path] = None,
    ploidy_regions: Optional[Path] = None,
    population_callset: Optional[Path] = None,
    pcr_indel_model: str = "CONSERVATIVE",
    phred_scaled_global_read_mismapping_rate: int = 45,
    pileup_detection: bool = False,
    ploidy: int = 2,
    read_filter: Optional[List[str]] = None,
    read_index: Optional[List[Path]] = None,
    read_validation_stringency: str = "SILENT",
    recover_all_dangling_branches: bool = False,
    reference_model_deletion_quality: int = 30,
    sample_name: Optional[str] = None,
    seconds_between_progress_updates: float = 10.0,
    sequence_dictionary: Optional[Path] = None,
    show_hidden: bool = False,
    sites_only_vcf_output: bool = False,
    smith_waterman: str = "FASTEST_AVAILABLE",
    smith_waterman_dangling_end_gap_extend_penalty: int = -6,
    smith_waterman_dangling_end_gap_open_penalty: int = -110,
    smith_waterman_dangling_end_match_value: int = 25,
    smith_waterman_dangling_end_mismatch_penalty: int = -50,
    smith_waterman_haplotype_to_reference_gap_extend_penalty: int = -11,
    smith_waterman_haplotype_to_reference_gap_open_penalty: int = -260,
    smith_waterman_haplotype_to_reference_match_value: int = 200,
    smith_waterman_haplotype_to_reference_mismatch_penalty: int = -150,
    smith_waterman_read_to_haplotype_gap_extend_penalty: int = -5,
    smith_waterman_read_to_haplotype_gap_open_penalty: int = -30,
    smith_waterman_read_to_haplotype_match_value: int = 10,
    smith_waterman_read_to_haplotype_mismatch_penalty: int = -15,
    soft_clip_low_quality_ends: bool = False,
    standard_min_confidence_threshold_for_calling: float = 30.0,
    tmp_dir: Optional[Path] = None,
    transform_dragen_mapping_quality: bool = False,
    use_filtered_reads_for_annotations: bool = False,
    use_jdk_deflater: bool = False,
    use_jdk_inflater: bool = False,
    use_new_qual_calculator: bool = True,
    use_pdhmm: bool = False,
    use_pdhmm_overlap_optimization: bool = False,
    use_posteriors_to_calculate_qual: bool = False,
    verbosity: str = "INFO",
    version: bool = False,
):
    """
    Call germline SNPs and indels via local re-assembly of haplotypes using GATK HaplotypeCaller.
    
    Parameters:
    - input_bams: List of input BAM/SAM/CRAM files containing reads.
    - output: Output VCF or GVCF file path.
    - reference: Reference sequence FASTA file path.
    - alleles: Optional VCF file with alleles to force-call.
    - annotate_with_num_discovered_alleles: Annotate records with number of discovered alternate alleles.
    - annotation: List of specific annotations to add.
    - annotation_group: List of annotation groups to add.
    - annotations_to_exclude: List of annotations to exclude.
    - arguments_file: List of argument files to include.
    - assembly_region_out: Output file for assembly regions.
    - assembly_region_padding: Padding bases around assembly regions.
    - base_quality_score_threshold: Bases below this quality are reduced to minimum.
    - cloud_index_prefetch_buffer: Size of cloud-only prefetch buffer (MB).
    - cloud_prefetch_buffer: Size of cloud-only prefetch buffer (MB).
    - contamination_fraction_to_filter: Fraction of contamination to filter.
    - dbsnp: dbSNP VCF file.
    - disable_bam_index_caching: Disable BAM index caching.
    - disable_sequence_dictionary_validation: Disable sequence dictionary validation.
    - dont_use_dragstr_pair_hmm_scores: Disable DRAGstr pair-hmm scores.
    - dont_use_soft_clipped_bases: Do not analyze soft clipped bases.
    - dragen_mode: Enable DRAGEN-GATK mode.
    - dragstr_het_hom_ratio: Het to hom prior ratio for DRAGstr.
    - dragstr_params_path: Path to DRAGstr model parameters.
    - enable_dynamic_read_disqualification_for_genotyping: Enable less strict read disqualification.
    - flow_order_for_annotations: Flow order for annotations.
    - founder_id: Samples representing population founders.
    - gcs_max_retries: Max retries for GCS bucket errors.
    - gcs_project_for_requester_pays: Project for requester pays buckets.
    - genotype_assignment_method: Method to assign genotypes.
    - graph_output: Output file for debug assembly graph.
    - help: Show help message.
    - heterozygosity: Heterozygosity prior.
    - heterozygosity_stdev: Std dev of heterozygosity.
    - indel_heterozygosity: Indel heterozygosity prior.
    - interval_merging_rule: Interval merging rule.
    - intervals: List of genomic intervals to operate on.
    - interval_exclusion_padding: Padding for excluded intervals.
    - interval_padding: Padding for included intervals.
    - interval_set_rule: Interval set merging approach.
    - inverted_read_filter: Inverted read filters.
    - keep_boundary_flows: Prevent spreading of boundary flows.
    - kmer_size: List of kmer sizes for assembler.
    - lenient: Lenient VCF processing.
    - likelihood_calculation_engine: Engine for likelihood calculation.
    - linked_de_bruijn_graph: Use linked De Bruijn graph.
    - mapping_quality_threshold_for_genotyping: Mapping quality threshold for genotyping.
    - max_alternate_alleles: Max alternate alleles to genotype.
    - max_assembly_region_size: Max assembly region size.
    - max_effective_depth_adjustment_for_frd: Max depth for FRD adjustment.
    - max_genotype_count: Max genotypes to consider.
    - max_mnp_distance: Max distance to merge phased substitutions into MNPs.
    - max_num_haplotypes_in_population: Max haplotypes to consider.
    - max_prob_propagation_distance: Max bases for probability propagation.
    - max_reads_per_alignment_start: Max reads per alignment start.
    - max_unpruned_variants: Max variants allowed in graph pruner.
    - max_variants_per_shard: Max records per VCF shard.
    - min_assembly_region_size: Min assembly region size.
    - min_base_quality_score: Min base quality for calling.
    - min_dangling_branch_length: Min dangling branch length to recover.
    - min_pruning: Min support to not prune paths.
    - native_pair_hmm_threads: Threads for native PairHMM.
    - native_pair_hmm_use_double_precision: Use double precision in PairHMM.
    - num_pruning_samples: Number of samples for pruning threshold.
    - num_reference_samples_if_no_call: Number of hom-ref genotypes to infer if no call.
    - output_mode: Output call mode.
    - pedigree: Pedigree file path.
    - ploidy_regions: Interval file specifying ploidy regions.
    - population_callset: Callset for genotype priors.
    - pcr_indel_model: PCR indel error model.
    - phred_scaled_global_read_mismapping_rate: Global read mismapping rate.
    - pileup_detection: Enable pileup-based haplotypes.
    - ploidy: Sample ploidy.
    - read_filter: Read filters to apply.
    - read_index: Indices for read inputs.
    - read_validation_stringency: Validation stringency for reads.
    - recover_all_dangling_branches: Recover all dangling branches.
    - reference_model_deletion_quality: Deletion quality in reference model.
    - sample_name: Single sample name from multi-sample BAM.
    - seconds_between_progress_updates: Interval for progress updates.
    - sequence_dictionary: Sequence dictionary file.
    - show_hidden: Show hidden arguments.
    - sites_only_vcf_output: Emit VCF without genotype fields.
    - smith_waterman: Smith-Waterman implementation.
    - smith_waterman_dangling_end_gap_extend_penalty: Gap extend penalty for dangling ends.
    - smith_waterman_dangling_end_gap_open_penalty: Gap open penalty for dangling ends.
    - smith_waterman_dangling_end_match_value: Match value for dangling ends.
    - smith_waterman_dangling_end_mismatch_penalty: Mismatch penalty for dangling ends.
    - smith_waterman_haplotype_to_reference_gap_extend_penalty: Gap extend penalty haplotype to reference.
    - smith_waterman_haplotype_to_reference_gap_open_penalty: Gap open penalty haplotype to reference.
    - smith_waterman_haplotype_to_reference_match_value: Match value haplotype to reference.
    - smith_waterman_haplotype_to_reference_mismatch_penalty: Mismatch penalty haplotype to reference.
    - smith_waterman_read_to_haplotype_gap_extend_penalty: Gap extend penalty read to haplotype.
    - smith_waterman_read_to_haplotype_gap_open_penalty: Gap open penalty read to haplotype.
    - smith_waterman_read_to_haplotype_match_value: Match value read to haplotype.
    - smith_waterman_read_to_haplotype_mismatch_penalty: Mismatch penalty read to haplotype.
    - soft_clip_low_quality_ends: Preserve low-quality read ends as softclips.
    - standard_min_confidence_threshold_for_calling: Minimum confidence threshold for calling.
    - tmp_dir: Temporary directory.
    - transform_dragen_mapping_quality: Transform DRAGEN mapping quality.
    - use_filtered_reads_for_annotations: Use contamination-filtered reads for annotations.
    - use_jdk_deflater: Use JDK deflater.
    - use_jdk_inflater: Use JDK inflater.
    - use_new_qual_calculator: Use new quality calculator (deprecated).
    - use_pdhmm: Use Partially Determined HMM.
    - use_pdhmm_overlap_optimization: PDHMM overlap optimization.
    - use_posteriors_to_calculate_qual: Use genotype posteriors to calculate QUAL.
    - verbosity: Logging verbosity level.
    - version: Show version information.
    """
    import shlex

    # Validate input files
    if not input_bams or len(input_bams) == 0:
        raise ValueError("At least one input BAM file must be provided.")
    for bam in input_bams:
        if not bam.exists():
            raise FileNotFoundError(f"Input BAM file not found: {bam}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference file not found: {reference}")
    if alleles is not None and not alleles.exists():
        raise FileNotFoundError(f"Alleles file not found: {alleles}")
    if dbsnp is not None and not dbsnp.exists():
        raise FileNotFoundError(f"dbSNP file not found: {dbsnp}")
    if pedigree is not None and not pedigree.exists():
        raise FileNotFoundError(f"Pedigree file not found: {pedigree}")
    if dragstr_params_path is not None and not dragstr_params_path.exists():
        raise FileNotFoundError(f"DRAGstr params file not found: {dragstr_params_path}")
    if contamination_fraction_to_filter < 0.0:
        raise ValueError("contamination_fraction_to_filter must be >= 0.0")
    if base_quality_score_threshold < 0 or base_quality_score_threshold > 127:
        raise ValueError("base_quality_score_threshold must be between 0 and 127")
    if min_base_quality_score < 0 or min_base_quality_score > 127:
        raise ValueError("min_base_quality_score must be between 0 and 127")
    if max_alternate_alleles < 1:
        raise ValueError("max_alternate_alleles must be >= 1")
    if max_genotype_count < 1:
        raise ValueError("max_genotype_count must be >= 1")
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")
    if standard_min_confidence_threshold_for_calling < 0.0:
        raise ValueError("standard_min_confidence_threshold_for_calling must be >= 0.0")
    if native_pair_hmm_threads < 1:
        raise ValueError("native_pair_hmm_threads must be >= 1")
    if gcs_max_retries < 0:
        raise ValueError("gcs_max_retries must be >= 0")
    if max_variants_per_shard < 0:
        raise ValueError("max_variants_per_shard must be >= 0")

    # Prepare command
    cmd = ["gatk", "--java-options", "-Xmx4g", "HaplotypeCaller"]

    # Required arguments
    cmd += ["-I"] + [str(bam) for bam in input_bams]
    cmd += ["-O", str(output)]
    cmd += ["-R", str(reference)]

    # Optional arguments
    if alleles:
        cmd += ["--alleles", str(alleles)]
    if annotate_with_num_discovered_alleles:
        cmd.append("--annotate-with-num-discovered-alleles")
    if annotation:
        for ann in annotation:
            cmd += ["-A", ann]
    if annotation_group:
        for ag in annotation_group:
            cmd += ["-G", ag]
    if annotations_to_exclude:
        for ax in annotations_to_exclude:
            cmd += ["-AX", ax]
    if arguments_file:
        for af in arguments_file:
            if not af.exists():
                raise FileNotFoundError(f"Arguments file not found: {af}")
            cmd += ["--arguments_file", str(af)]
    if assembly_region_out:
        cmd += ["--assembly-region-out", str(assembly_region_out)]
    if assembly_region_padding != 100:
        cmd += ["--assembly-region-padding", str(assembly_region_padding)]
    if base_quality_score_threshold != 18:
        cmd += ["--base-quality-score-threshold", str(base_quality_score_threshold)]
    if cloud_index_prefetch_buffer != -1:
        cmd += ["--cloud-index-prefetch-buffer", str(cloud_index_prefetch_buffer)]
    if cloud_prefetch_buffer != 40:
        cmd += ["--cloud-prefetch-buffer", str(cloud_prefetch_buffer)]
    if contamination_fraction_to_filter != 0.0:
        cmd += ["--contamination-fraction-to-filter", str(contamination_fraction_to_filter)]
    if dbsnp:
        cmd += ["--dbsnp", str(dbsnp)]
    if disable_bam_index_caching:
        cmd.append("--disable-bam-index-caching")
    if disable_sequence_dictionary_validation:
        cmd.append("--disable-sequence-dictionary-validation")
    if dont_use_dragstr_pair_hmm_scores:
        cmd.append("--dont-use-dragstr-pair-hmm-scores")
    if dont_use_soft_clipped_bases:
        cmd.append("--dont-use-soft-clipped-bases")
    if dragen_mode:
        cmd.append("--dragen-mode")
    if dragstr_het_hom_ratio != 2:
        cmd += ["--dragstr-het-hom-ratio", str(dragstr_het_hom_ratio)]
    if dragstr_params_path:
        cmd += ["--dragstr-params-path", str(dragstr_params_path)]
    if enable_dynamic_read_disqualification_for_genotyping:
        cmd.append("--enable-dynamic-read-disqualification-for-genotyping")
    if flow_order_for_annotations:
        for fo in flow_order_for_annotations:
            cmd += ["--flow-order-for-annotations", fo]
    if founder_id:
        for fid in founder_id:
            cmd += ["--founder-id", fid]
    if gcs_max_retries != 20:
        cmd += ["--gcs-max-retries", str(gcs_max_retries)]
    if gcs_project_for_requester_pays:
        cmd += ["--gcs-project-for-requester-pays", gcs_project_for_requester_pays]
    if genotype_assignment_method != "USE_PLS_TO_ASSIGN":
        cmd += ["--genotype-assignment-method", genotype_assignment_method]
    if graph_output:
        cmd += ["--graph-output", str(graph_output)]
    if help:
        cmd.append("--help")
    if heterozygosity != 0.001:
        cmd += ["--heterozygosity", str(heterozygosity)]
    if heterozygosity_stdev != 0.01:
        cmd += ["--heterozygosity-stdev", str(heterozygosity_stdev)]
    if indel_heterozygosity != 0.000125:
        cmd += ["--indel-heterozygosity", str(indel_heterozygosity)]
    if interval_merging_rule != "ALL":
        cmd += ["--interval-merging-rule", interval_merging_rule]
    if intervals:
        for interval in intervals:
            cmd += ["-L", interval]
    if interval_exclusion_padding != 0:
        cmd += ["--interval-exclusion-padding", str(interval_exclusion_padding)]
    if interval_padding != 0:
        cmd += ["--interval-padding", str(interval_padding)]
    if interval_set_rule != "UNION":
        cmd += ["--interval-set-rule", interval_set_rule]
    if inverted_read_filter:
        for irf in inverted_read_filter:
            cmd += ["--inverted-read-filter", irf]
    if keep_boundary_flows:
        cmd.append("--keep-boundary-flows")
    if kmer_size:
        for k in kmer_size:
            cmd += ["--kmer-size", str(k)]
    if lenient:
        cmd.append("--lenient")
    if likelihood_calculation_engine != "PairHMM":
        cmd += ["--likelihood-calculation-engine", likelihood_calculation_engine]
    if linked_de_bruijn_graph:
        cmd.append("--linked-de-bruijn-graph")
    if mapping_quality_threshold_for_genotyping != 20:
        cmd += ["--mapping-quality-threshold-for-genotyping", str(mapping_quality_threshold_for_genotyping)]
    if max_alternate_alleles != 6:
        cmd += ["--max-alternate-alleles", str(max_alternate_alleles)]
    if max_assembly_region_size != 300:
        cmd += ["--max-assembly-region-size", str(max_assembly_region_size)]
    if max_effective_depth_adjustment_for_frd != 0:
        cmd += ["--max-effective-depth-adjustment-for-frd", str(max_effective_depth_adjustment_for_frd)]
    if max_genotype_count != 1024:
        cmd += ["--max-genotype-count", str(max_genotype_count)]
    if max_mnp_distance != 0:
        cmd += ["--max-mnp-distance", str(max_mnp_distance)]
    if max_num_haplotypes_in_population != 128:
        cmd += ["--max-num-haplotypes-in-population", str(max_num_haplotypes_in_population)]
    if max_prob_propagation_distance != 50:
        cmd += ["--max-prob-propagation-distance", str(max_prob_propagation_distance)]
    if max_reads_per_alignment_start != 50:
        cmd += ["--max-reads-per-alignment-start", str(max_reads_per_alignment_start)]
    if max_unpruned_variants != 100:
        cmd += ["--max-unpruned-variants", str(max_unpruned_variants)]
    if max_variants_per_shard != 0:
        cmd += ["--max-variants-per-shard", str(max_variants_per_shard)]
    if min_assembly_region_size != 50:
        cmd += ["--min-assembly-region-size", str(min_assembly_region_size)]
    if min_base_quality_score != 10:
        cmd += ["--min-base-quality-score", str(min_base_quality_score)]
    if min_dangling_branch_length != 4:
        cmd += ["--min-dangling-branch-length", str(min_dangling_branch_length)]
    if min_pruning != 2:
        cmd += ["--min-pruning", str(min_pruning)]
    if native_pair_hmm_threads != 4:
        cmd += ["--native-pair-hmm-threads", str(native_pair_hmm_threads)]
    if native_pair_hmm_use_double_precision:
        cmd.append("--native-pair-hmm-use-double-precision")
    if num_pruning_samples != 1:
        cmd += ["--num-pruning-samples", str(num_pruning_samples)]
    if num_reference_samples_if_no_call != 0:
        cmd += ["--num-reference-samples-if-no-call", str(num_reference_samples_if_no_call)]
    if output_mode != "EMIT_VARIANTS_ONLY":
        cmd += ["--output-mode", output_mode]
    if pedigree:
        cmd += ["--pedigree", str(pedigree)]
    if ploidy_regions:
        cmd += ["--ploidy-regions", str(ploidy_regions)]
    if population_callset:
        cmd += ["--population-callset", str(population_callset)]
    if pcr_indel_model != "CONSERVATIVE":
        cmd += ["--pcr-indel-model", pcr_indel_model]
    if phred_scaled_global_read_mismapping_rate != 45:
        cmd += ["--phred-scaled-global-read-mismapping-rate", str(phred_scaled_global_read_mismapping_rate)]
    if pileup_detection:
        cmd.append("--pileup-detection")
    if ploidy != 2:
        cmd += ["--sample-ploidy", str(ploidy)]
    if read_filter:
        for rf in read_filter:
            cmd += ["--read-filter", rf]
    if read_index:
        for ri in read_index:
            cmd += ["--read-index", str(ri)]
    if read_validation_stringency != "SILENT":
        cmd += ["--read-validation-stringency", read_validation_stringency]
    if recover_all_dangling_branches:
        cmd.append("--recover-all-dangling-branches")
    if reference_model_deletion_quality != 30:
        cmd += ["--reference-model-deletion-quality", str(reference_model_deletion_quality)]
    if sample_name:
        cmd += ["--sample-name", sample_name]
    if seconds_between_progress_updates != 10.0:
        cmd += ["--seconds-between-progress-updates", str(seconds_between_progress_updates)]
    if sequence_dictionary:
        cmd += ["--sequence-dictionary", str(sequence_dictionary)]
    if show_hidden:
        cmd.append("--showHidden")
    if sites_only_vcf_output:
        cmd.append("--sites-only-vcf-output")
    if smith_waterman != "FASTEST_AVAILABLE":
        cmd += ["--smith-waterman", smith_waterman]
    if smith_waterman_dangling_end_gap_extend_penalty != -6:
        cmd += ["--smith-waterman-dangling-end-gap-extend-penalty", str(smith_waterman_dangling_end_gap_extend_penalty)]
    if smith_waterman_dangling_end_gap_open_penalty != -110:
        cmd += ["--smith-waterman-dangling-end-gap-open-penalty", str(smith_waterman_dangling_end_gap_open_penalty)]
    if smith_waterman_dangling_end_match_value != 25:
        cmd += ["--smith-waterman-dangling-end-match-value", str(smith_waterman_dangling_end_match_value)]
    if smith_waterman_dangling_end_mismatch_penalty != -50:
        cmd += ["--smith-waterman-dangling-end-mismatch-penalty", str(smith_waterman_dangling_end_mismatch_penalty)]
    if smith_waterman_haplotype_to_reference_gap_extend_penalty != -11:
        cmd += ["--smith-waterman-haplotype-to-reference-gap-extend-penalty", str(smith_waterman_haplotype_to_reference_gap_extend_penalty)]
    if smith_waterman_haplotype_to_reference_gap_open_penalty != -260:
        cmd += ["--smith-waterman-haplotype-to-reference-gap-open-penalty", str(smith_waterman_haplotype_to_reference_gap_open_penalty)]
    if smith_waterman_haplotype_to_reference_match_value != 200:
        cmd += ["--smith-waterman-haplotype-to-reference-match-value", str(smith_waterman_haplotype_to_reference_match_value)]
    if smith_waterman_haplotype_to_reference_mismatch_penalty != -150:
        cmd += ["--smith-waterman-haplotype-to-reference-mismatch-penalty", str(smith_waterman_haplotype_to_reference_mismatch_penalty)]
    if smith_waterman_read_to_haplotype_gap_extend_penalty != -5:
        cmd += ["--smith-waterman-read-to-haplotype-gap-extend-penalty", str(smith_waterman_read_to_haplotype_gap_extend_penalty)]
    if smith_waterman_read_to_haplotype_gap_open_penalty != -30:
        cmd += ["--smith-waterman-read-to-haplotype-gap-open-penalty", str(smith_waterman_read_to_haplotype_gap_open_penalty)]
    if smith_waterman_read_to_haplotype_match_value != 10:
        cmd += ["--smith-waterman-read-to-haplotype-match-value", str(smith_waterman_read_to_haplotype_match_value)]
    if smith_waterman_read_to_haplotype_mismatch_penalty != -15:
        cmd += ["--smith-waterman-read-to-haplotype-mismatch-penalty", str(smith_waterman_read_to_haplotype_mismatch_penalty)]
    if soft_clip_low_quality_ends:
        cmd.append("--soft-clip-low-quality-ends")
    if standard_min_confidence_threshold_for_calling != 30.0:
        cmd += ["--standard-min-confidence-threshold-for-calling", str(standard_min_confidence_threshold_for_calling)]
    if tmp_dir:
        if not tmp_dir.exists():
            raise FileNotFoundError(f"Temporary directory does not exist: {tmp_dir}")
        cmd += ["--tmp-dir", str(tmp_dir)]
    if transform_dragen_mapping_quality:
        cmd.append("--transform-dragen-mapping-quality")
    if use_filtered_reads_for_annotations:
        cmd.append("--use-filtered-reads-for-annotations")
    if use_jdk_deflater:
        cmd.append("--use-jdk-deflater")
    if use_jdk_inflater:
        cmd.append("--use-jdk-inflater")
    if not use_new_qual_calculator:
        cmd.append("--use-new-qual-calculator")
    if use_pdhmm:
        cmd.append("--use-pdhmm")
    if use_pdhmm_overlap_optimization:
        cmd.append("--use-pdhmm-overlap-optimization")
    if use_posteriors_to_calculate_qual:
        cmd.append("--use-posteriors-to-calculate-qual")
    if verbosity != "INFO":
        cmd += ["--verbosity", verbosity]
    if version:
        cmd.append("--version")

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
    if assembly_region_out:
        output_files.append(str(assembly_region_out))
    if graph_output:
        output_files.append(str(graph_output))
    if tmp_dir:
        output_files.append(str(tmp_dir))

    return {
        "command_executed": " ".join(shlex.quote(c) for c in cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()