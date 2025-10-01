from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional, Tuple
import subprocess

mcp = FastMCP()


@mcp.tool()
def freebayes(
    fasta_reference: Path,
    bam_files: List[Path] = [],
    bam_list: Optional[Path] = None,
    stdin: bool = False,
    targets: Optional[Path] = None,
    region: Optional[str] = None,
    samples: Optional[Path] = None,
    populations: Optional[Path] = None,
    cnv_map: Optional[Path] = None,
    vcf_output: Optional[Path] = None,
    gvcf: bool = False,
    gvcf_chunk: Optional[int] = None,
    gvcf_dont_use_chunk: Optional[bool] = None,
    variant_input: Optional[Path] = None,
    only_use_input_alleles: bool = False,
    haplotype_basis_alleles: Optional[Path] = None,
    report_all_haplotype_alleles: bool = False,
    report_monomorphic: bool = False,
    pvar: float = 0.0,
    strict_vcf: bool = False,
    theta: float = 0.001,
    ploidy: int = 2,
    pooled_discrete: bool = False,
    pooled_continuous: bool = False,
    use_reference_allele: bool = False,
    reference_quality: Optional[str] = None,  # format "MQ,BQ"
    use_best_n_alleles: int = 0,
    max_complex_gap: int = 3,
    haplotype_length: Optional[int] = None,
    min_repeat_size: int = 5,
    min_repeat_entropy: float = 1.0,
    no_partial_observations: bool = False,
    throw_away_snp_obs: bool = False,
    throw_away_indels_obs: bool = False,
    throw_away_mnp_obs: bool = False,
    throw_away_complex_obs: bool = False,
    dont_left_align_indels: bool = False,
    use_duplicate_reads: bool = False,
    min_mapping_quality: int = 1,
    min_base_quality: int = 0,
    min_supporting_allele_qsum: int = 0,
    min_supporting_mapping_qsum: int = 0,
    mismatch_base_quality_threshold: int = 10,
    read_mismatch_limit: Optional[int] = None,
    read_max_mismatch_fraction: float = 1.0,
    read_snp_limit: Optional[int] = None,
    read_indel_limit: Optional[int] = None,
    standard_filters: bool = False,
    min_alternate_fraction: float = 0.05,
    min_alternate_count: int = 2,
    min_alternate_qsum: int = 0,
    min_alternate_total: int = 1,
    min_coverage: int = 0,
    limit_coverage: Optional[int] = None,
    skip_coverage: Optional[int] = None,
    trim_complex_tail: bool = False,
    no_population_priors: bool = False,
    hwe_priors_or: bool = False,
    binomial_obs_priors_or: bool = False,
    allele_balance_priors_or: bool = False,
    observation_bias: Optional[Path] = None,
    base_quality_cap: Optional[int] = None,
    prob_contamination: float = 1e-8,
    legacy_gls: bool = False,
    contamination_estimates: Optional[Path] = None,
    report_genotype_likelihood_max: bool = False,
    genotyping_max_iterations: int = 1000,
    genotyping_max_banddepth: int = 6,
    posterior_integration_limits: Optional[Tuple[int, int]] = None,
    exclude_unobserved_genotypes: bool = False,
    genotype_variant_threshold: Optional[float] = None,
    use_mapping_quality: bool = False,
    harmonic_indel_quality: bool = False,
    read_dependence_factor: float = 0.9,
    genotype_qualities: bool = False,
    debug: bool = False,
    debug_verbose: bool = False,
) -> dict:
    """
    Run FreeBayes Bayesian haplotype-based polymorphism discovery on BAM files with a reference.

    Parameters:
    - fasta_reference: Reference FASTA file (required).
    - bam_files: List of BAM files to analyze.
    - bam_list: File containing list of BAM files.
    - stdin: Read BAM input from stdin.
    - targets: BED file to limit analysis to targets.
    - region: Region string <chrom>:<start>-<end> to limit analysis.
    - samples: File listing samples to analyze.
    - populations: File listing sample-population pairs.
    - cnv_map: Copy number variation map BED file.
    - vcf_output: Output VCF file path (default stdout).
    - gvcf: Write gVCF output.
    - gvcf_chunk: Emit gVCF record every NUM bases.
    - gvcf_dont_use_chunk: Emit gVCF record for all bases if true.
    - variant_input: Input VCF file with variants.
    - only_use_input_alleles: Only call alleles in input VCF.
    - haplotype_basis_alleles: VCF file for haplotype basis alleles.
    - report_all_haplotype_alleles: Report info about all haplotype alleles.
    - report_monomorphic: Report monomorphic loci.
    - pvar: Minimum polymorphism probability to report.
    - strict_vcf: Generate strict VCF format.
    - theta: Expected mutation rate (default 0.001).
    - ploidy: Default ploidy (default 2).
    - pooled_discrete: Model pooled samples with discrete genotypes.
    - pooled_continuous: Frequency-based pooled caller.
    - use_reference_allele: Include reference allele in analysis.
    - reference_quality: Mapping and base quality for reference allele as "MQ,BQ".
    - use_best_n_alleles: Evaluate only best N SNP alleles (0=all).
    - max_complex_gap: Max gap for haplotype calls (default 3).
    - haplotype_length: Haplotype length for clumping.
    - min_repeat_size: Minimum repeat size (default 5).
    - min_repeat_entropy: Minimum repeat entropy (default 1.0).
    - no_partial_observations: Exclude partial observations.
    - throw_away_snp_obs: Remove SNP observations.
    - throw_away_indels_obs: Remove indel observations.
    - throw_away_mnp_obs: Remove MNP observations.
    - throw_away_complex_obs: Remove complex allele observations.
    - dont_left_align_indels: Disable left-alignment of indels.
    - use_duplicate_reads: Include duplicate-marked alignments.
    - min_mapping_quality: Minimum mapping quality (default 1).
    - min_base_quality: Minimum base quality (default 0).
    - min_supporting_allele_qsum: Minimum sum of allele qualities (default 0).
    - min_supporting_mapping_qsum: Minimum sum of mapping qualities (default 0).
    - mismatch_base_quality_threshold: Base quality threshold for mismatches (default 10).
    - read_mismatch_limit: Max mismatches per read (None=unbounded).
    - read_max_mismatch_fraction: Max mismatch fraction per read (default 1.0).
    - read_snp_limit: Max SNP mismatches per read (None=unbounded).
    - read_indel_limit: Max indels per read (None=unbounded).
    - standard_filters: Use stringent filters (-m30 -q20 -R0 -S0).
    - min_alternate_fraction: Minimum fraction of alt observations (default 0.05).
    - min_alternate_count: Minimum count of alt observations (default 2).
    - min_alternate_qsum: Minimum quality sum of alt observations (default 0).
    - min_alternate_total: Minimum alt observations in population (default 1).
    - min_coverage: Minimum coverage to process site (default 0).
    - limit_coverage: Downsample coverage limit (None=no limit).
    - skip_coverage: Skip sites with coverage > N (None=no limit).
    - trim_complex_tail: Trim complex tails.
    - no_population_priors: Disable population priors.
    - hwe_priors_or: Disable HWE priors.
    - binomial_obs_priors_or: Disable binomial observation priors.
    - allele_balance_priors_or: Disable allele balance priors.
    - observation_bias: File with allele observation biases.
    - base_quality_cap: Cap base quality.
    - prob_contamination: Contamination estimate (default 1e-8).
    - legacy_gls: Use legacy genotype likelihoods.
    - contamination_estimates: File with per-sample contamination estimates.
    - report_genotype_likelihood_max: Report max likelihood genotypes.
    - genotyping_max_iterations: Max genotyping iterations (default 1000).
    - genotyping_max_banddepth: Max genotype banddepth (default 6).
    - posterior_integration_limits: Tuple (N,M) for posterior integration limits.
    - exclude_unobserved_genotypes: Skip genotyping unobserved genotypes.
    - genotype_variant_threshold: Limit posterior integration threshold.
    - use_mapping_quality: Use mapping quality in likelihoods.
    - harmonic_indel_quality: Use harmonic indel quality.
    - read_dependence_factor: Read dependence factor (default 0.9).
    - genotype_qualities: Calculate genotype qualities.
    - debug: Print debugging output.
    - debug_verbose: Print verbose debugging output.

    Returns:
        dict: command_executed, stdout, stderr, output_files (VCF output if specified)
    """
    # Validate paths
    if not fasta_reference.exists():
        raise FileNotFoundError(f"Reference FASTA file not found: {fasta_reference}")
    if bam_list is not None and not bam_list.exists():
        raise FileNotFoundError(f"BAM list file not found: {bam_list}")
    for bam in bam_files:
        if not bam.exists():
            raise FileNotFoundError(f"BAM file not found: {bam}")
    if targets is not None and not targets.exists():
        raise FileNotFoundError(f"Targets BED file not found: {targets}")
    if samples is not None and not samples.exists():
        raise FileNotFoundError(f"Samples file not found: {samples}")
    if populations is not None and not populations.exists():
        raise FileNotFoundError(f"Populations file not found: {populations}")
    if cnv_map is not None and not cnv_map.exists():
        raise FileNotFoundError(f"CNV map file not found: {cnv_map}")
    if variant_input is not None and not variant_input.exists():
        raise FileNotFoundError(f"Variant input VCF file not found: {variant_input}")
    if haplotype_basis_alleles is not None and not haplotype_basis_alleles.exists():
        raise FileNotFoundError(f"Haplotype basis alleles VCF file not found: {haplotype_basis_alleles}")
    if observation_bias is not None and not observation_bias.exists():
        raise FileNotFoundError(f"Observation bias file not found: {observation_bias}")
    if contamination_estimates is not None and not contamination_estimates.exists():
        raise FileNotFoundError(f"Contamination estimates file not found: {contamination_estimates}")

    # Validate numeric parameters
    if pvar < 0.0 or pvar > 1.0:
        raise ValueError("pvar must be between 0.0 and 1.0")
    if theta < 0.0:
        raise ValueError("theta must be non-negative")
    if ploidy < 1:
        raise ValueError("ploidy must be at least 1")
    if use_best_n_alleles < 0:
        raise ValueError("use_best_n_alleles must be >= 0")
    if max_complex_gap < -1:
        raise ValueError("max_complex_gap must be >= -1")
    if min_repeat_size < 0:
        raise ValueError("min_repeat_size must be >= 0")
    if min_repeat_entropy < 0.0:
        raise ValueError("min_repeat_entropy must be >= 0.0")
    if min_mapping_quality < 0:
        raise ValueError("min_mapping_quality must be >= 0")
    if min_base_quality < 0:
        raise ValueError("min_base_quality must be >= 0")
    if min_supporting_allele_qsum < 0:
        raise ValueError("min_supporting_allele_qsum must be >= 0")
    if min_supporting_mapping_qsum < 0:
        raise ValueError("min_supporting_mapping_qsum must be >= 0")
    if mismatch_base_quality_threshold < 0:
        raise ValueError("mismatch_base_quality_threshold must be >= 0")
    if read_mismatch_limit is not None and read_mismatch_limit < 0:
        raise ValueError("read_mismatch_limit must be >= 0")
    if not (0.0 <= read_max_mismatch_fraction <= 1.0):
        raise ValueError("read_max_mismatch_fraction must be between 0.0 and 1.0")
    if read_snp_limit is not None and read_snp_limit < 0:
        raise ValueError("read_snp_limit must be >= 0")
    if read_indel_limit is not None and read_indel_limit < 0:
        raise ValueError("read_indel_limit must be >= 0")
    if min_alternate_fraction < 0.0 or min_alternate_fraction > 1.0:
        raise ValueError("min_alternate_fraction must be between 0.0 and 1.0")
    if min_alternate_count < 0:
        raise ValueError("min_alternate_count must be >= 0")
    if min_alternate_qsum < 0:
        raise ValueError("min_alternate_qsum must be >= 0")
    if min_alternate_total < 0:
        raise ValueError("min_alternate_total must be >= 0")
    if min_coverage < 0:
        raise ValueError("min_coverage must be >= 0")
    if limit_coverage is not None and limit_coverage < 0:
        raise ValueError("limit_coverage must be >= 0")
    if skip_coverage is not None and skip_coverage < 0:
        raise ValueError("skip_coverage must be >= 0")
    if base_quality_cap is not None and base_quality_cap < 0:
        raise ValueError("base_quality_cap must be >= 0")
    if prob_contamination < 0.0 or prob_contamination > 1.0:
        raise ValueError("prob_contamination must be between 0.0 and 1.0")
    if genotyping_max_iterations < 1:
        raise ValueError("genotyping_max_iterations must be >= 1")
    if genotyping_max_banddepth < 1:
        raise ValueError("genotyping_max_banddepth must be >= 1")
    if posterior_integration_limits is not None:
        if len(posterior_integration_limits) != 2:
            raise ValueError("posterior_integration_limits must be a tuple of two integers")
        if posterior_integration_limits[0] < 0 or posterior_integration_limits[1] < 0:
            raise ValueError("posterior_integration_limits values must be >= 0")
    if genotype_variant_threshold is not None and genotype_variant_threshold <= 0:
        raise ValueError("genotype_variant_threshold must be > 0")
    if read_dependence_factor < 0.0 or read_dependence_factor > 1.0:
        raise ValueError("read_dependence_factor must be between 0.0 and 1.0")

    # Build command line
    cmd = ["freebayes"]

    # Required reference
    cmd += ["-f", str(fasta_reference)]

    # BAM inputs
    if stdin:
        cmd.append("-c")
    if bam_list:
        cmd += ["-L", str(bam_list)]
    if bam_files:
        for bam in bam_files:
            cmd += ["-b", str(bam)]

    # Targets and regions
    if targets:
        cmd += ["-t", str(targets)]
    if region:
        cmd += ["-r", region]

    # Samples and populations
    if samples:
        cmd += ["-s", str(samples)]
    if populations:
        cmd += ["--populations", str(populations)]

    # CNV map
    if cnv_map:
        cmd += ["-A", str(cnv_map)]

    # Output
    if vcf_output:
        cmd += ["-v", str(vcf_output)]
    if gvcf:
        cmd.append("--gvcf")
    if gvcf_chunk is not None:
        if gvcf_chunk < 1:
            raise ValueError("gvcf_chunk must be >= 1")
        cmd += ["--gvcf-chunk", str(gvcf_chunk)]
    if gvcf_dont_use_chunk is not None:
        cmd += ["-&", "true" if gvcf_dont_use_chunk else "false"]

    # Variant input and allele options
    if variant_input:
        cmd += ["-@", str(variant_input)]
    if only_use_input_alleles:
        cmd.append("-l")
    if haplotype_basis_alleles:
        cmd += ["--haplotype-basis-alleles", str(haplotype_basis_alleles)]
    if report_all_haplotype_alleles:
        cmd.append("--report-all-haplotype-alleles")
    if report_monomorphic:
        cmd.append("--report-monomorphic")
    if pvar > 0.0:
        cmd += ["-P", str(pvar)]
    if strict_vcf:
        cmd.append("--strict-vcf")

    # Population model
    cmd += ["-T", str(theta)]
    cmd += ["-p", str(ploidy)]
    if pooled_discrete:
        cmd.append("-J")
    if pooled_continuous:
        cmd.append("-K")

    # Reference allele
    if use_reference_allele:
        cmd.append("-Z")
    if reference_quality:
        # Validate format MQ,BQ
        parts = reference_quality.split(",")
        if len(parts) != 2:
            raise ValueError("reference_quality must be in format MQ,BQ")
        mq, bq = parts
        if not mq.isdigit() or not bq.isdigit():
            raise ValueError("reference_quality MQ and BQ must be integers")
        cmd += ["--reference-quality", reference_quality]

    # Allele scope
    if use_best_n_alleles > 0:
        cmd += ["-n", str(use_best_n_alleles)]
    if max_complex_gap != 3:
        cmd += ["-E", str(max_complex_gap)]
    if haplotype_length is not None:
        cmd += ["--haplotype-length", str(haplotype_length)]
    if min_repeat_size != 5:
        cmd += ["--min-repeat-size", str(min_repeat_size)]
    if min_repeat_entropy != 1.0:
        cmd += ["--min-repeat-entropy", str(min_repeat_entropy)]
    if no_partial_observations:
        cmd.append("--no-partial-observations")

    # Throw away observations
    if throw_away_snp_obs:
        cmd.append("-I")
    if throw_away_indels_obs:
        cmd.append("-i")
    if throw_away_mnp_obs:
        cmd.append("-X")
    if throw_away_complex_obs:
        cmd.append("-u")

    # Indel realignment
    if dont_left_align_indels:
        cmd.append("-O")

    # Input filters
    if use_duplicate_reads:
        cmd.append("-4")
    if min_mapping_quality != 1:
        cmd += ["-m", str(min_mapping_quality)]
    if min_base_quality != 0:
        cmd += ["-q", str(min_base_quality)]
    if min_supporting_allele_qsum != 0:
        cmd += ["-R", str(min_supporting_allele_qsum)]
    if min_supporting_mapping_qsum != 0:
        cmd += ["-Y", str(min_supporting_mapping_qsum)]
    if mismatch_base_quality_threshold != 10:
        cmd += ["-Q", str(mismatch_base_quality_threshold)]
    if read_mismatch_limit is not None:
        cmd += ["-U", str(read_mismatch_limit)]
    if read_max_mismatch_fraction != 1.0:
        cmd += ["-z", str(read_max_mismatch_fraction)]
    if read_snp_limit is not None:
        cmd += ["-$", str(read_snp_limit)]
    if read_indel_limit is not None:
        cmd += ["-e", str(read_indel_limit)]
    if standard_filters:
        cmd.append("-0")
    if min_alternate_fraction != 0.05:
        cmd += ["-F", str(min_alternate_fraction)]
    if min_alternate_count != 2:
        cmd += ["-C", str(min_alternate_count)]
    if min_alternate_qsum != 0:
        cmd += ["-3", str(min_alternate_qsum)]
    if min_alternate_total != 1:
        cmd += ["-G", str(min_alternate_total)]
    if min_coverage != 0:
        cmd += ["--min-coverage", str(min_coverage)]
    if limit_coverage is not None:
        cmd += ["--limit-coverage", str(limit_coverage)]
    if skip_coverage is not None:
        cmd += ["-g", str(skip_coverage)]
    if trim_complex_tail:
        cmd.append("--trim-complex-tail")

    # Population priors
    if no_population_priors:
        cmd.append("-k")

    # Mappability priors
    if hwe_priors_or:
        cmd.append("-w")
    if binomial_obs_priors_or:
        cmd.append("-V")
    if allele_balance_priors_or:
        cmd.append("-a")

    # Genotype likelihoods
    if observation_bias:
        cmd += ["--observation-bias", str(observation_bias)]
    if base_quality_cap is not None:
        cmd += ["--base-quality-cap", str(base_quality_cap)]
    if prob_contamination != 1e-8:
        cmd += ["--prob-contamination", str(prob_contamination)]
    if legacy_gls:
        cmd.append("--legacy-gls")
    if contamination_estimates:
        cmd += ["--contamination-estimates", str(contamination_estimates)]

    # Algorithmic features
    if report_genotype_likelihood_max:
        cmd.append("--report-genotype-likelihood-max")
    if genotyping_max_iterations != 1000:
        cmd += ["-B", str(genotyping_max_iterations)]
    if genotyping_max_banddepth != 6:
        cmd += ["--genotyping-max-banddepth", str(genotyping_max_banddepth)]
    if posterior_integration_limits is not None:
        cmd += ["-W", f"{posterior_integration_limits[0]},{posterior_integration_limits[1]}"]
    if exclude_unobserved_genotypes:
        cmd.append("-N")
    if genotype_variant_threshold is not None:
        cmd += ["-S", str(genotype_variant_threshold)]
    if use_mapping_quality:
        cmd.append("-j")
    if harmonic_indel_quality:
        cmd.append("-H")
    if read_dependence_factor != 0.9:
        cmd += ["-D", str(read_dependence_factor)]
    if genotype_qualities:
        cmd.append("-=")

    # Debugging
    if debug:
        cmd.append("-d")
    if debug_verbose:
        cmd.append("-dd")

    # Execute command
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"FreeBayes execution failed with return code {e.returncode}",
        }

    output_files = []
    if vcf_output:
        output_files.append(str(vcf_output))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


if __name__ == '__main__':
    mcp.run()