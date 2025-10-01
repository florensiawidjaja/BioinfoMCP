from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()

def _validate_file_path(path: str, must_exist: bool = True) -> Path:
    p = Path(path)
    if must_exist and not p.exists():
        raise FileNotFoundError(f"File not found: {path}")
    return p

def _validate_output_path(path: Optional[str]) -> Optional[Path]:
    if path is None:
        return None
    p = Path(path)
    if p.exists() and not p.is_file():
        raise ValueError(f"Output path exists and is not a file: {path}")
    return p

def _build_common_options(
    collapse: Optional[str] = None,
    apply_filters: Optional[str] = None,
    no_version: bool = False,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    targets: Optional[str] = None,
    targets_file: Optional[str] = None,
    targets_overlap: Optional[str] = None,
    threads: int = 0,
    verbosity: int = 1,
    write_index: Optional[str] = None,
) -> List[str]:
    opts = []
    if collapse:
        if collapse not in {"snps", "indels", "both", "all", "some", "none", "id"}:
            raise ValueError(f"Invalid collapse value: {collapse}")
        opts += ["-c", collapse]
    if apply_filters:
        opts += ["-f", apply_filters]
    if no_version:
        opts.append("--no-version")
    if output:
        opts += ["-o", output]
    if output_type:
        if not output_type[0] in {"b", "u", "z", "v"}:
            raise ValueError(f"Invalid output-type value: {output_type}")
        opts += ["-O", output_type]
    if regions:
        opts += ["-r", regions]
    if regions_file:
        opts += ["-R", regions_file]
    if regions_overlap:
        if regions_overlap not in {"pos", "record", "variant", "0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        opts += ["--regions-overlap", regions_overlap]
    if samples:
        opts += ["-s", samples]
    if samples_file:
        opts += ["-S", samples_file]
    if targets:
        opts += ["-t", targets]
    if targets_file:
        opts += ["-T", targets_file]
    if targets_overlap:
        if targets_overlap not in {"pos", "record", "variant", "0", "1", "2"}:
            raise ValueError(f"Invalid targets-overlap value: {targets_overlap}")
        opts += ["--targets-overlap", targets_overlap]
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        opts += ["--threads", str(threads)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        opts += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        opts += ["-W", write_index]
    return opts

@mcp.tool()
def bcftools_annotate(
    file: str,
    annotations: Optional[str] = None,
    columns: Optional[str] = None,
    columns_file: Optional[str] = None,
    exclude: Optional[str] = None,
    force: bool = False,
    header_lines: Optional[str] = None,
    set_id: Optional[str] = None,
    include: Optional[str] = None,
    keep_sites: bool = False,
    merge_logic: Optional[str] = None,
    mark_sites: Optional[str] = None,
    min_overlap: Optional[str] = None,
    no_version: bool = False,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    pair_logic: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    rename_annots: Optional[str] = None,
    rename_chrs: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    single_overlaps: bool = False,
    threads: int = 0,
    remove: Optional[str] = None,
    verbosity: int = 1,
    write_index: Optional[str] = None,
):
    """
    Add or remove annotations in VCF/BCF files using bcftools annotate.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "annotate"]
    if annotations:
        ann_path = _validate_file_path(annotations)
        cmd += ["-a", str(ann_path)]
    if columns:
        cmd += ["-c", columns]
    if columns_file:
        cf_path = _validate_file_path(columns_file)
        cmd += ["-C", str(cf_path)]
    if exclude:
        cmd += ["-e", exclude]
    if force:
        cmd.append("--force")
    if header_lines:
        hl_path = _validate_file_path(header_lines)
        cmd += ["-h", str(hl_path)]
    if set_id:
        cmd += ["-I", set_id]
    if include:
        cmd += ["-i", include]
    if keep_sites:
        cmd.append("-k")
    if merge_logic:
        cmd += ["-l", merge_logic]
    if mark_sites:
        cmd += ["-m", mark_sites]
    if min_overlap:
        cmd += ["--min-overlap", min_overlap]
    if no_version:
        cmd.append("--no-version")
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if pair_logic:
        if pair_logic not in {"snps", "indels", "both", "all", "some", "exact", "id"}:
            raise ValueError(f"Invalid pair-logic value: {pair_logic}")
        cmd += ["--pair-logic", pair_logic]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if rename_annots:
        ra_path = _validate_file_path(rename_annots)
        cmd += ["--rename-annots", str(ra_path)]
    if rename_chrs:
        rc_path = _validate_file_path(rename_chrs)
        cmd += ["--rename-chrs", str(rc_path)]
    if samples:
        cmd += ["-s", samples]
    if samples_file:
        sf_path = _validate_file_path(samples_file)
        cmd += ["-S", str(sf_path)]
    if single_overlaps:
        cmd.append("--single-overlaps")
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if remove:
        cmd += ["-x", remove]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools annotate failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_call(
    file: str,
    no_version: bool = False,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    ploidy: Optional[str] = None,
    ploidy_file: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    targets: Optional[str] = None,
    targets_file: Optional[str] = None,
    targets_overlap: Optional[str] = None,
    threads: int = 0,
    write_index: Optional[str] = None,
    keep_alts: bool = False,
    keep_unseen_allele: bool = False,
    format_fields: Optional[str] = None,
    prior_freqs: Optional[str] = None,
    group_samples: Optional[str] = None,
    gvcf: Optional[str] = None,
    insert_missed: Optional[int] = None,
    keep_masked_ref: bool = False,
    skip_variants: Optional[str] = None,
    variants_only: bool = False,
    consensus_caller: bool = False,
    constrain: Optional[str] = None,
    multiallelic_caller: bool = False,
    novel_rate: Optional[str] = None,
    pval_threshold: Optional[float] = None,
    prior: Optional[float] = None,
    chromosome_X: bool = False,
    chromosome_Y: bool = False,
    verbosity: int = 1,
):
    """
    SNP/indel calling from mpileup output using bcftools call.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "call"]
    if no_version:
        cmd.append("--no-version")
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if ploidy:
        cmd += ["--ploidy", ploidy]
    if ploidy_file:
        pf_path = _validate_file_path(ploidy_file)
        cmd += ["--ploidy-file", str(pf_path)]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if samples:
        cmd += ["-s", samples]
    if samples_file:
        sf_path = _validate_file_path(samples_file)
        cmd += ["-S", str(sf_path)]
    if targets:
        cmd += ["-t", targets]
    if targets_file:
        tf_path = _validate_file_path(targets_file)
        cmd += ["-T", str(tf_path)]
    if targets_overlap:
        if targets_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid targets-overlap value: {targets_overlap}")
        cmd += ["--targets-overlap", targets_overlap]
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]
    if keep_alts:
        cmd.append("-A")
    if keep_unseen_allele:
        cmd.append("-*")
    if format_fields:
        cmd += ["-f", format_fields]
    if prior_freqs:
        cmd += ["-F", prior_freqs]
    if group_samples:
        if group_samples != "-":
            gs_path = _validate_file_path(group_samples)
            cmd += ["-G", str(gs_path)]
        else:
            cmd += ["-G", "-"]
    if gvcf:
        cmd += ["-g", gvcf]
    if insert_missed is not None:
        if insert_missed < 0:
            raise ValueError("insert_missed must be non-negative")
        cmd += ["-i", str(insert_missed)]
    if keep_masked_ref:
        cmd.append("-M")
    if skip_variants:
        if skip_variants not in {"snps", "indels"}:
            raise ValueError(f"Invalid skip-variants value: {skip_variants}")
        cmd += ["-V", skip_variants]
    if variants_only:
        cmd.append("-v")
    if consensus_caller and multiallelic_caller:
        raise ValueError("Options -c and -m are mutually exclusive")
    if consensus_caller:
        cmd.append("-c")
    if constrain:
        if constrain not in {"alleles", "trio"}:
            raise ValueError(f"Invalid constrain value: {constrain}")
        cmd += ["-C", constrain]
    if multiallelic_caller:
        cmd.append("-m")
    if novel_rate:
        cmd += ["-n", novel_rate]
    if pval_threshold is not None:
        if pval_threshold < 0.0:
            raise ValueError("pval_threshold must be non-negative")
        cmd += ["-p", str(pval_threshold)]
    if prior is not None:
        if prior < 0.0:
            raise ValueError("prior must be non-negative")
        cmd += ["-P", str(prior)]
    if chromosome_X:
        cmd.append("-X")
    if chromosome_Y:
        cmd.append("-Y")
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools call failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_view(
    file: str,
    drop_genotypes: bool = False,
    header_only: bool = False,
    no_header: bool = False,
    with_header: bool = False,
    compression_level: Optional[int] = None,
    no_version: bool = False,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    threads: int = 0,
    verbosity: int = 1,
    write_index: Optional[str] = None,
    trim_unseen_alleles: int = 0,
    trim_alt_alleles: bool = False,
    force_samples: bool = False,
    no_update: bool = False,
    min_pq: Optional[int] = None,
    min_ac: Optional[int] = None,
    max_ac: Optional[int] = None,
    exclude: Optional[str] = None,
    apply_filters: Optional[str] = None,
    genotype: Optional[str] = None,
    include: Optional[str] = None,
    known: bool = False,
    min_alleles: Optional[int] = None,
    max_alleles: Optional[int] = None,
    novel: bool = False,
    phased: bool = False,
    exclude_phased: bool = False,
    min_af: Optional[float] = None,
    max_af: Optional[float] = None,
    uncalled: bool = False,
    exclude_uncalled: bool = False,
    types: Optional[str] = None,
    exclude_types: Optional[str] = None,
    private: bool = False,
    exclude_private: bool = False,
):
    """
    View, subset and filter VCF or BCF files by position and filtering expression.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "view"]
    if drop_genotypes:
        cmd.append("-G")
    if header_only:
        cmd.append("-h")
    if no_header:
        cmd.append("-H")
    if with_header:
        cmd.append("--with-header")
    if compression_level is not None:
        if not (0 <= compression_level <= 9):
            raise ValueError("compression_level must be between 0 and 9")
        cmd += ["-l", str(compression_level)]
    if no_version:
        cmd.append("--no-version")
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if samples:
        cmd += ["-s", samples]
    if samples_file:
        sf_path = _validate_file_path(samples_file)
        cmd += ["-S", str(sf_path)]
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]
    if trim_unseen_alleles not in {0, 1, 2}:
        raise ValueError("trim_unseen_alleles must be 0, 1, or 2")
    if trim_unseen_alleles == 1:
        cmd.append("-A")
    elif trim_unseen_alleles == 2:
        cmd.append("-AA")
    if trim_alt_alleles:
        cmd.append("-a")
    if force_samples:
        cmd.append("--force-samples")
    if no_update:
        cmd.append("-I")
    if min_pq is not None:
        if min_pq < 0:
            raise ValueError("min_pq must be non-negative")
        cmd += ["-q", str(min_pq)]
    if min_ac is not None:
        if min_ac < 0:
            raise ValueError("min_ac must be non-negative")
        cmd += ["-c", str(min_ac)]
    if max_ac is not None:
        if max_ac < 0:
            raise ValueError("max_ac must be non-negative")
        cmd += ["-C", str(max_ac)]
    if exclude:
        cmd += ["-e", exclude]
    if apply_filters:
        cmd += ["-f", apply_filters]
    if genotype:
        cmd += ["-g", genotype]
    if include:
        cmd += ["-i", include]
    if known:
        cmd.append("-k")
    if min_alleles is not None:
        if min_alleles < 0:
            raise ValueError("min_alleles must be non-negative")
        cmd += ["-m", str(min_alleles)]
    if max_alleles is not None:
        if max_alleles < 0:
            raise ValueError("max_alleles must be non-negative")
        cmd += ["-M", str(max_alleles)]
    if novel:
        cmd.append("-n")
    if phased:
        cmd.append("-p")
    if exclude_phased:
        cmd.append("-P")
    if min_af is not None:
        if not (0.0 <= min_af <= 1.0):
            raise ValueError("min_af must be between 0 and 1")
        cmd += ["-q", str(min_af)]
    if max_af is not None:
        if not (0.0 <= max_af <= 1.0):
            raise ValueError("max_af must be between 0 and 1")
        cmd += ["-Q", str(max_af)]
    if uncalled:
        cmd.append("-u")
    if exclude_uncalled:
        cmd.append("-U")
    if types:
        cmd += ["-v", types]
    if exclude_types:
        cmd += ["-V", exclude_types]
    if private:
        cmd.append("-x")
    if exclude_private:
        cmd.append("-X")

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools view failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_index(
    file: str,
    csi: bool = True,
    force: bool = False,
    min_shift: int = 14,
    output: Optional[str] = None,
    tbi: bool = False,
    threads: int = 0,
    verbosity: int = 1,
):
    """
    Create index for bgzip compressed VCF/BCF files for random access.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "index"]
    if csi and not tbi:
        cmd.append("-c")
    if force:
        cmd.append("-f")
    if min_shift < 0:
        raise ValueError("min_shift must be non-negative")
    cmd += ["-m", str(min_shift)]
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if tbi:
        cmd.append("-t")
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        else:
            # Default index file name
            if tbi:
                idx_file = file_path.with_suffix(file_path.suffix + ".tbi")
            else:
                idx_file = file_path.with_suffix(file_path.suffix + ".csi")
            if idx_file.exists():
                output_files.append(str(idx_file.resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools index failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_concat(
    files: List[str],
    allow_overlaps: bool = False,
    compact_ps: bool = False,
    rm_dups: Optional[str] = None,
    file_list: Optional[str] = None,
    ligate: bool = False,
    ligate_force: bool = False,
    ligate_warn: bool = False,
    no_version: bool = False,
    naive: bool = False,
    naive_force: bool = False,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    min_pq: Optional[int] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    threads: int = 0,
    verbosity: int = 1,
    write_index: Optional[str] = None,
):
    """
    Concatenate or combine VCF/BCF files with bcftools concat.
    """
    if file_list:
        fl_path = _validate_file_path(file_list)
    else:
        for f in files:
            _validate_file_path(f)
    cmd = ["bcftools", "concat"]
    if allow_overlaps:
        cmd.append("-a")
    if compact_ps:
        cmd.append("-c")
    if rm_dups:
        if rm_dups not in {"snps", "indels", "both", "all", "exact"}:
            raise ValueError(f"Invalid rm_dups value: {rm_dups}")
        cmd += ["-d", rm_dups]
    if file_list:
        cmd += ["-f", str(fl_path)]
    if ligate:
        cmd.append("-l")
    if ligate_force:
        cmd.append("--ligate-force")
    if ligate_warn:
        cmd.append("--ligate-warn")
    if no_version:
        cmd.append("--no-version")
    if naive:
        cmd.append("-n")
    if naive_force:
        cmd.append("--naive-force")
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if min_pq is not None:
        if min_pq < 0:
            raise ValueError("min_pq must be non-negative")
        cmd += ["-q", str(min_pq)]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]

    if not file_list:
        cmd += files

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools concat failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_query(
    file: str,
    exclude: Optional[str] = None,
    force_samples: bool = False,
    format: Optional[str] = None,
    print_filtered: Optional[str] = None,
    print_header: bool = False,
    include: Optional[str] = None,
    list_samples: bool = False,
    disable_automatic_newline: bool = False,
    output: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    allow_undef_tags: bool = False,
    vcf_list: Optional[str] = None,
    verbosity: int = 1,
):
    """
    Extract fields from VCF or BCF files and output in user-defined format using bcftools query.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "query"]
    if exclude:
        cmd += ["-e", exclude]
    if force_samples:
        cmd.append("--force-samples")
    if format:
        cmd += ["-f", format]
    if print_filtered:
        cmd += ["-F", print_filtered]
    if print_header:
        cmd.append("-H")
    if include:
        cmd += ["-i", include]
    if list_samples:
        cmd.append("-l")
    if disable_automatic_newline:
        cmd.append("-N")
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if samples:
        cmd += ["-s", samples]
    if samples_file:
        sf_path = _validate_file_path(samples_file)
        cmd += ["-S", str(sf_path)]
    if allow_undef_tags:
        cmd.append("-u")
    if vcf_list:
        vl_path = _validate_file_path(vcf_list)
        cmd += ["-v", str(vl_path)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools query failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_stats(
    file1: str,
    file2: Optional[str] = None,
    af_bins: Optional[str] = None,
    af_tag: Optional[str] = None,
    all_contigs: bool = False,
    nrecords: bool = False,
    stats: bool = False,
    exclude: Optional[str] = None,
    exons: Optional[str] = None,
    apply_filters: Optional[str] = None,
    fasta_ref: Optional[str] = None,
    include: Optional[str] = None,
    split_by_id: bool = False,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    samples: Optional[str] = None,
    samples_file: Optional[str] = None,
    targets: Optional[str] = None,
    targets_file: Optional[str] = None,
    targets_overlap: Optional[str] = None,
    user_tstv: Optional[str] = None,
    verbosity: int = 1,
):
    """
    Produce VCF/BCF stats using bcftools stats.
    """
    file1_path = _validate_file_path(file1)
    cmd = ["bcftools", "stats"]
    if file2:
        file2_path = _validate_file_path(file2)
    if af_bins:
        cmd += ["--af-bins", af_bins]
    if af_tag:
        cmd += ["--af-tag", af_tag]
    if all_contigs:
        cmd.append("-a")
    if nrecords:
        cmd.append("-n")
    if stats:
        cmd.append("-s")
    if exclude:
        cmd += ["-e", exclude]
    if exons:
        exons_path = _validate_file_path(exons)
        cmd += ["-E", str(exons_path)]
    if apply_filters:
        cmd += ["-f", apply_filters]
    if fasta_ref:
        fasta_path = _validate_file_path(fasta_ref)
        cmd += ["-F", str(fasta_path)]
    if include:
        cmd += ["-i", include]
    if split_by_id:
        cmd.append("-I")
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if samples:
        cmd += ["-s", samples]
    if samples_file:
        sf_path = _validate_file_path(samples_file)
        cmd += ["-S", str(sf_path)]
    if targets:
        cmd += ["-t", targets]
    if targets_file:
        tf_path = _validate_file_path(targets_file)
        cmd += ["-T", str(tf_path)]
    if targets_overlap:
        if targets_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid targets-overlap value: {targets_overlap}")
        cmd += ["--targets-overlap", targets_overlap]
    if user_tstv:
        cmd += ["-u", user_tstv]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]

    cmd.append(str(file1_path))
    if file2:
        cmd.append(str(file2_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools stats failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_sort(
    file: str,
    max_mem: Optional[str] = None,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    temp_dir: Optional[str] = None,
    verbosity: int = 1,
    write_index: Optional[str] = None,
):
    """
    Sort VCF/BCF files using bcftools sort.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", "sort"]
    if max_mem:
        cmd += ["-m", max_mem]
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if temp_dir:
        temp_path = Path(temp_dir)
        cmd += ["-T", str(temp_path)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools sort failed with exit code {e.returncode}",
        }

@mcp.tool()
def bcftools_plugin(
    plugin_name: str,
    file: str,
    plugin_options: Optional[List[str]] = None,
    exclude: Optional[str] = None,
    include: Optional[str] = None,
    regions: Optional[str] = None,
    regions_file: Optional[str] = None,
    regions_overlap: Optional[str] = None,
    output: Optional[str] = None,
    output_type: Optional[str] = None,
    threads: int = 0,
    verbosity: int = 1,
    write_index: Optional[str] = None,
):
    """
    Run a bcftools plugin on a VCF/BCF file.
    """
    file_path = _validate_file_path(file)
    cmd = ["bcftools", f"+{plugin_name}"]
    if exclude:
        cmd += ["-e", exclude]
    if include:
        cmd += ["-i", include]
    if regions:
        cmd += ["-r", regions]
    if regions_file:
        rf_path = _validate_file_path(regions_file)
        cmd += ["-R", str(rf_path)]
    if regions_overlap:
        if regions_overlap not in {"0", "1", "2"}:
            raise ValueError(f"Invalid regions-overlap value: {regions_overlap}")
        cmd += ["--regions-overlap", regions_overlap]
    if output:
        out_path = Path(output)
        cmd += ["-o", str(out_path)]
    if output_type:
        cmd += ["-O", output_type]
    if threads < 0:
        raise ValueError("threads must be >= 0")
    if threads > 0:
        cmd += ["--threads", str(threads)]
    if verbosity < 0:
        raise ValueError("verbosity must be >= 0")
    if verbosity != 1:
        cmd += ["-v", str(verbosity)]
    if write_index:
        if write_index not in {"tbi", "csi"}:
            raise ValueError(f"Invalid write-index format: {write_index}")
        cmd += ["-W", write_index]
    if plugin_options:
        cmd += plugin_options

    cmd.append(str(file_path))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        if output:
            output_files.append(str(Path(output).resolve()))
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"bcftools plugin {plugin_name} failed with exit code {e.returncode}",
        }

if __name__ == '__main__':
    mcp.run()