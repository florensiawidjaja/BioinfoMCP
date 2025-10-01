from fastmcp import FastMCP
from typing import Optional, List, Literal
from pathlib import Path
import subprocess

mcp = FastMCP()

def validate_file_paths(files: List[Path]):
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"Input file does not exist: {f}")

def validate_kmers(kmers: List[int]):
    if not kmers:
        return
    if any(k % 2 == 0 for k in kmers):
        raise ValueError("All k-mer sizes must be odd integers.")
    if any(k >= 128 for k in kmers):
        raise ValueError("All k-mer sizes must be less than 128.")
    if sorted(kmers) != kmers:
        raise ValueError("k-mer sizes must be listed in ascending order.")

def build_library_args(
    # Single library
    interlaced_12: Optional[List[Path]] = None,
    forward_1: Optional[List[Path]] = None,
    reverse_2: Optional[List[Path]] = None,
    merged: Optional[List[Path]] = None,
    unpaired_s: Optional[List[Path]] = None,
    # Multiple single-read libraries
    single_reads: Optional[List[List[Path]]] = None,  # --s1 ... --s9 ...
    # Multiple paired-end libraries
    pe_interlaced_12: Optional[List[List[Path]]] = None,  # --pe1-12 ... --pe9-12
    pe_left_1: Optional[List[List[Path]]] = None,  # --pe1-1 ... --pe9-1
    pe_right_2: Optional[List[List[Path]]] = None,  # --pe1-2 ... --pe9-2
    pe_merged_m: Optional[List[List[Path]]] = None,  # --pe1-m ... --pe9-m
    pe_unpaired_s: Optional[List[List[Path]]] = None,  # --pe1-s ... --pe9-s
    pe_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,  # --pe1-fr etc
    # Mate-pair libraries
    mp_interlaced_12: Optional[List[List[Path]]] = None,
    mp_left_1: Optional[List[List[Path]]] = None,
    mp_right_2: Optional[List[List[Path]]] = None,
    mp_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,
    # High-quality mate-pair libraries
    hqmp_interlaced_12: Optional[List[List[Path]]] = None,
    hqmp_left_1: Optional[List[List[Path]]] = None,
    hqmp_right_2: Optional[List[List[Path]]] = None,
    hqmp_unpaired_s: Optional[List[List[Path]]] = None,
    hqmp_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,
    # Hybrid assembly long reads
    pacbio: Optional[List[Path]] = None,
    nanopore: Optional[List[Path]] = None,
    sanger: Optional[List[Path]] = None,
    # Contigs
    trusted_contigs: Optional[List[Path]] = None,
    untrusted_contigs: Optional[List[Path]] = None,
) -> List[str]:
    """
    Build command line arguments for all library input options.
    Lists of lists are used for multiple libraries (index 0 = library 1).
    """
    args = []

    # Validate and add single library inputs
    if interlaced_12:
        validate_file_paths(interlaced_12)
        for f in interlaced_12:
            args.extend(["--12", str(f)])
    if forward_1:
        validate_file_paths(forward_1)
        for f in forward_1:
            args.extend(["-1", str(f)])
    if reverse_2:
        validate_file_paths(reverse_2)
        for f in reverse_2:
            args.extend(["-2", str(f)])
    if merged:
        validate_file_paths(merged)
        for f in merged:
            args.extend(["--merged", str(f)])
    if unpaired_s:
        validate_file_paths(unpaired_s)
        for f in unpaired_s:
            args.extend(["-s", str(f)])

    # Multiple single-read libraries --s1 ... --s9
    if single_reads:
        if len(single_reads) > 9:
            raise ValueError("Maximum 9 single-read libraries supported (--s1 to --s9).")
        for i, lib_files in enumerate(single_reads, start=1):
            validate_file_paths(lib_files)
            for f in lib_files:
                args.extend([f"--s{i}", str(f)])

    # Multiple paired-end libraries
    # Each is a list of libraries, each library is a list of files
    def add_multi_lib_files(prefix: str, libs: Optional[List[List[Path]]], suffix: str):
        if libs:
            if len(libs) > 9:
                raise ValueError(f"Maximum 9 libraries supported for {prefix}{suffix}.")
            for i, lib_files in enumerate(libs, start=1):
                validate_file_paths(lib_files)
                for f in lib_files:
                    args.extend([f"--{prefix}{i}{suffix}", str(f)])

    add_multi_lib_files("pe", pe_interlaced_12, "-12")
    add_multi_lib_files("pe", pe_left_1, "-1")
    add_multi_lib_files("pe", pe_right_2, "-2")
    add_multi_lib_files("pe", pe_merged_m, "-m")
    add_multi_lib_files("pe", pe_unpaired_s, "-s")

    # Paired-end orientation
    if pe_orientation:
        if len(pe_orientation) > 9:
            raise ValueError("Maximum 9 paired-end libraries orientation supported.")
        for i, orient in enumerate(pe_orientation, start=1):
            if orient not in ("fr", "rf", "ff"):
                raise ValueError(f"Invalid paired-end orientation: {orient}")
            args.append(f"--pe{i}-{orient}")

    # Mate-pair libraries
    add_multi_lib_files("mp", mp_interlaced_12, "-12")
    add_multi_lib_files("mp", mp_left_1, "-1")
    add_multi_lib_files("mp", mp_right_2, "-2")

    if mp_orientation:
        if len(mp_orientation) > 9:
            raise ValueError("Maximum 9 mate-pair libraries orientation supported.")
        for i, orient in enumerate(mp_orientation, start=1):
            if orient not in ("fr", "rf", "ff"):
                raise ValueError(f"Invalid mate-pair orientation: {orient}")
            args.append(f"--mp{i}-{orient}")

    # High-quality mate-pair libraries
    add_multi_lib_files("hqmp", hqmp_interlaced_12, "-12")
    add_multi_lib_files("hqmp", hqmp_left_1, "-1")
    add_multi_lib_files("hqmp", hqmp_right_2, "-2")
    add_multi_lib_files("hqmp", hqmp_unpaired_s, "-s")

    if hqmp_orientation:
        if len(hqmp_orientation) > 9:
            raise ValueError("Maximum 9 high-quality mate-pair libraries orientation supported.")
        for i, orient in enumerate(hqmp_orientation, start=1):
            if orient not in ("fr", "rf", "ff"):
                raise ValueError(f"Invalid high-quality mate-pair orientation: {orient}")
            args.append(f"--hqmp{i}-{orient}")

    # Hybrid assembly long reads
    if pacbio:
        validate_file_paths(pacbio)
        for f in pacbio:
            args.extend(["--pacbio", str(f)])
    if nanopore:
        validate_file_paths(nanopore)
        for f in nanopore:
            args.extend(["--nanopore", str(f)])
    if sanger:
        validate_file_paths(sanger)
        for f in sanger:
            args.extend(["--sanger", str(f)])

    # Contigs
    if trusted_contigs:
        validate_file_paths(trusted_contigs)
        for f in trusted_contigs:
            args.extend(["--trusted-contigs", str(f)])
    if untrusted_contigs:
        validate_file_paths(untrusted_contigs)
        for f in untrusted_contigs:
            args.extend(["--untrusted-contigs", str(f)])

    return args

@mcp.tool()
def spades(
    output_dir: Path,
    # Running modes (flags)
    isolate: bool = False,
    sc: bool = False,
    meta: bool = False,
    plasmid: bool = False,
    metaplasmid: bool = False,
    metaviral: bool = False,
    bio: bool = False,
    rna: bool = False,
    rnaviral: bool = False,
    corona: bool = False,
    iontorrent: bool = False,
    sewage: bool = False,
    # Basic options
    test: bool = False,
    # Pipeline options
    only_error_correction: bool = False,
    only_assembler: bool = False,
    careful: bool = False,
    continue_run: bool = False,
    restart_from: Optional[Literal["ec","as","mc","last"]] = None,
    restart_from_k: Optional[int] = None,  # for k<int> restart
    checkpoints: Literal["none","all","last"] = "none",
    disable_gzip_output: bool = False,
    # Input data - single library
    interlaced_12: Optional[List[Path]] = None,
    forward_1: Optional[List[Path]] = None,
    reverse_2: Optional[List[Path]] = None,
    merged: Optional[List[Path]] = None,
    unpaired_s: Optional[List[Path]] = None,
    # Multiple single-read libraries
    single_reads: Optional[List[List[Path]]] = None,
    # Multiple paired-end libraries
    pe_interlaced_12: Optional[List[List[Path]]] = None,
    pe_left_1: Optional[List[List[Path]]] = None,
    pe_right_2: Optional[List[List[Path]]] = None,
    pe_merged_m: Optional[List[List[Path]]] = None,
    pe_unpaired_s: Optional[List[List[Path]]] = None,
    pe_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,
    # Mate-pair libraries
    mp_interlaced_12: Optional[List[List[Path]]] = None,
    mp_left_1: Optional[List[List[Path]]] = None,
    mp_right_2: Optional[List[List[Path]]] = None,
    mp_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,
    # High-quality mate-pair libraries
    hqmp_interlaced_12: Optional[List[List[Path]]] = None,
    hqmp_left_1: Optional[List[List[Path]]] = None,
    hqmp_right_2: Optional[List[List[Path]]] = None,
    hqmp_unpaired_s: Optional[List[List[Path]]] = None,
    hqmp_orientation: Optional[List[Literal["fr","rf","ff"]]] = None,
    # Hybrid assembly long reads
    pacbio: Optional[List[Path]] = None,
    nanopore: Optional[List[Path]] = None,
    sanger: Optional[List[Path]] = None,
    # Contigs
    trusted_contigs: Optional[List[Path]] = None,
    untrusted_contigs: Optional[List[Path]] = None,
    # Other input
    assembly_graph: Optional[Path] = None,
    # Advanced options
    threads: int = 16,
    memory: int = 250,
    tmp_dir: Optional[Path] = None,
    kmers: Optional[List[int]] = None,
    cov_cutoff: str = "off",  # "off", "auto" or positive float as string
    phred_offset: Optional[int] = None,  # 33 or 64
    custom_hmms: Optional[Path] = None,
    gfa11: bool = False,
    # Dataset YAML file alternative input
    dataset: Optional[Path] = None,
):
    """
    Run SPAdes assembler with comprehensive options for isolate, single-cell, metagenomic,
    plasmid, RNA, viral, corona, IonTorrent, sewage modes and multiple library types.
    """
    # Validate output_dir
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
    if not output_dir.is_dir():
        raise ValueError(f"Output directory path is not a directory: {output_dir}")

    # Validate mutually exclusive modes
    mode_flags = [isolate, sc, meta, plasmid, metaplasmid, metaviral, bio, rna, rnaviral, corona, iontorrent, sewage]
    if sum(mode_flags) > 1:
        raise ValueError("Only one running mode flag can be set among isolate, sc, meta, plasmid, metaplasmid, metaviral, bio, rna, rnaviral, corona, iontorrent, sewage.")

    # Validate incompatible options
    if isolate and (only_error_correction or careful):
        raise ValueError("--isolate is not compatible with --only-error-correction or --careful options.")
    if meta and careful:
        raise ValueError("--meta mode does not support --careful option.")
    if rna and (only_error_correction or careful):
        raise ValueError("--rna mode is not compatible with --only-error-correction or --careful options.")
    if rnaviral and (only_error_correction or careful):
        raise ValueError("--rnaviral mode is not compatible with --only-error-correction or --careful options.")
    if plasmid and sc:
        raise ValueError("--plasmid mode is not compatible with --sc mode.")
    if bio and any(mode_flags) and bio not in mode_flags:
        raise ValueError("--bio mode is not compatible with any other modes.")

    # Validate restart_from_k if restart_from is k<int>
    if restart_from_k is not None:
        if restart_from is not None:
            raise ValueError("Specify either restart_from or restart_from_k, not both.")
        if restart_from_k <= 0:
            raise ValueError("restart_from_k must be a positive integer.")

    # Validate cov_cutoff
    if cov_cutoff != "off" and cov_cutoff != "auto":
        try:
            cov_val = float(cov_cutoff)
            if cov_val <= 0:
                raise ValueError("cov_cutoff must be positive float, 'auto' or 'off'.")
        except Exception:
            raise ValueError("cov_cutoff must be positive float, 'auto' or 'off'.")

    # Validate phred_offset
    if phred_offset is not None and phred_offset not in (33, 64):
        raise ValueError("phred_offset must be either 33 or 64 if specified.")

    # Validate threads and memory
    if threads <= 0:
        raise ValueError("threads must be positive integer.")
    if memory <= 0:
        raise ValueError("memory must be positive integer.")

    # Validate kmers
    if kmers is not None:
        validate_kmers(kmers)

    # Validate dataset YAML file
    if dataset is not None:
        if any([
            interlaced_12, forward_1, reverse_2, merged, unpaired_s,
            single_reads, pe_interlaced_12, pe_left_1, pe_right_2, pe_merged_m, pe_unpaired_s, pe_orientation,
            mp_interlaced_12, mp_left_1, mp_right_2, mp_orientation,
            hqmp_interlaced_12, hqmp_left_1, hqmp_right_2, hqmp_unpaired_s, hqmp_orientation,
            pacbio, nanopore, sanger, trusted_contigs, untrusted_contigs
        ]):
            raise ValueError("--dataset option cannot be used with any other input data options.")
        if not dataset.exists():
            raise FileNotFoundError(f"Dataset YAML file does not exist: {dataset}")

    # Build command line
    cmd = ["spades.py"]

    # Running modes flags
    if isolate:
        cmd.append("--isolate")
    if sc:
        cmd.append("--sc")
    if meta:
        cmd.append("--meta")
    if plasmid:
        cmd.append("--plasmid")
    if metaplasmid:
        cmd.append("--metaplasmid")
    if metaviral:
        cmd.append("--metaviral")
    if bio:
        cmd.append("--bio")
    if rna:
        cmd.append("--rna")
    if rnaviral:
        cmd.append("--rnaviral")
    if corona:
        cmd.append("--corona")
    if iontorrent:
        cmd.append("--iontorrent")
    if sewage:
        cmd.append("--sewage")

    # Basic options
    if test:
        cmd.append("--test")

    # Pipeline options
    if only_error_correction:
        cmd.append("--only-error-correction")
    if only_assembler:
        cmd.append("--only-assembler")
    if careful:
        cmd.append("--careful")
    if continue_run:
        cmd.append("--continue")
    if restart_from is not None:
        if restart_from_k is not None:
            # restart from k<int>
            cmd.append(f"--restart-from k{restart_from_k}")
        else:
            cmd.append(f"--restart-from {restart_from}")

    # Checkpoints
    if checkpoints not in ("none", "all", "last"):
        raise ValueError("checkpoints must be one of 'none', 'all', or 'last'.")
    if checkpoints != "none":
        cmd.extend(["--checkpoints", checkpoints])

    # Disable gzip output
    if disable_gzip_output:
        cmd.append("--disable-gzip-output")

    # Input data or dataset
    if dataset is not None:
        cmd.extend(["--dataset", str(dataset)])
    else:
        lib_args = build_library_args(
            interlaced_12=interlaced_12,
            forward_1=forward_1,
            reverse_2=reverse_2,
            merged=merged,
            unpaired_s=unpaired_s,
            single_reads=single_reads,
            pe_interlaced_12=pe_interlaced_12,
            pe_left_1=pe_left_1,
            pe_right_2=pe_right_2,
            pe_merged_m=pe_merged_m,
            pe_unpaired_s=pe_unpaired_s,
            pe_orientation=pe_orientation,
            mp_interlaced_12=mp_interlaced_12,
            mp_left_1=mp_left_1,
            mp_right_2=mp_right_2,
            mp_orientation=mp_orientation,
            hqmp_interlaced_12=hqmp_interlaced_12,
            hqmp_left_1=hqmp_left_1,
            hqmp_right_2=hqmp_right_2,
            hqmp_unpaired_s=hqmp_unpaired_s,
            hqmp_orientation=hqmp_orientation,
            pacbio=pacbio,
            nanopore=nanopore,
            sanger=sanger,
            trusted_contigs=trusted_contigs,
            untrusted_contigs=untrusted_contigs,
        )
        cmd.extend(lib_args)

    # Other input
    if assembly_graph is not None:
        if not assembly_graph.exists():
            raise FileNotFoundError(f"Assembly graph file does not exist: {assembly_graph}")
        cmd.extend(["--assembly-graph", str(assembly_graph)])

    # Advanced options
    cmd.extend(["-t", str(threads)])
    cmd.extend(["-m", str(memory)])
    if tmp_dir is not None:
        if not tmp_dir.exists():
            tmp_dir.mkdir(parents=True, exist_ok=True)
        cmd.extend(["--tmp-dir", str(tmp_dir)])
    if kmers is not None and len(kmers) > 0:
        kmer_str = ",".join(str(k) for k in kmers)
        cmd.extend(["-k", kmer_str])
    if cov_cutoff != "off":
        cmd.extend(["--cov-cutoff", cov_cutoff])
    if phred_offset is not None:
        cmd.extend(["--phred-offset", str(phred_offset)])
    if custom_hmms is not None:
        if not custom_hmms.exists():
            raise FileNotFoundError(f"Custom HMMs file or directory does not exist: {custom_hmms}")
        cmd.extend(["--custom-hmms", str(custom_hmms)])
    if gfa11:
        cmd.append("--gfa11")

    # Output directory (required)
    cmd.extend(["-o", str(output_dir)])

    # Run subprocess
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = proc.stdout
        stderr = proc.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "error": f"SPAdes execution failed with return code {e.returncode}",
            "output_files": []
        }

    # Collect output files - typical SPAdes output files in output_dir
    output_files = []
    contigs = output_dir / "contigs.fasta"
    scaffolds = output_dir / "scaffolds.fasta"
    assembly_graph_gfa = output_dir / "assembly_graph.gfa"
    if contigs.exists():
        output_files.append(str(contigs))
    if scaffolds.exists():
        output_files.append(str(scaffolds))
    if assembly_graph_gfa.exists():
        output_files.append(str(assembly_graph_gfa))

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }

if __name__ == '__main__':
    mcp.run()