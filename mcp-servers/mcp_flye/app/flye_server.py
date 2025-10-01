from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional
import subprocess

mcp = FastMCP()

@mcp.tool()
def flye(
    input_type: str,
    input_files: List[Path],
    out_dir: Path,
    genome_size: Optional[str] = None,
    threads: int = 1,
    iterations: int = 2,
    meta: bool = False,
    polish_target: bool = False,
    min_overlap: Optional[str] = None,
    keep_haplotypes: bool = False,
    debug: bool = False,
    scaffold: bool = False,
    resume: bool = False,
    resume_from: Optional[str] = None,
    stop_after: Optional[str] = None,
    read_error: Optional[float] = None,
    extra_params: Optional[str] = None,
    deterministic: bool = False,
):
    """
    Flye assembler for long reads.

    Parameters:
    - input_type: One of --pacbio-raw, --pacbio-corr, --pacbio-hifi, --nano-raw, --nano-corr, --nano-hq (without dashes)
    - input_files: List of input read files (at least one)
    - out_dir: Output directory path (required)
    - genome_size: Estimated genome size (optional)
    - threads: Number of threads to use (default 1)
    - iterations: Number of assembly iterations (default 2)
    - meta: Enable metagenome mode (default False)
    - polish_target: Enable polish target mode (default False)
    - min_overlap: Minimum overlap size (optional)
    - keep_haplotypes: Keep haplotypes (default False)
    - debug: Enable debug mode (default False)
    - scaffold: Enable scaffolding (default False)
    - resume: Resume previous run (default False)
    - resume_from: Resume from specific step (optional)
    - stop_after: Stop after specific step (optional)
    - read_error: Read error rate (float, optional)
    - extra_params: Extra parameters as string (optional)
    - deterministic: Enable deterministic mode (default False)
    """
    # Validate input_type
    valid_input_types = {
        "pacbio-raw": "--pacbio-raw",
        "pacbio-corr": "--pacbio-corr",
        "pacbio-hifi": "--pacbio-hifi",
        "nano-raw": "--nano-raw",
        "nano-corr": "--nano-corr",
        "nano-hq": "--nano-hq"
    }
    if input_type not in valid_input_types:
        raise ValueError(f"Invalid input_type '{input_type}'. Must be one of {list(valid_input_types.keys())}")

    # Validate input_files
    if not input_files or len(input_files) == 0:
        raise ValueError("At least one input file must be provided in input_files")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file does not exist: {f}")

    # Validate out_dir
    if not isinstance(out_dir, Path):
        raise TypeError("out_dir must be a pathlib.Path object")
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    # Validate threads
    if threads < 1:
        raise ValueError("threads must be >= 1")

    # Validate iterations
    if iterations < 1:
        raise ValueError("iterations must be >= 1")

    # Validate read_error if provided
    if read_error is not None:
        if not (0.0 <= read_error <= 1.0):
            raise ValueError("read_error must be between 0.0 and 1.0")

    # Build command
    cmd = ["flye"]
    cmd.append(valid_input_types[input_type])
    for f in input_files:
        cmd.append(str(f))
    cmd.extend(["--out-dir", str(out_dir)])
    if genome_size:
        cmd.extend(["--genome-size", genome_size])
    cmd.extend(["--threads", str(threads)])
    cmd.extend(["--iterations", str(iterations)])
    if meta:
        cmd.append("--meta")
    if polish_target:
        cmd.append("--polish-target")
    if min_overlap:
        cmd.extend(["--min-overlap", min_overlap])
    if keep_haplotypes:
        cmd.append("--keep-haplotypes")
    if debug:
        cmd.append("--debug")
    if scaffold:
        cmd.append("--scaffold")
    if resume:
        cmd.append("--resume")
    if resume_from:
        cmd.extend(["--resume-from", resume_from])
    if stop_after:
        cmd.extend(["--stop-after", stop_after])
    if read_error is not None:
        cmd.extend(["--read-error", str(read_error)])
    if extra_params:
        # Split extra_params by spaces to allow multiple extra params
        extra_params_split = extra_params.strip().split()
        cmd.extend(extra_params_split)
    if deterministic:
        cmd.append("--deterministic")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": [],
            "error": f"Flye execution failed with return code {e.returncode}"
        }

    # Collect output files - Flye outputs multiple files in out_dir, but we cannot enumerate all.
    # Return the out_dir path as output location.
    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(out_dir)]
    }


if __name__ == '__main__':
    mcp.run()