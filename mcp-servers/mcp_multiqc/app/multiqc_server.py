from fastmcp import FastMCP
from pathlib import Path
import subprocess
from typing import Optional

mcp = FastMCP()

@mcp.tool()
def multiqc(
    analysis_directory: Optional[Path] = None,
    outdir: Optional[Path] = None,
    filename: str = "multiqc_report.html",
    force: bool = False,
    config_file: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    no_data_dir: bool = False,
    no_report: bool = False,
    no_plots: bool = False,
    no_config: bool = False,
    no_title: bool = False,
    title: Optional[str] = None,
    ignore_dirs: Optional[str] = None,
    ignore_samples: Optional[str] = None,
    exclude_modules: Optional[str] = None,
    include_modules: Optional[str] = None,
    verbose: bool = False,
):
    """
    Run MultiQC v1.29 to aggregate results from bioinformatics analyses.

    Parameters:
    - analysis_directory: Directory to scan for analysis results (default: current directory)
    - outdir: Output directory for the MultiQC report (default: current directory)
    - filename: Name of the output report file (default: multiqc_report.html)
    - force: Overwrite existing output files
    - config_file: Path to a custom MultiQC config file
    - data_dir: Path to a directory containing MultiQC data files
    - no_data_dir: Do not use the MultiQC data directory
    - no_report: Do not generate the HTML report
    - no_plots: Do not generate plots
    - no_config: Do not load config files
    - no_title: Do not add a title to the report
    - title: Custom title for the report
    - ignore_dirs: Comma-separated list of directories to ignore
    - ignore_samples: Comma-separated list of samples to ignore
    - exclude_modules: Comma-separated list of modules to exclude
    - include_modules: Comma-separated list of modules to include
    - verbose: Enable verbose output
    """
    import shlex

    # Validate paths
    if analysis_directory is not None:
        if not analysis_directory.exists() or not analysis_directory.is_dir():
            raise ValueError(f"Analysis directory '{analysis_directory}' does not exist or is not a directory.")
    else:
        analysis_directory = Path.cwd()

    if outdir is not None:
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)
    else:
        outdir = Path.cwd()

    if config_file is not None and not config_file.exists():
        raise ValueError(f"Config file '{config_file}' does not exist.")

    if data_dir is not None and not data_dir.exists():
        raise ValueError(f"Data directory '{data_dir}' does not exist.")

    # Build command
    cmd = ["multiqc"]

    # Add analysis directory
    cmd.append(str(analysis_directory))

    # Output directory
    cmd.extend(["-o", str(outdir)])

    # Filename
    if filename:
        cmd.extend(["-n", filename])

    # Flags
    if force:
        cmd.append("-f")
    if config_file:
        cmd.extend(["-c", str(config_file)])
    if data_dir:
        cmd.extend(["--data-dir", str(data_dir)])
    if no_data_dir:
        cmd.append("--no-data-dir")
    if no_report:
        cmd.append("--no-report")
    if no_plots:
        cmd.append("--no-plots")
    if no_config:
        cmd.append("--no-config")
    if no_title:
        cmd.append("--no-title")
    if title:
        cmd.extend(["-t", title])
    if ignore_dirs:
        cmd.extend(["--ignore-dir", ignore_dirs])
    if ignore_samples:
        cmd.extend(["--ignore-samples", ignore_samples])
    if exclude_modules:
        cmd.extend(["--exclude", exclude_modules])
    if include_modules:
        cmd.extend(["--include", include_modules])
    if verbose:
        cmd.append("-v")

    # Run subprocess
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # Collect output files: the main report file in outdir
        output_report = outdir / filename
        output_files = []
        if output_report.exists():
            output_files.append(str(output_report.resolve()))
        return {
            "command_executed": " ".join(shlex.quote(c) for c in cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files,
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(shlex.quote(c) for c in cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"MultiQC failed with exit code {e.returncode}"
        }


if __name__ == '__main__':
    mcp.run()