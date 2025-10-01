from fastmcp import FastMCP
from typing import List, Optional
from pathlib import Path
import subprocess
import tempfile

mcp = FastMCP()


@mcp.tool()
def macs3_hmmratac(
    input_files: List[Path],
    format: str = "BAMPE",
    outdir: Path = Path("."),
    name: str = "NA",
    blacklist: Optional[Path] = None,
    modelonly: bool = False,
    model: str = "NA",
    training: str = "NA",
    min_frag_p: float = 0.001,
    cutoff_analysis_only: bool = False,
    cutoff_analysis_max: int = 100,
    cutoff_analysis_steps: int = 100,
    hmm_type: str = "gaussian",
    upper: int = 20,
    lower: int = 10,
    prescan_cutoff: float = 1.2,
) -> dict:
    """
    HMMRATAC peak calling algorithm for ATAC-seq data based on Hidden Markov Model.
    Processes paired-end BAMPE or BEDPE input files to identify accessible chromatin regions.
    Outputs narrowPeak format files with accessible regions.

    Parameters:
    - input_files: List of input BAMPE or BEDPE files (gzipped allowed). All must be same format.
    - format: Format of input files, either "BAMPE" or "BEDPE". Default "BAMPE".
    - outdir: Directory to write output files. Default current directory.
    - name: Prefix name for output files. Default "NA".
    - blacklist: Optional BED file of blacklisted regions to exclude fragments.
    - modelonly: If True, only generate HMM model JSON file and quit. Default False.
    - model: JSON file of pre-trained HMM model to use instead of training. Default "NA".
    - training: BED file of custom training regions for HMM training. Default "NA".
    - min_frag_p: Minimum fragment probability threshold (0-1) to include fragments. Default 0.001.
    - cutoff_analysis_only: If True, only run cutoff analysis report and quit. Default False.
    - cutoff_analysis_max: Max cutoff score for cutoff analysis. Default 100.
    - cutoff_analysis_steps: Number of steps for cutoff analysis resolution. Default 100.
    - hmm_type: Emission type for HMM: "gaussian" (default) or "poisson".
    - upper: Upper fold change cutoff for training sites. Default 20.
    - lower: Lower fold change cutoff for training sites. Default 10.
    - prescan_cutoff: Fold change cutoff for prescanning candidate regions (>1). Default 1.2.

    Returns:
    A dict with keys:
    - command_executed: The full command line executed.
    - stdout: Standard output from the command.
    - stderr: Standard error from the command.
    - output_files: List of output file paths generated.
    """
    # Validate input files
    if not input_files or len(input_files) == 0:
        raise ValueError("At least one input file must be provided in input_files.")
    for f in input_files:
        if not f.exists():
            raise FileNotFoundError(f"Input file does not exist: {f}")
    # Validate format
    format_upper = format.upper()
    if format_upper not in ("BAMPE", "BEDPE"):
        raise ValueError(f"Invalid format '{format}'. Must be 'BAMPE' or 'BEDPE'.")
    # Validate outdir
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    # Validate blacklist file if provided
    if blacklist is not None and not blacklist.exists():
        raise FileNotFoundError(f"Blacklist file does not exist: {blacklist}")
    # Validate min_frag_p
    if not (0 <= min_frag_p <= 1):
        raise ValueError(f"min_frag_p must be between 0 and 1, got {min_frag_p}")
    # Validate hmm_type
    hmm_type_lower = hmm_type.lower()
    if hmm_type_lower not in ("gaussian", "poisson"):
        raise ValueError(f"hmm_type must be 'gaussian' or 'poisson', got {hmm_type}")
    # Validate prescan_cutoff
    if prescan_cutoff <= 1:
        raise ValueError(f"prescan_cutoff must be > 1, got {prescan_cutoff}")
    # Validate upper and lower cutoffs
    if lower < 0:
        raise ValueError(f"lower cutoff must be >= 0, got {lower}")
    if upper <= lower:
        raise ValueError(f"upper cutoff must be greater than lower cutoff, got upper={upper}, lower={lower}")
    # Validate cutoff_analysis_max and cutoff_analysis_steps
    if cutoff_analysis_max < 0:
        raise ValueError(f"cutoff_analysis_max must be >= 0, got {cutoff_analysis_max}")
    if cutoff_analysis_steps <= 0:
        raise ValueError(f"cutoff_analysis_steps must be > 0, got {cutoff_analysis_steps}")
    # Validate training file if provided
    if training != "NA":
        training_path = Path(training)
        if not training_path.exists():
            raise FileNotFoundError(f"Training regions file does not exist: {training_path}")

    # Build command line
    cmd = ["macs3", "hmmratac"]
    # Input files
    for f in input_files:
        cmd.extend(["-i", str(f)])
    # Format
    cmd.extend(["-f", format_upper])
    # Output directory
    cmd.extend(["--outdir", str(outdir)])
    # Name prefix
    cmd.extend(["-n", name])
    # Blacklist
    if blacklist is not None:
        cmd.extend(["-e", str(blacklist)])
    # modelonly
    if modelonly:
        cmd.append("--modelonly")
    # model
    if model != "NA":
        cmd.extend(["--model", model])
    # training regions
    if training != "NA":
        cmd.extend(["-t", training])
    # min_frag_p
    cmd.extend(["--min-frag-p", str(min_frag_p)])
    # cutoff_analysis_only
    if cutoff_analysis_only:
        cmd.append("--cutoff-analysis-only")
    # cutoff_analysis_max
    cmd.extend(["--cutoff-analysis-max", str(cutoff_analysis_max)])
    # cutoff_analysis_steps
    cmd.extend(["--cutoff-analysis-steps", str(cutoff_analysis_steps)])
    # hmm_type
    cmd.extend(["--hmm-type", hmm_type_lower])
    # upper cutoff
    cmd.extend(["-u", str(upper)])
    # lower cutoff
    cmd.extend(["-l", str(lower)])
    # prescan cutoff
    cmd.extend(["-c", str(prescan_cutoff)])

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
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": [],
            "error": f"Command failed with return code {e.returncode}",
        }

    # Determine output files
    # The main output is a narrowPeak file named {name}_peaks.narrowPeak in outdir
    peak_file = outdir / f"{name}_peaks.narrowPeak"
    output_files = []
    if peak_file.exists():
        output_files.append(str(peak_file))

    # Also if modelonly or model json is generated, it will be {name}_model.json in outdir
    model_json = outdir / f"{name}_model.json"
    if modelonly or (model != "NA"):
        if model_json.exists():
            output_files.append(str(model_json))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


if __name__ == "__main__":
    mcp.run()