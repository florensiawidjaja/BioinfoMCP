from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, Union
import subprocess
import multiprocessing

mcp = FastMCP()

@mcp.tool()
def correctGCBias(
    bamfile: str,
    effectiveGenomeSize: int,
    genome: str,
    GCbiasFrequenciesFile: str,
    correctedFile: str,
    binSize: int = 50,
    region: Optional[str] = None,
    numberOfProcessors: Union[int, str] = 1,
    verbose: bool = False,
) -> dict:
    """
    Corrects GC-bias in a sorted BAM file using the method by Benjamini & Speed (2012).
    Requires the output of computeGCBias and a genome file in 2bit format.
    
    Parameters:
    - bamfile: Path to the sorted BAM file to correct.
    - effectiveGenomeSize: Effective genome size (mappable portion of the genome).
    - genome: Genome file in 2bit format.
    - GCbiasFrequenciesFile: Output file from computeGCBias with observed and expected read frequencies.
    - correctedFile: Name of the corrected output file (.bam, .bw, or .bg).
    - binSize: Size of bins for bigWig/bedGraph output (default 50).
    - region: Optional genomic region to limit operation (format: chr or chr:start:end).
    - numberOfProcessors: Number of processors to use; can be int or "max"/"max/2" (default 1).
    - verbose: If True, show processing messages.
    
    Returns:
    A dictionary with command executed, stdout, stderr, and list of output files.
    """
    # Validate input files
    bam_path = Path(bamfile)
    if not bam_path.is_file():
        raise FileNotFoundError(f"BAM file not found: {bamfile}")
    genome_path = Path(genome)
    if not genome_path.is_file():
        raise FileNotFoundError(f"Genome 2bit file not found: {genome}")
    freq_path = Path(GCbiasFrequenciesFile)
    if not freq_path.is_file():
        raise FileNotFoundError(f"GCbias frequencies file not found: {GCbiasFrequenciesFile}")

    # Validate correctedFile extension
    corrected_path = Path(correctedFile)
    if corrected_path.suffix not in [".bam", ".bw", ".bg"]:
        raise ValueError("correctedFile must end with .bam, .bw, or .bg")

    # Validate effectiveGenomeSize
    if effectiveGenomeSize <= 0:
        raise ValueError("effectiveGenomeSize must be a positive integer")

    # Validate binSize
    if binSize <= 0:
        raise ValueError("binSize must be a positive integer")

    # Validate numberOfProcessors
    max_cpus = multiprocessing.cpu_count()
    if isinstance(numberOfProcessors, str):
        if numberOfProcessors == "max":
            nproc = max_cpus
        elif numberOfProcessors == "max/2":
            nproc = max_cpus // 2 if max_cpus > 1 else 1
        else:
            raise ValueError("numberOfProcessors string must be 'max' or 'max/2'")
    elif isinstance(numberOfProcessors, int):
        if numberOfProcessors < 1:
            raise ValueError("numberOfProcessors must be at least 1")
        nproc = min(numberOfProcessors, max_cpus)
    else:
        raise TypeError("numberOfProcessors must be int or str ('max'/'max/2')")

    # Build command line
    cmd = [
        "correctGCBias",
        "-b", str(bam_path),
        "--effectiveGenomeSize", str(effectiveGenomeSize),
        "-g", str(genome_path),
        "--GCbiasFrequenciesFile", str(freq_path),
        "-o", str(corrected_path),
        "--binSize", str(binSize),
        "-p", str(nproc)
    ]
    if region:
        cmd.extend(["-r", region])
    if verbose:
        cmd.append("-v")

    # Run subprocess
    try:
        completed = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "output_files": [str(corrected_path)]
    }


if __name__ == '__main__':
    mcp.run()