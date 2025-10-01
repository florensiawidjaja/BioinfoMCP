from fastmcp import FastMCP
from pathlib import Path
import subprocess
from typing import Optional

mcp = FastMCP()

@mcp.tool()
def seqtk_seq(
    input_fasta: str,
    output_fasta: str,
    length: int = 0,
    trim_left: int = 0,
    trim_right: int = 0,
    rev_comp: bool = False,
):
    """
    Extract or trim sequences from FASTA/FASTQ files using seqtk seq command.
    Parameters:
    - input_fasta: Path to input FASTA/FASTQ file.
    - output_fasta: Path to output FASTA/FASTQ file.
    - length: Length to truncate sequences to (0 means no truncation).
    - trim_left: Number of bases to trim from the left.
    - trim_right: Number of bases to trim from the right.
    - rev_comp: Whether to reverse complement the sequences.
    """
    # Validate input file
    input_path = Path(input_fasta)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file {input_fasta} does not exist")

    # Validate output path parent directory
    output_path = Path(output_fasta)
    if not output_path.parent.exists():
        raise FileNotFoundError(f"Output directory {output_path.parent} does not exist")

    # Validate numeric parameters
    if length < 0:
        raise ValueError("length must be >= 0")
    if trim_left < 0:
        raise ValueError("trim_left must be >= 0")
    if trim_right < 0:
        raise ValueError("trim_right must be >= 0")

    # Build command
    cmd = ["seqtk", "seq"]
    if length > 0:
        cmd.append(f"-L{length}")
    if trim_left > 0:
        cmd.append(f"-b{trim_left}")
    if trim_right > 0:
        cmd.append(f"-e{trim_right}")
    if rev_comp:
        cmd.append("-r")
    cmd.append(str(input_path))

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        # Write output to file
        with open(output_path, "w") as f:
            f.write(result.stdout)

        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"seqtk seq failed with return code {e.returncode}",
        }

if __name__ == '__main__':
    mcp.run()