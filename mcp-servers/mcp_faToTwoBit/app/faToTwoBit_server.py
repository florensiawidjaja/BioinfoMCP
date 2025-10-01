from fastmcp import FastMCP
from pathlib import Path
import subprocess
from typing import List, Optional

mcp = FastMCP()

@mcp.tool()
def faToTwoBit(
    input_fasta: List[Path],
    output_2bit: Path,
    long: bool = False,
    noMask: bool = False,
    stripVersion: bool = False,
    ignoreDups: bool = False,
    namePrefix: Optional[str] = None,
):
    """
    Convert DNA sequences from fasta format to 2bit format.

    Parameters:
    - input_fasta: List of input fasta files (at least one required).
    - output_2bit: Output 2bit file path.
    - long: Use 64-bit offsets for index to allow >4Gb sequences (not compatible with older code).
    - noMask: Ignore lower-case masking in fasta file.
    - stripVersion: Strip off version number after '.' for GenBank accessions.
    - ignoreDups: Convert first sequence only if duplicate sequence names exist.
    - namePrefix: Add prefix to start of sequence names in 2bit file.
    """
    # Validate input fasta files
    if len(input_fasta) == 0:
        raise ValueError("At least one input fasta file must be provided.")
    for fasta_file in input_fasta:
        if not fasta_file.is_file():
            raise FileNotFoundError(f"Input fasta file not found: {fasta_file}")

    # Validate output path parent directory exists
    if not output_2bit.parent.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output_2bit.parent}")

    # Build command
    cmd = ["faToTwoBit"]
    if long:
        cmd.append("-long")
    if noMask:
        cmd.append("-noMask")
    if stripVersion:
        cmd.append("-stripVersion")
    if ignoreDups:
        cmd.append("-ignoreDups")
    if namePrefix is not None:
        if not isinstance(namePrefix, str) or len(namePrefix) == 0:
            raise ValueError("namePrefix must be a non-empty string if provided.")
        cmd.append(f"-namePrefix={namePrefix}")

    # Add input fasta files
    cmd.extend(str(f) for f in input_fasta)
    # Add output file
    cmd.append(str(output_2bit))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output_2bit)]
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else "",
            "output_files": []
        }

if __name__ == '__main__':
    mcp.run()