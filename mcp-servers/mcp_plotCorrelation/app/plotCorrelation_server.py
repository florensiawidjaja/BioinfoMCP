from fastmcp import FastMCP
from pathlib import Path
import subprocess
from typing import Literal

mcp = FastMCP()

@mcp.tool()
def plotCorrelation(
    corData: Path,
    corMethod: Literal['pearson', 'spearman', 'kendall'],
    whatToPlot: Literal['heatmap', 'scatterplot', 'correlation'],
    output: Path,
) -> dict:
    """
    Generate a correlation plot from a matrix file.

    Parameters:
    - corData: Path to the input matrix file (e.g., matrix.gz).
    - corMethod: Correlation method to use. Allowed values: 'pearson', 'spearman', 'kendall'.
    - whatToPlot: Type of plot to generate. Allowed values: 'heatmap', 'scatterplot', 'correlation'.
    - output: Path to the output plot file (e.g., plot.png).
    """
    # Validate input file
    if not corData.exists():
        raise FileNotFoundError(f"Input correlation data file not found: {corData}")

    # Validate output directory
    output_dir = output.parent
    if not output_dir.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")

    # Construct command
    cmd = [
        "plotCorrelation",
        "-in", str(corData),
        "-c", corMethod,
        "-p", whatToPlot,
        "-o", str(output)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": [str(output)]
        }
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": []
        }

if __name__ == '__main__':
    mcp.run()