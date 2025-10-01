import pytest
from fastmcp import FastMCP, Client
from pathlib import Path
import subprocess
from typing import Optional, List
import os

@pytest.fixture
def mcp_server():
    mcp = FastMCP("TestServer")

    @mcp.tool()
    def fastqc(
        input_files: List[Path],
        outdir: Path = Path("fastqc_out"),
        threads: int = 1,
        extract: bool = False,
        noextract: bool = False,
        quiet: bool = False,
        contaminations: Optional[Path] = None,
        format: Optional[str] = None,
        kmers: Optional[int] = None,
        no_group: bool = False,
        no_warn: bool = False,
    ):
        """
        Run FastQC quality control checks on raw sequence data files.

        Parameters:
        - input_files: List of input sequence files (fastq, fasta, bam, sam).
        - outdir: Output directory for FastQC results (default: fastqc_out).
        - threads: Number of files to process in parallel (default: 1).
        - extract: Extract the zipped output files after analysis.
        - noextract: Do not extract the zipped output files.
        - quiet: Suppress all FastQC output.
        - contaminations: Path to a contaminations file to use.
        - format: Force input file format (fastq, bam, sam).
        - kmers: Number of most frequent kmers to report (default is 5 if set).
        - no_group: Disable grouping of bases for per base sequence content.
        - no_warn: Suppress warnings about input file format.
        """

        # Validate input files
        if len(input_files) == 0:
            raise ValueError("At least one input file must be provided.")
        for f in input_files:
            if not f.exists():
                raise FileNotFoundError(f"Input file does not exist: {f}")

        # Validate threads
        if threads < 1:
            raise ValueError("threads must be >= 1")

        # Validate format if provided
        valid_formats = {"fastq", "bam", "sam"}
        if format is not None and format not in valid_formats:
            raise ValueError(f"format must be one of {valid_formats}")

        # Validate kmers if provided
        if kmers is not None and kmers < 1:
            raise ValueError("kmers must be >= 1")

        # Validate outdir path
        outdir = Path(outdir)
        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)

        # Build command
        cmd = ["fastqc"]
        cmd += [str(f) for f in input_files]
        cmd += ["--outdir", str(outdir)]
        cmd += ["--threads", str(threads)]
        if extract:
            cmd.append("--extract")
        if noextract:
            cmd.append("--noextract")
        if quiet:
            cmd.append("--quiet")
        if contaminations is not None:
            contaminations = Path(contaminations)
            if not contaminations.exists():
                raise FileNotFoundError(f"Contaminations file does not exist: {contaminations}")
            cmd += ["--contaminations", str(contaminations)]
        if format is not None:
            cmd += ["--format", format]
        if kmers is not None:
            cmd += ["--kmers", str(kmers)]
        if no_group:
            cmd.append("--nogroup")
        if no_warn:
            cmd.append("--nowarn")

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Collect output files: FastQC creates *_fastqc.html and *_fastqc.zip for each input file
            output_files = []
            for f in input_files:
                base = f.stem
                html_file = outdir / f"{base}_fastqc.html"
                zip_file = outdir / f"{base}_fastqc.zip"
                if html_file.exists():
                    output_files.append(str(html_file))
                if zip_file.exists():
                    output_files.append(str(zip_file))
        except subprocess.CalledProcessError as e:
            return {
                "command_executed": " ".join(cmd),
                "stdout": e.stdout,
                "stderr": e.stderr,
                "output_files": [],
                "error": f"FastQC failed with return code {e.returncode}"
            }

        return {
            "command_executed": " ".join(cmd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_files": output_files
        }
    return mcp

async def test_tool_functionality(mcp_server):
    # Pass the server directly to the Client constructor
    async with Client(mcp_server) as client:
        result = await client.call_tool("fastqc", {"input_files": ["/Users/florensiawidjaja/Documents/BioInfoMCP/SRR097977.fastq"]})
        assert result.data == "Hello, World!"

# if __name__ == '__main__':
#     mcp.run()