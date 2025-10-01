from fastmcp import FastMCP, Client
from pathlib import Path
from typing import Optional, List
import subprocess
import pytest

@pytest.fixture
def mcp_server():
    mcp = FastMCP("TestServer")


    def _validate_path_exists(path: Path, param_name: str):
        if not path.exists():
            raise FileNotFoundError(f"The file or directory specified for '{param_name}' does not exist: {path}")

    def _list_to_str(items: Optional[List[str]]) -> str:
        if items is None:
            return ""
        return " ".join(items)

    @mcp.tool()
    def star_genome_generate(
        genomeDir: str,
        genomeFastaFiles: List[str],
        runThreadN: int = 1,
        sjdbGTFfile: Optional[str] = None,
        sjdbOverhang: int = 100,
        genomeSAindexNbases: int = 14,
        genomeChrBinNbits: int = 18,
        genomeSAsparseD: int = 1,
        genomeSuffixLengthMax: int = -1,
        sjdbFileChrStartEnd: Optional[List[str]] = None,
        sjdbGTFchrPrefix: Optional[str] = None,
        sjdbGTFfeatureExon: str = "exon",
        sjdbGTFtagExonParentTranscript: str = "transcript_id",
        sjdbGTFtagExonParentGene: str = "gene_id",
        sjdbGTFtagExonParentGeneName: Optional[str] = None,
        sjdbGTFtagExonParentGeneType: Optional[str] = None,
        sjdbScore: int = 2,
        genomeTransformType: Optional[str] = None,
        genomeTransformVCF: Optional[str] = None,
        limitGenomeGenerateRAM: int = 31000000000,
        outFileNamePrefix: str = "./",
        outTmpDir: Optional[str] = None,
        outTmpKeep: str = "None",
    ):
        """
        Generate genome indices for STAR aligner using reference genome FASTA files and optional annotations.
        """
        # Validate inputs
        genome_dir_path = Path(genomeDir)
        if not genome_dir_path.exists():
            raise FileNotFoundError(f"Genome directory '{genomeDir}' does not exist.")
        if not genome_dir_path.is_dir():
            raise ValueError(f"Genome directory '{genomeDir}' is not a directory.")
        for fasta_file in genomeFastaFiles:
            _validate_path_exists(Path(fasta_file), "genomeFastaFiles")
        if sjdbGTFfile:
            _validate_path_exists(Path(sjdbGTFfile), "sjdbGTFfile")
        if sjdbFileChrStartEnd:
            for f in sjdbFileChrStartEnd:
                _validate_path_exists(Path(f), "sjdbFileChrStartEnd")
        if genomeTransformVCF:
            _validate_path_exists(Path(genomeTransformVCF), "genomeTransformVCF")
        if outTmpDir:
            tmp_dir_path = Path(outTmpDir)
            if not tmp_dir_path.exists():
                raise FileNotFoundError(f"Temporary directory '{outTmpDir}' does not exist.")

        # Build command
        cmd = [
            "STAR",
            "--runMode", "genomeGenerate",
            "--genomeDir", genomeDir,
            "--runThreadN", str(runThreadN),
            "--genomeSAindexNbases", str(genomeSAindexNbases),
            "--genomeChrBinNbits", str(genomeChrBinNbits),
            "--genomeSAsparseD", str(genomeSAsparseD),
            "--sjdbScore", str(sjdbScore),
            "--limitGenomeGenerateRAM", str(limitGenomeGenerateRAM),
            "--outFileNamePrefix", outFileNamePrefix,
            "--outTmpKeep", outTmpKeep,
        ]

        # Add genome fasta files
        cmd.append("--genomeFastaFiles")
        cmd.extend(genomeFastaFiles)

        # Optional parameters
        if sjdbGTFfile:
            cmd.extend(["--sjdbGTFfile", sjdbGTFfile])
        if sjdbOverhang > 0:
            cmd.extend(["--sjdbOverhang", str(sjdbOverhang)])
        if sjdbFileChrStartEnd:
            cmd.append("--sjdbFileChrStartEnd")
            cmd.extend(sjdbFileChrStartEnd)
        if sjdbGTFchrPrefix:
            cmd.extend(["--sjdbGTFchrPrefix", sjdbGTFchrPrefix])
        if sjdbGTFfeatureExon:
            cmd.extend(["--sjdbGTFfeatureExon", sjdbGTFfeatureExon])
        if sjdbGTFtagExonParentTranscript:
            cmd.extend(["--sjdbGTFtagExonParentTranscript", sjdbGTFtagExonParentTranscript])
        if sjdbGTFtagExonParentGene:
            cmd.extend(["--sjdbGTFtagExonParentGene", sjdbGTFtagExonParentGene])
        if sjdbGTFtagExonParentGeneName:
            cmd.extend(["--sjdbGTFtagExonParentGeneName", sjdbGTFtagExonParentGeneName])
        if sjdbGTFtagExonParentGeneType:
            cmd.extend(["--sjdbGTFtagExonParentGeneType", sjdbGTFtagExonParentGeneType])
        if genomeSuffixLengthMax != -1:
            cmd.extend(["--genomeSuffixLengthMax", str(genomeSuffixLengthMax)])
        if genomeTransformType:
            cmd.extend(["--genomeTransformType", genomeTransformType])
        if genomeTransformVCF:
            cmd.extend(["--genomeTransformVCF", genomeTransformVCF])
        if outTmpDir:
            cmd.extend(["--outTmpDir", outTmpDir])

        # Run subprocess
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Collect output files: genomeDir contents (list all files)
            output_files = [str(p) for p in genome_dir_path.glob("**/*") if p.is_file()]
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
                "error": f"STAR genomeGenerate failed with return code {e.returncode}",
            }

    @mcp.tool()
    def star_align_reads(
        genomeDir: str,
        readFilesIn: List[str],
        runThreadN: int = 1,
        readFilesCommand: Optional[str] = None,
        outFileNamePrefix: str = "./",
        outSAMtype: Optional[List[str]] = None,
        outSAMattributes: Optional[List[str]] = None,
        outFilterMultimapNmax: int = 10,
        outFilterMismatchNmax: int = 10,
        outFilterMismatchNoverReadLmax: float = 1.0,
        outFilterType: str = "Normal",
        outSAMunmapped: str = "None",
        outSAMstrandField: Optional[str] = None,
        outSAMorder: str = "Paired",
        outSAMprimaryFlag: str = "OneBestScore",
        outSAMmapqUnique: int = 255,
        outSAMmultNmax: int = -1,
        outReadsUnmapped: str = "None",
        genomeLoad: str = "NoSharedMemory",
        sjdbGTFfile: Optional[str] = None,
        sjdbFileChrStartEnd: Optional[List[str]] = None,
        sjdbOverhang: int = 100,
        twopassMode: str = "None",
        limitBAMsortRAM: int = 0,
        outBAMcompression: int = 1,
        outBAMsortingThreadN: int = 0,
        outTmpDir: Optional[str] = None,
        outTmpKeep: str = "None",
        readFilesManifest: Optional[str] = None,
        outSAMattrRGline: Optional[List[str]] = None,
        runRNGseed: int = 777,
    ):
        """
        Run STAR alignment of reads to a genome index directory.
        Supports single or paired-end reads, compressed input, and various output options.
        """
        # Validate inputs
        genome_dir_path = Path(genomeDir)
        _validate_path_exists(genome_dir_path, "genomeDir")

        # Validate read files or manifest
        if readFilesManifest:
            manifest_path = Path(readFilesManifest)
            _validate_path_exists(manifest_path, "readFilesManifest")
        else:
            for rf in readFilesIn:
                _validate_path_exists(Path(rf), "readFilesIn")

        # Build command
        cmd = [
            "STAR",
            "--runMode", "alignReads",
            "--genomeDir", genomeDir,
            "--runThreadN", str(runThreadN),
            "--outFileNamePrefix", outFileNamePrefix,
            "--outFilterMultimapNmax", str(outFilterMultimapNmax),
            "--outFilterMismatchNmax", str(outFilterMismatchNmax),
            "--outFilterMismatchNoverReadLmax", str(outFilterMismatchNoverReadLmax),
            "--outFilterType", outFilterType,
            "--outSAMunmapped", outSAMunmapped,
            "--outSAMorder", outSAMorder,
            "--outSAMprimaryFlag", outSAMprimaryFlag,
            "--outSAMmapqUnique", str(outSAMmapqUnique),
            "--outSAMmultNmax", str(outSAMmultNmax),
            "--outReadsUnmapped", outReadsUnmapped,
            "--genomeLoad", genomeLoad,
            "--runRNGseed", str(runRNGseed),
            "--outTmpKeep", outTmpKeep,
        ]

        # Input reads or manifest
        if readFilesManifest:
            cmd.extend(["--readFilesManifest", readFilesManifest])
        else:
            cmd.append("--readFilesIn")
            cmd.extend(readFilesIn)

        # Optional parameters
        if readFilesCommand:
            cmd.extend(["--readFilesCommand", readFilesCommand])
        if outSAMtype:
            cmd.append("--outSAMtype")
            cmd.extend(outSAMtype)
        if outSAMattributes:
            cmd.append("--outSAMattributes")
            cmd.extend(outSAMattributes)
        if sjdbGTFfile:
            cmd.extend(["--sjdbGTFfile", sjdbGTFfile])
        if sjdbFileChrStartEnd:
            cmd.append("--sjdbFileChrStartEnd")
            cmd.extend(sjdbFileChrStartEnd)
        if sjdbOverhang > 0:
            cmd.extend(["--sjdbOverhang", str(sjdbOverhang)])
        if limitBAMsortRAM > 0:
            cmd.extend(["--limitBAMsortRAM", str(limitBAMsortRAM)])
        if outBAMcompression >= -2 and outBAMcompression <= 10:
            cmd.extend(["--outBAMcompression", str(outBAMcompression)])
        if outBAMsortingThreadN >= 0:
            cmd.extend(["--outBAMsortingThreadN", str(outBAMsortingThreadN)])
        if outTmpDir:
            tmp_dir_path = Path(outTmpDir)
            if not tmp_dir_path.exists():
                raise FileNotFoundError(f"Temporary directory '{outTmpDir}' does not exist.")
            cmd.extend(["--outTmpDir", outTmpDir])
        if outSAMattrRGline:
            # Join RG lines with commas surrounded by spaces as required
            rg_line = " , ".join(outSAMattrRGline)
            cmd.extend(["--outSAMattrRGline", rg_line])
        if outSAMstrandField:
            cmd.extend(["--outSAMstrandField", outSAMstrandField])
        if twopassMode != "None":
            cmd.extend(["--twopassMode", twopassMode])

        # Run subprocess
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Collect output files based on prefix directory
            prefix_path = Path(outFileNamePrefix)
            if prefix_path.is_dir():
                output_files = [str(p) for p in prefix_path.glob("**/*") if p.is_file()]
            else:
                # If prefix is a file prefix, collect files starting with prefix
                parent_dir = prefix_path.parent if prefix_path.parent.exists() else Path(".")
                output_files = [str(p) for p in parent_dir.glob(f"{prefix_path.name}*") if p.is_file()]
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
                "error": f"STAR alignment failed with return code {e.returncode}",
            }

    @mcp.tool()
    def star_remove_genome_from_shared_memory(
        genomeDir: str,
    ):
        """
        Remove a loaded genome from shared memory using STAR --genomeLoad Remove.
        """
        genome_dir_path = Path(genomeDir)
        _validate_path_exists(genome_dir_path, "genomeDir")

        cmd = [
            "STAR",
            "--genomeDir", genomeDir,
            "--genomeLoad", "Remove",
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
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
                "error": f"STAR genomeRemove failed with return code {e.returncode}",
            }
    return mcp
    
async def test_tool_functionality(mcp_server):
    # Pass the server directly to the Client constructor
    async with Client(mcp_server) as client:
        # result = await client.call_tool("fastqc", {"input_files": ["/Users/florensiawidjaja/Documents/BioInfoMCP/SRR097977.fastq"]})
        result = await client.call_tool("fastqc", {"input_files": ["/Users/florensiawidjaja/Documents/BioInfoMCP/SRR097977.fastq"]})
        assert result.data == "Hello, World!"

# if __name__ == '__main__':
    # mcp.run()