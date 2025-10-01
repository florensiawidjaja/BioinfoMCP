from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


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
    genomeTransformVCF: Optional[str] = None,
    genomeTransformType: str = "None",
    sjdbFileChrStartEnd: Optional[List[str]] = None,
    sjdbGTFchrPrefix: Optional[str] = None,
    sjdbGTFfeatureExon: str = "exon",
    sjdbGTFtagExonParentTranscript: str = "transcript_id",
    sjdbGTFtagExonParentGene: str = "gene_id",
    sjdbGTFtagExonParentGeneName: Optional[List[str]] = None,
    sjdbGTFtagExonParentGeneType: Optional[List[str]] = None,
    sjdbScore: int = 2,
    sjdbInsertSave: str = "Basic",
    limitGenomeGenerateRAM: int = 31000000000,
    parametersFiles: Optional[str] = None,
    sysShell: Optional[str] = None,
):
    """
    Generate genome indices for STAR aligner.
    This step processes reference genome FASTA files and optional annotations (GTF/GFF)
    to create genome indices used for read mapping.
    """
    # Validate paths
    genome_dir_path = Path(genomeDir)
    if not genome_dir_path.exists():
        raise FileNotFoundError(f"Genome directory {genomeDir} does not exist.")
    if not genome_dir_path.is_dir():
        raise ValueError(f"Genome directory {genomeDir} is not a directory.")
    for fasta in genomeFastaFiles:
        fasta_path = Path(fasta)
        if not fasta_path.is_file():
            raise FileNotFoundError(f"Genome FASTA file {fasta} does not exist.")
    if sjdbGTFfile:
        gtf_path = Path(sjdbGTFfile)
        if not gtf_path.is_file():
            raise FileNotFoundError(f"GTF file {sjdbGTFfile} does not exist.")
    if genomeTransformVCF:
        vcf_path = Path(genomeTransformVCF)
        if not vcf_path.is_file():
            raise FileNotFoundError(f"VCF file {genomeTransformVCF} does not exist.")
    if sjdbFileChrStartEnd:
        for sjdb_file in sjdbFileChrStartEnd:
            sjdb_path = Path(sjdb_file)
            if not sjdb_path.is_file():
                raise FileNotFoundError(f"SJDB file {sjdb_file} does not exist.")

    # Validate genomeTransformType
    valid_transform_types = {"None", "Haploid", "Diploid"}
    if genomeTransformType not in valid_transform_types:
        raise ValueError(f"Invalid genomeTransformType: {genomeTransformType}. Must be one of {valid_transform_types}.")

    # Validate sjdbInsertSave
    valid_sjdbInsertSave = {"Basic", "All"}
    if sjdbInsertSave not in valid_sjdbInsertSave:
        raise ValueError(f"Invalid sjdbInsertSave: {sjdbInsertSave}. Must be one of {valid_sjdbInsertSave}.")

    # Build command
    cmd = ["STAR"]
    cmd += ["--runMode", "genomeGenerate"]
    cmd += ["--genomeDir", str(genome_dir_path)]
    cmd += ["--genomeFastaFiles"] + genomeFastaFiles
    cmd += ["--runThreadN", str(runThreadN)]
    cmd += ["--sjdbOverhang", str(sjdbOverhang)]
    cmd += ["--genomeSAindexNbases", str(genomeSAindexNbases)]
    cmd += ["--genomeChrBinNbits", str(genomeChrBinNbits)]
    cmd += ["--genomeSAsparseD", str(genomeSAsparseD)]
    cmd += ["--genomeSuffixLengthMax", str(genomeSuffixLengthMax)]
    cmd += ["--sjdbScore", str(sjdbScore)]
    cmd += ["--sjdbInsertSave", sjdbInsertSave]
    cmd += ["--limitGenomeGenerateRAM", str(limitGenomeGenerateRAM)]

    if sjdbGTFfile:
        cmd += ["--sjdbGTFfile", sjdbGTFfile]
    if sjdbGTFchrPrefix:
        cmd += ["--sjdbGTFchrPrefix", sjdbGTFchrPrefix]
    if sjdbGTFfeatureExon:
        cmd += ["--sjdbGTFfeatureExon", sjdbGTFfeatureExon]
    if sjdbGTFtagExonParentTranscript:
        cmd += ["--sjdbGTFtagExonParentTranscript", sjdbGTFtagExonParentTranscript]
    if sjdbGTFtagExonParentGene:
        cmd += ["--sjdbGTFtagExonParentGene", sjdbGTFtagExonParentGene]
    if sjdbGTFtagExonParentGeneName:
        cmd += ["--sjdbGTFtagExonParentGeneName"] + sjdbGTFtagExonParentGeneName
    if sjdbGTFtagExonParentGeneType:
        cmd += ["--sjdbGTFtagExonParentGeneType"] + sjdbGTFtagExonParentGeneType
    if sjdbFileChrStartEnd:
        cmd += ["--sjdbFileChrStartEnd"] + sjdbFileChrStartEnd
    if genomeTransformVCF:
        cmd += ["--genomeTransformVC", genomeTransformVCF]
    if genomeTransformType != "None":
        cmd += ["--genomeTransformType", genomeTransformType]
    if parametersFiles:
        cmd += ["--parametersFiles", parametersFiles]
    if sysShell:
        cmd += ["--sysShell", sysShell]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"STAR genomeGenerate failed with exit code {e.returncode}",
            "output_files": []
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": [str(genome_dir_path)]
    }


@mcp.tool()
def star_align_reads(
    genomeDir: str,
    readFilesIn: List[str],
    runThreadN: int = 1,
    readFilesCommand: Optional[str] = None,
    outFileNamePrefix: str = "./",
    outSAMtype: List[str] = ["SAM"],
    outSAMstrandField: Optional[str] = None,
    outSAMattributes: Optional[List[str]] = None,
    outFilterMultimapNmax: int = 10,
    outFilterMismatchNmax: int = 10,
    outFilterMismatchNoverLmax: float = 0.3,
    outFilterMismatchNoverReadLmax: float = 1.0,
    outFilterScoreMin: int = 0,
    outFilterScoreMinOverLread: float = 0.66,
    outFilterMatchNmin: int = 0,
    outFilterMatchNminOverLread: float = 0.66,
    outFilterIntronMotifs: Optional[str] = None,
    outSAMunmapped: Optional[str] = None,
    outSAMorder: str = "Paired",
    outSAMprimaryFlag: str = "OneBestScore",
    outSAMmapqUnique: int = 255,
    outSAMflagOR: int = 0,
    outSAMflagAND: int = 65535,
    outSAMattrRGline: Optional[str] = None,
    outSAMheaderHD: Optional[str] = None,
    outSAMheaderPG: Optional[str] = None,
    outSAMheaderCommentFile: Optional[str] = None,
    outSAMfilter: Optional[str] = None,
    outSAMmultNmax: int = -1,
    outSAMtlen: int = 1,
    outBAMcompression: int = 1,
    outBAMsortingThreadN: int = 0,
    outBAMsortingBinsN: int = 50,
    genomeLoad: str = "NoSharedMemory",
    sjdbGTFfile: Optional[str] = None,
    sjdbFileChrStartEnd: Optional[List[str]] = None,
    sjdbOverhang: int = 100,
    twopassMode: str = "None",
    twopass1readsN: int = -1,
    varVCFfile: Optional[str] = None,
    waspOutputMode: Optional[str] = None,
    quantMode: Optional[List[str]] = None,
    quantTranscriptomeBAMcompression: int = 1,
    quantTranscriptomeSAMoutput: Optional[str] = None,
    chimSegmentMin: int = 0,
    chimOutType: Optional[List[str]] = None,
    chimScoreMin: int = 0,
    chimScoreDropMax: int = 20,
    chimScoreSeparation: int = 10,
    chimScoreJunctionNonGTAG: int = -1,
    chimJunctionOverhangMin: int = 20,
    chimSegmentReadGapMax: int = 0,
    chimFilter: str = "banGenomicN",
    chimMainSegmentMultNmax: int = 10,
    chimMultimapNmax: int = 0,
    chimMultimapScoreRange: int = 1,
    chimNonchimScoreDropMin: int = 20,
    outReadsUnmapped: Optional[str] = None,
    readMapNumber: int = -1,
    readFilesType: str = "Fastx",
    readFilesSAMattrKeep: Optional[str] = None,
    readFilesPrefix: Optional[str] = None,
    readMatesLengthsIn: str = "NotEqual",
    readNameSeparator: str = "/",
    readQualityScoreBase: int = 33,
    soloType: Optional[List[str]] = None,
    soloCBwhitelist: Optional[List[str]] = None,
    soloFeatures: Optional[List[str]] = None,
    soloMultiMappers: Optional[str] = None,
    soloUMIdedup: Optional[str] = None,
    soloUMIfiltering: Optional[str] = None,
    soloOutFileNames: Optional[List[str]] = None,
    soloCellFilter: Optional[str] = None,
    soloStrand: str = "Forward",
    soloCBtype: str = "Sequence",
    soloCBstart: int = 1,
    soloCBlen: int = 16,
    soloUMIstart: int = 17,
    soloUMIlen: int = 10,
    soloBarcodeReadLength: int = 1,
    soloBarcodeMate: int = 0,
    soloCBposition: Optional[List[str]] = None,
    soloUMIposition: Optional[str] = None,
    soloAdapterSequence: Optional[str] = None,
    soloAdapterMismatchesNmax: int = 1,
    soloCBmatchWLtype: str = "1MM multi",
    soloInputSAMattrBarcodeSeq: Optional[List[str]] = None,
    soloInputSAMattrBarcodeQual: Optional[List[str]] = None,
):
    """
    Run STAR mapping of RNA-seq reads to a genome.
    Supports single-end and paired-end reads, compressed inputs, annotations on the fly,
    2-pass mapping, chimeric detection, WASP filtering, and STARsolo single-cell RNA-seq.
    """
    # Validate genomeDir
    genome_dir_path = Path(genomeDir)
    if not genome_dir_path.exists() or not genome_dir_path.is_dir():
        raise FileNotFoundError(f"Genome directory {genomeDir} does not exist or is not a directory.")

    # Validate read files
    if not readFilesIn or len(readFilesIn) == 0:
        raise ValueError("At least one read file must be provided in readFilesIn.")
    for rf in readFilesIn:
        rf_path = Path(rf)
        if not rf_path.is_file():
            raise FileNotFoundError(f"Read file {rf} does not exist.")

    # Validate outSAMtype values
    valid_outSAMtypes = {"SAM", "BAM", "Unsorted", "SortedByCoordinate", "None"}
    for ost in outSAMtype:
        if ost not in valid_outSAMtypes:
            raise ValueError(f"Invalid outSAMtype value: {ost}. Must be one of {valid_outSAMtypes}.")

    # Validate genomeLoad
    valid_genomeLoad = {"NoSharedMemory", "LoadAndKeep", "LoadAndRemove", "LoadAndExit", "Remove"}
    if genomeLoad not in valid_genomeLoad:
        raise ValueError(f"Invalid genomeLoad value: {genomeLoad}. Must be one of {valid_genomeLoad}.")

    # Validate twopassMode
    valid_twopassMode = {"None", "Basic"}
    if twopassMode not in valid_twopassMode:
        raise ValueError(f"Invalid twopassMode value: {twopassMode}. Must be one of {valid_twopassMode}.")

    # Validate outFilterIntronMotifs
    valid_outFilterIntronMotifs = {None, "None", "RemoveNoncanonical", "RemoveNoncanonicalUnannotated"}
    if outFilterIntronMotifs not in valid_outFilterIntronMotifs:
        raise ValueError(f"Invalid outFilterIntronMotifs value: {outFilterIntronMotifs}. Must be one of {valid_outFilterIntronMotifs}.")

    # Validate outSAMunmapped
    valid_outSAMunmapped = {None, "None", "Within", "KeepPairs"}
    if outSAMunmapped not in valid_outSAMunmapped:
        raise ValueError(f"Invalid outSAMunmapped value: {outSAMunmapped}. Must be one of {valid_outSAMunmapped}.")

    # Validate outReadsUnmapped
    valid_outReadsUnmapped = {None, "None", "Fastx"}
    if outReadsUnmapped not in valid_outReadsUnmapped:
        raise ValueError(f"Invalid outReadsUnmapped value: {outReadsUnmapped}. Must be one of {valid_outReadsUnmapped}.")

    # Build command
    cmd = ["STAR"]
    cmd += ["--runMode", "alignReads"]
    cmd += ["--genomeDir", str(genome_dir_path)]
    cmd += ["--runThreadN", str(runThreadN)]

    # readFilesIn can be multiple files, possibly paired-end separated by space
    # STAR expects read1 [read2], so join with space
    # For multiple files per mate, comma separated lists are used
    # We join all readFilesIn by space, assuming user provides correct order
    cmd += ["--readFilesIn"] + readFilesIn

    if readFilesCommand:
        cmd += ["--readFilesCommand", readFilesCommand]

    if outFileNamePrefix:
        cmd += ["--outFileNamePrefix", outFileNamePrefix]

    if outSAMtype:
        cmd += ["--outSAMtype"] + outSAMtype

    if outSAMstrandField:
        cmd += ["--outSAMstrandField", outSAMstrandField]

    if outSAMattributes:
        cmd += ["--outSAMattributes"] + outSAMattributes

    cmd += ["--outFilterMultimapNmax", str(outFilterMultimapNmax)]
    cmd += ["--outFilterMismatchNmax", str(outFilterMismatchNmax)]
    cmd += ["--outFilterMismatchNoverLmax", str(outFilterMismatchNoverLmax)]
    cmd += ["--outFilterMismatchNoverReadLmax", str(outFilterMismatchNoverReadLmax)]
    cmd += ["--outFilterScoreMin", str(outFilterScoreMin)]
    cmd += ["--outFilterScoreMinOverLread", str(outFilterScoreMinOverLread)]
    cmd += ["--outFilterMatchNmin", str(outFilterMatchNmin)]
    cmd += ["--outFilterMatchNminOverLread", str(outFilterMatchNminOverLread)]

    if outFilterIntronMotifs:
        cmd += ["--outFilterIntronMotifs", outFilterIntronMotifs]

    if outSAMunmapped:
        cmd += ["--outSAMunmapped", outSAMunmapped]

    cmd += ["--outSAMorder", outSAMorder]
    cmd += ["--outSAMprimaryFlag", outSAMprimaryFlag]
    cmd += ["--outSAMmapqUnique", str(outSAMmapqUnique)]
    cmd += ["--outSAMflagOR", str(outSAMflagOR)]
    cmd += ["--outSAMflagAND", str(outSAMflagAND)]

    if outSAMattrRGline:
        cmd += ["--outSAMattrRGline", outSAMattrRGline]
    if outSAMheaderHD:
        cmd += ["--outSAMheaderHD", outSAMheaderHD]
    if outSAMheaderPG:
        cmd += ["--outSAMheaderPG", outSAMheaderPG]
    if outSAMheaderCommentFile:
        cmd += ["--outSAMheaderCommentFile", outSAMheaderCommentFile]
    if outSAMfilter:
        cmd += ["--outSAMfilter", outSAMfilter]

    cmd += ["--outSAMmultNmax", str(outSAMmultNmax)]
    cmd += ["--outSAMtlen", str(outSAMtlen)]
    cmd += ["--outBAMcompression", str(outBAMcompression)]
    cmd += ["--outBAMsortingThreadN", str(outBAMsortingThreadN)]
    cmd += ["--outBAMsortingBinsN", str(outBAMsortingBinsN)]

    cmd += ["--genomeLoad", genomeLoad]

    if sjdbGTFfile:
        gtf_path = Path(sjdbGTFfile)
        if not gtf_path.is_file():
            raise FileNotFoundError(f"GTF file {sjdbGTFfile} does not exist.")
        cmd += ["--sjdbGTFfile", sjdbGTFfile]

    if sjdbFileChrStartEnd:
        for sjdb_file in sjdbFileChrStartEnd:
            sjdb_path = Path(sjdb_file)
            if not sjdb_path.is_file():
                raise FileNotFoundError(f"SJDB file {sjdb_file} does not exist.")
        cmd += ["--sjdbFileChrStartEnd"] + sjdbFileChrStartEnd

    cmd += ["--sjdbOverhang", str(sjdbOverhang)]

    cmd += ["--twopassMode", twopassMode]
    cmd += ["--twopass1readsN", str(twopass1readsN)]

    if varVCFfile:
        vcf_path = Path(varVCFfile)
        if not vcf_path.is_file():
            raise FileNotFoundError(f"VCF file {varVCFfile} does not exist.")
        cmd += ["--varVCFfile", varVCFfile]

    if waspOutputMode:
        cmd += ["--waspOutputMode", waspOutputMode]

    if quantMode:
        cmd += ["--quantMode"] + quantMode

    cmd += ["--quantTranscriptomeBAMcompression", str(quantTranscriptomeBAMcompression)]

    if quantTranscriptomeSAMoutput:
        cmd += ["--quantTranscriptomeSAMoutput", quantTranscriptomeSAMoutput]

    cmd += ["--chimSegmentMin", str(chimSegmentMin)]

    if chimOutType:
        cmd += ["--chimOutType"] + chimOutType

    cmd += ["--chimScoreMin", str(chimScoreMin)]
    cmd += ["--chimScoreDropMax", str(chimScoreDropMax)]
    cmd += ["--chimScoreSeparation", str(chimScoreSeparation)]
    cmd += ["--chimScoreJunctionNonGTAG", str(chimScoreJunctionNonGTAG)]
    cmd += ["--chimJunctionOverhangMin", str(chimJunctionOverhangMin)]
    cmd += ["--chimSegmentReadGapMax", str(chimSegmentReadGapMax)]
    cmd += ["--chimFilter", chimFilter]
    cmd += ["--chimMainSegmentMultNmax", str(chimMainSegmentMultNmax)]
    cmd += ["--chimMultimapNmax", str(chimMultimapNmax)]
    cmd += ["--chimMultimapScoreRange", str(chimMultimapScoreRange)]
    cmd += ["--chimNonchimScoreDropMin", str(chimNonchimScoreDropMin)]

    if outReadsUnmapped:
        cmd += ["--outReadsUnmapped", outReadsUnmapped]

    cmd += ["--readMapNumber", str(readMapNumber)]
    cmd += ["--readFilesType", readFilesType]

    if readFilesSAMattrKeep:
        cmd += ["--readFilesSAMattrKeep", readFilesSAMattrKeep]

    if readFilesPrefix:
        cmd += ["--readFilesPrefix", readFilesPrefix]

    cmd += ["--readMatesLengthsIn", readMatesLengthsIn]
    cmd += ["--readNameSeparator", readNameSeparator]
    cmd += ["--readQualityScoreBase", str(readQualityScoreBase)]

    # STARsolo options
    if soloType:
        cmd += ["--soloType"] + soloType
    if soloCBwhitelist:
        cmd += ["--soloCBwhitelist"] + soloCBwhitelist
    if soloFeatures:
        cmd += ["--soloFeatures"] + soloFeatures
    if soloMultiMappers:
        cmd += ["--soloMultiMappers", soloMultiMappers]
    if soloUMIdedup:
        cmd += ["--soloUMIdedup", soloUMIdedup]
    if soloUMIfiltering:
        cmd += ["--soloUMIfiltering", soloUMIfiltering]
    if soloOutFileNames:
        cmd += ["--soloOutFileNames"] + soloOutFileNames
    if soloCellFilter:
        cmd += ["--soloCellFilter", soloCellFilter]
    if soloStrand:
        cmd += ["--soloStrand", soloStrand]
    if soloCBtype:
        cmd += ["--soloCBtype", soloCBtype]
    cmd += ["--soloCBstart", str(soloCBstart)]
    cmd += ["--soloCBlen", str(soloCBlen)]
    cmd += ["--soloUMIstart", str(soloUMIstart)]
    cmd += ["--soloUMIlen", str(soloUMIlen)]
    cmd += ["--soloBarcodeReadLength", str(soloBarcodeReadLength)]
    cmd += ["--soloBarcodeMate", str(soloBarcodeMate)]
    if soloCBposition:
        cmd += ["--soloCBposition"] + soloCBposition
    if soloUMIposition:
        cmd += ["--soloUMIposition", soloUMIposition]
    if soloAdapterSequence:
        cmd += ["--soloAdapterSequence", soloAdapterSequence]
    cmd += ["--soloAdapterMismatchesNmax", str(soloAdapterMismatchesNmax)]
    cmd += ["--soloCBmatchWLtype", soloCBmatchWLtype]
    if soloInputSAMattrBarcodeSeq:
        cmd += ["--soloInputSAMattrBarcodeSeq"] + soloInputSAMattrBarcodeSeq
    if soloInputSAMattrBarcodeQual:
        cmd += ["--soloInputSAMattrBarcodeQual"] + soloInputSAMattrBarcodeQual

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "error": f"STAR alignReads failed with exit code {e.returncode}",
            "output_files": []
        }

    # Collect output files based on outFileNamePrefix and outSAMtype
    output_files = []
    prefix_path = Path(outFileNamePrefix)
    # STAR default output files if prefix is "./"
    # Aligned.out.sam or BAM files
    if "SAM" in outSAMtype:
        sam_file = prefix_path / "Aligned.out.sam"
        output_files.append(str(sam_file))
    if "BAM" in outSAMtype:
        # BAM files can be Unsorted or SortedByCoordinate or both
        if "Unsorted" in outSAMtype:
            bam_unsorted = prefix_path / "Aligned.out.bam"
            output_files.append(str(bam_unsorted))
        if "SortedByCoordinate" in outSAMtype:
            bam_sorted = prefix_path / "Aligned.sortedByCoord.out.bam"
            output_files.append(str(bam_sorted))
    # Log files
    log_out = prefix_path / "Log.out"
    log_final_out = prefix_path / "Log.final.out"
    log_progress_out = prefix_path / "Log.progress.out"
    output_files.extend([str(log_out), str(log_final_out), str(log_progress_out)])

    # SJ.out.tab splice junctions file
    sj_out = prefix_path / "SJ.out.tab"
    output_files.append(str(sj_out))

    # Unmapped reads files if requested
    if outReadsUnmapped == "Fastx":
        unmapped_mate1 = prefix_path / "Unmapped.out.mate1"
        unmapped_mate2 = prefix_path / "Unmapped.out.mate2"
        output_files.extend([str(unmapped_mate1), str(unmapped_mate2)])

    return {
        "command_executed": " ".join(cmd),
        "stdout": stdout,
        "stderr": stderr,
        "output_files": output_files
    }


if __name__ == '__main__':
    mcp.run()