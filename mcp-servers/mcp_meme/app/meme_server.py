from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List, Literal
import subprocess

mcp = FastMCP()


@mcp.tool()
def meme(
    primary_sequence_file: str,
    output_dir: str = "meme_out",
    output_dir_overwrite: Optional[str] = None,
    text_output: bool = False,
    brief: int = 1000,
    objfun: Literal["classic", "de", "se", "cd", "ce", "nc"] = "classic",
    test: Literal["mhg", "mbn", "mrs"] = "mhg",
    use_llr: bool = False,
    neg_control_file: Optional[str] = None,
    shuf_kmer: int = 2,
    hsfrac: float = 0.5,
    cefrac: float = 0.25,
    searchsize: int = 100000,
    norand: bool = False,
    csites: int = 1000,
    seed: int = 0,
    alph_file: Optional[str] = None,
    dna: bool = False,
    rna: bool = False,
    protein: bool = False,
    revcomp: bool = False,
    pal: bool = False,
    mod: Literal["oops", "zoops", "anr"] = "zoops",
    nmotifs: int = 1,
    evt: float = 10.0,
    time_limit: Optional[int] = None,
    nsites: Optional[int] = None,
    minsites: int = 2,
    maxsites: Optional[int] = None,
    wn_sites: float = 0.8,
    w: Optional[int] = None,
    minw: int = 8,
    maxw: int = 50,
    allw: bool = False,
    nomatrim: bool = False,
    wg: int = 11,
    ws: int = 1,
    noendgaps: bool = False,
    bfile: Optional[str] = None,
    markov_order: int = 0,
    psp_file: Optional[str] = None,
    maxiter: int = 50,
    distance: float = 0.001,
    prior: Literal["dirichlet", "dmix", "mega", "megap", "addone"] = "dirichlet",
    b: float = 0.01,
    plib: Optional[str] = None,
    spfuzz: Optional[float] = None,
    spmap: Literal["uni", "pam"] = "uni",
    cons: Optional[List[str]] = None,
    np: Optional[str] = None,
    maxsize: int = 0,
    nostatus: bool = False,
    sf: bool = False,
    verbose: bool = False,
) -> dict:
    """
    Run MEME motif discovery on a primary sequence file.

    Args:
        primary_sequence_file: Path to primary sequences in FASTA format or 'stdin'.
        output_dir: Directory to create for output files (incompatible with output_dir_overwrite).
        output_dir_overwrite: Directory to create or overwrite for output files (incompatible with output_dir).
        text_output: Output text format only to stdout.
        brief: Reduce output size if more than this many sequences (default 1000).
        objfun: Objective function to use for motif selection.
        test: Statistical test for motif enrichment (only for de or se objfun).
        use_llr: Use log-likelihood ratio method for evaluating EM starting points (only for classic objfun).
        neg_control_file: Control sequences file in FASTA format.
        shuf_kmer: k-mer size for shuffling primary sequences to create control sequences (1-6).
        hsfrac: Fraction of primary sequences held out for estimating motif parameters.
        cefrac: Fraction of sequence length defining central region for central enrichment.
        searchsize: Max letters used in motif search; 0 means no limit.
        norand: Do not randomize input sequence order before sampling.
        csites: Max number of sites used for computing E-value with classic objfun.
        seed: Random seed for shuffling and sampling.
        alph_file: Alphabet definition file (incompatible with dna, rna, protein).
        dna: Use standard DNA alphabet.
        rna: Use standard RNA alphabet.
        protein: Use standard protein alphabet.
        revcomp: Consider both strands for complementable alphabets.
        pal: Only look for palindromes in complementable alphabets.
        mod: Motif site distribution model.
        nmotifs: Number of motifs to find.
        evt: Stop if last motif E-value > evt.
        time_limit: Stop if estimated run time exceeds this (seconds).
        nsites: Exact number of motif occurrences (overrides minsites and maxsites).
        minsites: Minimum number of motif occurrences.
        maxsites: Maximum number of motif occurrences.
        wn_sites: Weight bias towards motifs with expected number of sites [0..1).
        w: Exact motif width.
        minw: Minimum motif width.
        maxw: Maximum motif width.
        allw: Find starting points for all widths from minw to maxw.
        nomatrim: Do not trim motif width using multiple alignments.
        wg: Gap opening cost for motif trimming.
        ws: Gap extension cost for motif trimming.
        noendgaps: Do not count end gaps in motif trimming.
        bfile: Markov background model file.
        markov_order: Maximum order of Markov model to read or create.
        psp_file: Position-specific priors file.
        maxiter: Maximum EM iterations per starting point.
        distance: EM convergence threshold.
        prior: Type of prior to use.
        b: Strength of prior on model parameters.
        plib: Dirichlet mixtures prior library file.
        spfuzz: Fuzziness parameter for sequence to theta mapping.
        spmap: Mapping function for estimating theta.
        cons: List of consensus sequences to override starting points.
        np: Number of processors or MPI command string.
        maxsize: Maximum allowed dataset size in letters (0 means no limit).
        nostatus: Suppress status messages.
        sf: Print sequence file name as given.
        verbose: Print extensive status messages.

    Returns:
        dict: Contains command_executed, stdout, stderr, and output_files list.
    """
    # Validate input file
    if primary_sequence_file != "stdin":
        primary_path = Path(primary_sequence_file)
        if not primary_path.is_file():
            raise FileNotFoundError(f"Primary sequence file not found: {primary_sequence_file}")

    # Validate mutually exclusive output directory options
    if output_dir and output_dir_overwrite:
        raise ValueError("Options output_dir (-o) and output_dir_overwrite (-oc) are mutually exclusive.")

    # Validate shuf_kmer range
    if not (1 <= shuf_kmer <= 6):
        raise ValueError("shuf_kmer must be between 1 and 6.")

    # Validate wn_sites range
    if not (0 <= wn_sites < 1):
        raise ValueError("wn_sites must be in the range [0..1).")

    # Validate prior option
    if prior not in {"dirichlet", "dmix", "mega", "megap", "addone"}:
        raise ValueError("Invalid prior option.")

    # Validate objfun and test compatibility
    if objfun not in {"classic", "de", "se", "cd", "ce", "nc"}:
        raise ValueError("Invalid objfun option.")
    if objfun not in {"de", "se"} and test != "mhg":
        raise ValueError("Option -test only valid with objfun 'de' or 'se'.")

    # Validate alphabet options exclusivity
    alph_opts = sum([bool(alph_file), dna, rna, protein])
    if alph_opts > 1:
        raise ValueError("Only one of alph_file, dna, rna, protein options can be specified.")

    # Validate motif width options
    if w is not None:
        if w < 1:
            raise ValueError("Motif width (-w) must be positive.")
        if w < minw or w > maxw:
            raise ValueError("Motif width (-w) must be between minw and maxw.")

    # Validate maxsites if given
    if maxsites is not None and maxsites < 1:
        raise ValueError("maxsites must be positive if specified.")

    # Validate evt positive
    if evt <= 0:
        raise ValueError("evt must be positive.")

    # Validate maxiter positive
    if maxiter < 1:
        raise ValueError("maxiter must be positive.")

    # Validate distance positive
    if distance <= 0:
        raise ValueError("distance must be positive.")

    # Validate spmap
    if spmap not in {"uni", "pam"}:
        raise ValueError("spmap must be 'uni' or 'pam'.")

    # Validate cons list if given
    if cons is not None:
        if not isinstance(cons, list):
            raise ValueError("cons must be a list of consensus sequences.")
        for c in cons:
            if not isinstance(c, str):
                raise ValueError("Each consensus sequence must be a string.")

    # Build command line
    cmd = ["meme"]

    # Primary sequence file
    if primary_sequence_file == "stdin":
        cmd.append("-")
    else:
        cmd.append(str(primary_sequence_file))

    # Output directory options
    if output_dir_overwrite:
        cmd.extend(["-oc", output_dir_overwrite])
    else:
        cmd.extend(["-o", output_dir])

    # Text output
    if text_output:
        cmd.append("-text")

    # Brief
    if brief != 1000:
        cmd.extend(["-brief", str(brief)])

    # Objective function
    if objfun != "classic":
        cmd.extend(["-objfun", objfun])

    # Test (only for de or se)
    if objfun in {"de", "se"} and test != "mhg":
        cmd.extend(["-test", test])

    # Use LLR
    if use_llr:
        cmd.append("-use_llr")

    # Control sequences
    if neg_control_file:
        neg_path = Path(neg_control_file)
        if not neg_path.is_file():
            raise FileNotFoundError(f"Control sequence file not found: {neg_control_file}")
        cmd.extend(["-neg", neg_control_file])

    # Shuffle kmer
    if shuf_kmer != 2:
        cmd.extend(["-shuf", str(shuf_kmer)])

    # hsfrac
    if hsfrac != 0.5:
        cmd.extend(["-hsfrac", str(hsfrac)])

    # cefrac
    if cefrac != 0.25:
        cmd.extend(["-cefrac", str(cefrac)])

    # searchsize
    if searchsize != 100000:
        cmd.extend(["-searchsize", str(searchsize)])

    # norand
    if norand:
        cmd.append("-norand")

    # csites
    if csites != 1000:
        cmd.extend(["-csites", str(csites)])

    # seed
    if seed != 0:
        cmd.extend(["-seed", str(seed)])

    # Alphabet options
    if alph_file:
        alph_path = Path(alph_file)
        if not alph_path.is_file():
            raise FileNotFoundError(f"Alphabet file not found: {alph_file}")
        cmd.extend(["-alph", alph_file])
    else:
        if dna:
            cmd.append("-dna")
        elif rna:
            cmd.append("-rna")
        elif protein:
            cmd.append("-protein")

    # Strands & palindromes
    if revcomp:
        cmd.append("-revcomp")
    if pal:
        cmd.append("-pal")

    # Motif site distribution model
    if mod != "zoops":
        cmd.extend(["-mod", mod])

    # Number of motifs
    if nmotifs != 1:
        cmd.extend(["-nmotifs", str(nmotifs)])

    # evt
    if evt != 10.0:
        cmd.extend(["-evt", str(evt)])

    # time limit
    if time_limit is not None:
        if time_limit < 1:
            raise ValueError("time_limit must be positive if specified.")
        cmd.extend(["-time", str(time_limit)])

    # nsites, minsites, maxsites
    if nsites is not None:
        if nsites < 1:
            raise ValueError("nsites must be positive if specified.")
        cmd.extend(["-nsites", str(nsites)])
    else:
        if minsites != 2:
            cmd.extend(["-minsites", str(minsites)])
        if maxsites is not None:
            cmd.extend(["-maxsites", str(maxsites)])

    # wn_sites
    if wn_sites != 0.8:
        cmd.extend(["-wnsites", str(wn_sites)])

    # Motif width options
    if w is not None:
        cmd.extend(["-w", str(w)])
    else:
        if minw != 8:
            cmd.extend(["-minw", str(minw)])
        if maxw != 50:
            cmd.extend(["-maxw", str(maxw)])

    # allw
    if allw:
        cmd.append("-allw")

    # nomatrim
    if nomatrim:
        cmd.append("-nomatrim")

    # wg, ws, noendgaps
    if wg != 11:
        cmd.extend(["-wg", str(wg)])
    if ws != 1:
        cmd.extend(["-ws", str(ws)])
    if noendgaps:
        cmd.append("-noendgaps")

    # Background model
    if bfile:
        bfile_path = Path(bfile)
        if not bfile_path.is_file():
            raise FileNotFoundError(f"Background model file not found: {bfile}")
        cmd.extend(["-bfile", bfile])
    if markov_order != 0:
        cmd.extend(["-markov_order", str(markov_order)])

    # Position-specific priors
    if psp_file:
        psp_path = Path(psp_file)
        if not psp_path.is_file():
            raise FileNotFoundError(f"Position-specific priors file not found: {psp_file}")
        cmd.extend(["-psp", psp_file])

    # EM algorithm
    if maxiter != 50:
        cmd.extend(["-maxiter", str(maxiter)])
    if distance != 0.001:
        cmd.extend(["-distance", str(distance)])

    # Prior
    if prior != "dirichlet":
        cmd.extend(["-prior", prior])
    if b != 0.01:
        cmd.extend(["-b", str(b)])

    # Dirichlet mixtures prior library
    if plib:
        plib_path = Path(plib)
        if not plib_path.is_file():
            raise FileNotFoundError(f"Dirichlet mixtures prior library file not found: {plib}")
        cmd.extend(["-plib", plib])

    # spfuzz
    if spfuzz is not None:
        if spfuzz < 0:
            raise ValueError("spfuzz must be non-negative if specified.")
        cmd.extend(["-spfuzz", str(spfuzz)])

    # spmap
    if spmap != "uni":
        cmd.extend(["-spmap", spmap])

    # Consensus sequences
    if cons:
        for cseq in cons:
            cmd.extend(["-cons", cseq])

    # Parallel processors
    if np:
        cmd.extend(["-p", np])

    # maxsize
    if maxsize != 0:
        cmd.extend(["-maxsize", str(maxsize)])

    # nostatus
    if nostatus:
        cmd.append("-nostatus")

    # sf
    if sf:
        cmd.append("-sf")

    # verbose
    if verbose:
        cmd.append("-V")

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
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": f"MEME execution failed with return code {e.returncode}",
        }

    # Determine output directory path
    out_dir_path = Path(output_dir_overwrite if output_dir_overwrite else output_dir)

    # Collect output files if output directory exists
    output_files = []
    if out_dir_path.is_dir():
        # Collect known output files
        known_files = [
            "meme.html",
            "meme.txt",
            "meme.xml",
        ]
        # Add logo files (logoN.png, logoN.eps, logo_rcN.png, logo_rcN.eps)
        # We will glob for logo*.png and logo*.eps files
        output_files.extend([str(p) for p in out_dir_path.glob("logo*.png")])
        output_files.extend([str(p) for p in out_dir_path.glob("logo*.eps")])
        # Add known files if exist
        for fname in known_files:
            fpath = out_dir_path / fname
            if fpath.is_file():
                output_files.append(str(fpath))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


if __name__ == '__main__':
    mcp.run()