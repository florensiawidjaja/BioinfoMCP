from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()

def _validate_file_path(path: Path, must_exist: bool = True):
    if must_exist and not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    return path.resolve()

def _validate_positive_int(value: int, name: str):
    if value < 0:
        raise ValueError(f"{name} must be non-negative integer, got {value}")
    return value

def _validate_positive_float(value: float, name: str):
    if value < 0:
        raise ValueError(f"{name} must be non-negative float, got {value}")
    return value

@mcp.tool()
def mafft(
    input: Path,
    auto: bool = False,
    sixmerpair: bool = True,
    globalpair: bool = False,
    localpair: bool = False,
    genafpair: bool = False,
    fastapair: bool = False,
    weighti: float = 2.7,
    retree: int = 2,
    maxiterate: int = 0,
    fft: bool = True,
    nofft: bool = False,
    noscore: bool = False,
    memsave: bool = False,
    parttree: bool = False,
    dpparttree: bool = False,
    fastaparttree: bool = False,
    partsize: int = 50,
    groupsize: Optional[int] = None,
    op: float = 1.53,
    ep: float = 0.123,
    lop: float = -2.00,
    lep: float = 0.1,
    lexp: float = -0.1,
    LOP: float = -6.00,
    LEXP: float = 0.00,
    bl: int = 62,
    jtt: Optional[int] = None,
    tm: Optional[int] = None,
    aamatrix: Optional[Path] = None,
    fmodel: bool = False,
    clustalout: bool = False,
    inputorder: bool = True,
    reorder: bool = False,
    treeout: bool = False,
    quiet: bool = False,
    nuc: bool = False,
    amino: bool = False,
    seed: Optional[List[Path]] = None,
):
    """
    Run the main MAFFT multiple sequence alignment with flexible options.

    Parameters:
    - input: Input FASTA file with sequences to align.
    - auto: Automatically select strategy (L-INS-i, FFT-NS-i, FFT-NS-2).
    - sixmerpair: Use 6mer distance calculation (default True).
    - globalpair: Use Needleman-Wunsch global pairwise alignment.
    - localpair: Use Smith-Waterman local pairwise alignment.
    - genafpair: Use generalized affine gap local pairwise alignment.
    - fastapair: Use FASTA pairwise alignment (requires fasta34).
    - weighti: Weighting factor for consistency term (default 2.7).
    - retree: Number of guide tree rebuilds (default 2).
    - maxiterate: Number of iterative refinement cycles (default 0).
    - fft: Use FFT approximation in group-to-group alignment (default True).
    - nofft: Disable FFT approximation (default False).
    - noscore: Do not check alignment score in refinement (default False).
    - memsave: Use Myers-Miller algorithm for memory saving.
    - parttree: Use PartTree fast tree-building method.
    - dpparttree: Use DP-based PartTree method.
    - fastaparttree: Use FASTA-based PartTree method.
    - partsize: Number of partitions in PartTree (default 50).
    - groupsize: Max number of sequences per alignment group.
    - op: Gap opening penalty at group-to-group alignment (default 1.53).
    - ep: Offset value for group-to-group alignment (default 0.123).
    - lop: Gap opening penalty at local pairwise alignment (default -2.00).
    - lep: Offset value at local pairwise alignment (default 0.1).
    - lexp: Gap extension penalty at local pairwise alignment (default -0.1).
    - LOP: Gap opening penalty to skip alignment (default -6.00).
    - LEXP: Gap extension penalty to skip alignment (default 0.00).
    - bl: BLOSUM matrix number (30,45,62,80; default 62).
    - jtt: JTT PAM number matrix (overrides bl if set).
    - tm: Transmembrane PAM number matrix (overrides bl if set).
    - aamatrix: Path to user-defined amino acid scoring matrix file.
    - fmodel: Incorporate AA/nuc composition info into scoring matrix.
    - clustalout: Output alignment in Clustal format instead of FASTA.
    - inputorder: Output order same as input (default True).
    - reorder: Output order aligned (overrides inputorder).
    - treeout: Output guide tree to input.tree file.
    - quiet: Suppress progress reporting.
    - nuc: Assume nucleotide sequences.
    - amino: Assume amino acid sequences.
    - seed: List of seed alignment files to preserve alignment within seeds.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)
    if aamatrix is not None:
        aamatrix = _validate_file_path(aamatrix)

    if seed:
        seed = [_validate_file_path(s) for s in seed]

    # Validate numeric parameters
    weighti = float(weighti)
    retree = _validate_positive_int(retree, "retree")
    maxiterate = _validate_positive_int(maxiterate, "maxiterate")
    partsize = _validate_positive_int(partsize, "partsize")
    if groupsize is not None:
        groupsize = _validate_positive_int(groupsize, "groupsize")

    # Validate bl parameter
    if bl not in (30, 45, 62, 80):
        raise ValueError("bl must be one of 30, 45, 62, or 80")

    # Validate mutually exclusive amino/nuc
    if amino and nuc:
        raise ValueError("Cannot specify both --amino and --nuc")

    # Build command line
    cmd = ["mafft"]

    # Algorithm options
    if auto:
        cmd.append("--auto")
    if sixmerpair:
        cmd.append("--6merpair")
    if globalpair:
        cmd.append("--globalpair")
    if localpair:
        cmd.append("--localpair")
    if genafpair:
        cmd.append("--genafpair")
    if fastapair:
        cmd.append("--fastapair")

    # Parameters with values
    cmd.extend(["--weighti", str(weighti)])
    cmd.extend(["--retree", str(retree)])
    cmd.extend(["--maxiterate", str(maxiterate)])

    if fft and not nofft:
        cmd.append("--fft")
    if nofft:
        cmd.append("--nofft")
    if noscore:
        cmd.append("--noscore")
    if memsave:
        cmd.append("--memsave")
    if parttree:
        cmd.append("--parttree")
    if dpparttree:
        cmd.append("--dpparttree")
    if fastaparttree:
        cmd.append("--fastaparttree")

    cmd.extend(["--partsize", str(partsize)])
    if groupsize is not None:
        cmd.extend(["--groupsize", str(groupsize)])

    cmd.extend(["--op", str(op)])
    cmd.extend(["--ep", str(ep)])
    cmd.extend(["--lop", str(lop)])
    cmd.extend(["--lep", str(lep)])
    cmd.extend(["--lexp", str(lexp)])
    cmd.extend(["--LOP", str(LOP)])
    cmd.extend(["--LEXP", str(LEXP)])

    # Matrix selection
    if jtt is not None:
        if jtt <= 0:
            raise ValueError("jtt must be > 0")
        cmd.extend(["--jtt", str(jtt)])
    elif tm is not None:
        if tm <= 0:
            raise ValueError("tm must be > 0")
        cmd.extend(["--tm", str(tm)])
    else:
        cmd.extend(["--bl", str(bl)])

    if aamatrix is not None:
        cmd.extend(["--aamatrix", str(aamatrix)])

    if fmodel:
        cmd.append("--fmodel")

    # Output options
    if clustalout:
        cmd.append("--clustalout")
    if inputorder and not reorder:
        cmd.append("--inputorder")
    if reorder:
        cmd.append("--reorder")
    if treeout:
        cmd.append("--treeout")
    if quiet:
        cmd.append("--quiet")

    # Input type
    if nuc:
        cmd.append("--nuc")
    if amino:
        cmd.append("--amino")

    # Seed alignments
    if seed:
        for s in seed:
            cmd.extend(["--seed", str(s)])

    cmd.append(str(input))

    # Output file: we will create a temporary output file
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        # Run command redirecting stdout to output_path
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        # Read stderr
        stderr = proc.stderr
        stdout = ""  # stdout redirected to file

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        # Clean up output file if error
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"MAFFT failed with exit code {e.returncode}",
        }

@mcp.tool()
def linsi(
    input: Path,
):
    """
    Run MAFFT L-INS-i method (accuracy-oriented, local pairwise alignment with iterative refinement).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["linsi", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"LINSI failed with exit code {e.returncode}",
        }

@mcp.tool()
def ginsi(
    input: Path,
):
    """
    Run MAFFT G-INS-i method (accuracy-oriented, global pairwise alignment with iterative refinement).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["ginsi", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"GINSI failed with exit code {e.returncode}",
        }

@mcp.tool()
def einsi(
    input: Path,
):
    """
    Run MAFFT E-INS-i method (accuracy-oriented, suitable for sequences with large unalignable regions).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["einsi", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"EINSI failed with exit code {e.returncode}",
        }

@mcp.tool()
def fftnsi(
    input: Path,
):
    """
    Run MAFFT FFT-NS-i method (speed-oriented, iterative refinement with two cycles).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["fftnsi", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"FFTNSI failed with exit code {e.returncode}",
        }

@mcp.tool()
def fftns(
    input: Path,
):
    """
    Run MAFFT FFT-NS-2 method (speed-oriented, fast progressive method).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["fftns", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"FFTNS failed with exit code {e.returncode}",
        }

@mcp.tool()
def nwnsi(
    input: Path,
):
    """
    Run MAFFT NW-NS-i method (speed-oriented, iterative refinement without FFT approximation, two cycles).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["nwnsi", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"NWNSI failed with exit code {e.returncode}",
        }

@mcp.tool()
def nwns(
    input: Path,
):
    """
    Run MAFFT NW-NS-2 method (speed-oriented, fast progressive method without FFT approximation).

    Parameters:
    - input: Input FASTA file with sequences to align.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    input = _validate_file_path(input)

    cmd = ["nwns", str(input)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"NWNS failed with exit code {e.returncode}",
        }

@mcp.tool()
def mafft_profile(
    group1: Path,
    group2: Path,
):
    """
    Run MAFFT profile alignment: align two groups of sequences (FASTA format).

    Parameters:
    - group1: FASTA file of first group.
    - group2: FASTA file of second group.

    Returns:
    Dict with command executed, stdout, stderr, and output file path.
    """
    group1 = _validate_file_path(group1)
    group2 = _validate_file_path(group2)

    cmd = ["mafft-profile", str(group1), str(group2)]

    import tempfile
    import os

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_out:
        output_path = Path(tmp_out.name)

    try:
        with open(output_path, "w") as out_fh:
            proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True, check=True)
        stderr = proc.stderr
        stdout = ""

        return {
            "command_executed": " ".join(cmd),
            "stdout": stdout,
            "stderr": stderr,
            "output_files": [str(output_path)],
        }
    except subprocess.CalledProcessError as e:
        if output_path.exists():
            os.unlink(output_path)
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout if e.stdout else "",
            "stderr": e.stderr if e.stderr else str(e),
            "output_files": [],
            "error": True,
            "error_message": f"mafft-profile failed with exit code {e.returncode}",
        }

if __name__ == '__main__':
    mcp.run()