from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import subprocess

mcp = FastMCP()


@mcp.tool()
def bwa_index(
    in_db_fasta: Path,
    p: Optional[str] = None,
    a: str = "is",
):
    """
    Index database sequences in the FASTA format using bwa index.
    -p STR: Prefix of the output database [default: same as db filename]
    -a STR: Algorithm for constructing BWT index. Options: 'is' (default), 'bwtsw'.
    """
    if not in_db_fasta.exists():
        raise FileNotFoundError(f"Input fasta file {in_db_fasta} does not exist")
    if a not in ("is", "bwtsw"):
        raise ValueError("Parameter 'a' must be either 'is' or 'bwtsw'")

    cmd = ["bwa", "index"]
    if p:
        cmd += ["-p", p]
    cmd += ["-a", a]
    cmd.append(str(in_db_fasta))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_files = []
        prefix = p if p else in_db_fasta.with_suffix("").name
        # BWA index creates multiple files with extensions: .amb, .ann, .bwt, .pac, .sa
        for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
            f = Path(prefix + ext)
            if f.exists():
                output_files.append(str(f.resolve()))
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
            "error": f"bwa index failed with return code {e.returncode}",
        }


@mcp.tool()
def bwa_mem(
    db_prefix: Path,
    reads_fq: Path,
    mates_fq: Optional[Path] = None,
    a: bool = False,
    C: bool = False,
    H: bool = False,
    M: bool = False,
    p: bool = False,
    t: int = 1,
    k: int = 19,
    w: int = 100,
    d: int = 100,
    r: float = 1.5,
    c: int = 10000,
    A: int = 1,
    B: int = 4,
    O: int = 6,
    E: int = 1,
    L: int = 5,
    U: int = 9,
    R: Optional[str] = None,
    v: int = 3,
    T: int = 30,
):
    """
    Align 70bp-1Mbp query sequences with the BWA-MEM algorithm.
    Supports single-end, paired-end, and interleaved paired-end reads.
    Parameters correspond to bwa mem options.
    """
    if not db_prefix.exists():
        raise FileNotFoundError(f"Database prefix {db_prefix} does not exist")
    if not reads_fq.exists():
        raise FileNotFoundError(f"Reads file {reads_fq} does not exist")
    if mates_fq and not mates_fq.exists():
        raise FileNotFoundError(f"Mates file {mates_fq} does not exist")
    if t < 1:
        raise ValueError("Number of threads 't' must be >= 1")
    if k < 1:
        raise ValueError("Minimum seed length 'k' must be >= 1")
    if w < 1:
        raise ValueError("Band width 'w' must be >= 1")
    if d < 0:
        raise ValueError("Off-diagonal X-dropoff 'd' must be >= 0")
    if r <= 0:
        raise ValueError("Trigger re-seeding ratio 'r' must be > 0")
    if c < 0:
        raise ValueError("Discard MEM occurrence 'c' must be >= 0")
    if A < 0 or B < 0 or O < 0 or E < 0 or L < 0 or U < 0:
        raise ValueError("Scoring penalties must be non-negative")
    if v < 0:
        raise ValueError("Verbose level 'v' must be >= 0")
    if T < 0:
        raise ValueError("Minimum output alignment score 'T' must be >= 0")

    cmd = ["bwa", "mem"]
    if a:
        cmd.append("-a")
    if C:
        cmd.append("-C")
    if H:
        cmd.append("-H")
    if M:
        cmd.append("-M")
    if p:
        cmd.append("-p")
    cmd += ["-t", str(t)]
    cmd += ["-k", str(k)]
    cmd += ["-w", str(w)]
    cmd += ["-d", str(d)]
    cmd += ["-r", str(r)]
    cmd += ["-c", str(c)]
    cmd += ["-A", str(A)]
    cmd += ["-B", str(B)]
    cmd += ["-O", str(O)]
    cmd += ["-E", str(E)]
    cmd += ["-L", str(L)]
    cmd += ["-U", str(U)]
    if R:
        # Replace literal \t with tab character
        R_fixed = R.replace("\\t", "\t")
        cmd += ["-R", R_fixed]
    cmd += ["-v", str(v)]
    cmd += ["-T", str(T)]
    cmd.append(str(db_prefix))
    cmd.append(str(reads_fq))
    if mates_fq and not p:
        cmd.append(str(mates_fq))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # bwa mem outputs SAM to stdout
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
            "error": f"bwa mem failed with return code {e.returncode}",
        }


@mcp.tool()
def bwa_aln(
    in_db_fasta: Path,
    in_query_fq: Path,
    n: float = 0.04,
    o: int = 1,
    e: int = -1,
    d: int = 16,
    i: int = 5,
    l: Optional[int] = None,
    k: int = 2,
    t: int = 1,
    M: int = 3,
    O: int = 11,
    E: int = 4,
    R: int = 0,
    c: bool = False,
    N: bool = False,
    q: int = 0,
    I: bool = False,
    B: int = 0,
    b: bool = False,
    zero: bool = False,
    one: bool = False,
    two: bool = False,
):
    """
    Find the SA coordinates of the input reads using bwa aln (BWA-backtrack).
    Parameters correspond to bwa aln options.
    """
    if not in_db_fasta.exists():
        raise FileNotFoundError(f"Input fasta file {in_db_fasta} does not exist")
    if not in_query_fq.exists():
        raise FileNotFoundError(f"Input query file {in_query_fq} does not exist")
    if n < 0:
        raise ValueError("Maximum edit distance 'n' must be non-negative")
    if o < 0:
        raise ValueError("Maximum number of gap opens 'o' must be non-negative")
    if e < -1:
        raise ValueError("Maximum number of gap extensions 'e' must be >= -1")
    if d < 0:
        raise ValueError("Disallow long deletion 'd' must be non-negative")
    if i < 0:
        raise ValueError("Disallow indel near ends 'i' must be non-negative")
    if l is not None and l < 1:
        raise ValueError("Seed length 'l' must be positive or None")
    if k < 0:
        raise ValueError("Maximum edit distance in seed 'k' must be non-negative")
    if t < 1:
        raise ValueError("Number of threads 't' must be >= 1")
    if M < 0 or O < 0 or E < 0 or R < 0 or q < 0 or B < 0:
        raise ValueError("Penalty and threshold parameters must be non-negative")

    cmd = ["bwa", "aln"]
    cmd += ["-n", str(n)]
    cmd += ["-o", str(o)]
    cmd += ["-e", str(e)]
    cmd += ["-d", str(d)]
    cmd += ["-i", str(i)]
    if l is not None:
        cmd += ["-l", str(l)]
    cmd += ["-k", str(k)]
    cmd += ["-t", str(t)]
    cmd += ["-M", str(M)]
    cmd += ["-O", str(O)]
    cmd += ["-E", str(E)]
    cmd += ["-R", str(R)]
    if c:
        cmd.append("-c")
    if N:
        cmd.append("-N")
    cmd += ["-q", str(q)]
    if I:
        cmd.append("-I")
    if B > 0:
        cmd += ["-B", str(B)]
    if b:
        cmd.append("-b")
    if zero:
        cmd.append("-0")
    if one:
        cmd.append("-1")
    if two:
        cmd.append("-2")
    cmd.append(str(in_db_fasta))
    cmd.append(str(in_query_fq))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # bwa aln outputs .sai to stdout
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
            "error": f"bwa aln failed with return code {e.returncode}",
        }


@mcp.tool()
def bwa_samse(
    in_db_fasta: Path,
    in_sai: Path,
    in_fq: Path,
    n: int = 3,
    r: Optional[str] = None,
):
    """
    Generate alignments in the SAM format given single-end reads using bwa samse.
    -n INT: Maximum number of alignments to output in XA tag [3]
    -r STR: Specify the read group header line (e.g. '@RG\\tID:foo\\tSM:bar')
    """
    if not in_db_fasta.exists():
        raise FileNotFoundError(f"Input fasta file {in_db_fasta} does not exist")
    if not in_sai.exists():
        raise FileNotFoundError(f"Input sai file {in_sai} does not exist")
    if not in_fq.exists():
        raise FileNotFoundError(f"Input fastq file {in_fq} does not exist")
    if n < 0:
        raise ValueError("Maximum number of alignments 'n' must be non-negative")

    cmd = ["bwa", "samse"]
    cmd += ["-n", str(n)]
    if r:
        r_fixed = r.replace("\\t", "\t")
        cmd += ["-r", r_fixed]
    cmd.append(str(in_db_fasta))
    cmd.append(str(in_sai))
    cmd.append(str(in_fq))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # bwa samse outputs SAM to stdout
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
            "error": f"bwa samse failed with return code {e.returncode}",
        }


@mcp.tool()
def bwa_sampe(
    in_db_fasta: Path,
    in1_sai: Path,
    in2_sai: Path,
    in1_fq: Path,
    in2_fq: Path,
    a: int = 500,
    o: int = 100000,
    n: int = 3,
    N: int = 10,
    P: bool = False,
    r: Optional[str] = None,
):
    """
    Generate alignments in the SAM format given paired-end reads using bwa sampe.
    -a INT: Maximum insert size for proper pair [500]
    -o INT: Maximum occurrences of a read for pairing [100000]
    -n INT: Max alignments in XA tag for properly paired reads [3]
    -N INT: Max alignments in XA tag for discordant pairs [10]
    -P: Load entire FM-index into memory
    -r STR: Specify the read group header line
    """
    for f in [in_db_fasta, in1_sai, in2_sai, in1_fq, in2_fq]:
        if not f.exists():
            raise FileNotFoundError(f"Input file {f} does not exist")
    if a < 0 or o < 0 or n < 0 or N < 0:
        raise ValueError("Parameters a, o, n, N must be non-negative")

    cmd = ["bwa", "sampe"]
    cmd += ["-a", str(a)]
    cmd += ["-o", str(o)]
    if P:
        cmd.append("-P")
    cmd += ["-n", str(n)]
    cmd += ["-N", str(N)]
    if r:
        r_fixed = r.replace("\\t", "\t")
        cmd += ["-r", r_fixed]
    cmd.append(str(in_db_fasta))
    cmd.append(str(in1_sai))
    cmd.append(str(in2_sai))
    cmd.append(str(in1_fq))
    cmd.append(str(in2_fq))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # bwa sampe outputs SAM to stdout
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
            "error": f"bwa sampe failed with return code {e.returncode}",
        }


@mcp.tool()
def bwa_bwasw(
    in_db_fasta: Path,
    in_fq: Path,
    mate_fq: Optional[Path] = None,
    a: int = 1,
    b: int = 3,
    q: int = 5,
    r: int = 2,
    t: int = 1,
    w: int = 33,
    T: int = 37,
    c: float = 5.5,
    z: int = 1,
    s: int = 3,
    N: int = 5,
):
    """
    Align query sequences using bwa bwasw (BWA-SW algorithm).
    Supports single-end and paired-end (Illumina short-insert) reads.
    """
    if not in_db_fasta.exists():
        raise FileNotFoundError(f"Input fasta file {in_db_fasta} does not exist")
    if not in_fq.exists():
        raise FileNotFoundError(f"Input fastq file {in_fq} does not exist")
    if mate_fq and not mate_fq.exists():
        raise FileNotFoundError(f"Mate fastq file {mate_fq} does not exist")
    if t < 1:
        raise ValueError("Number of threads 't' must be >= 1")
    if w < 1:
        raise ValueError("Band width 'w' must be >= 1")
    if T < 0:
        raise ValueError("Minimum score threshold 'T' must be >= 0")
    if c < 0:
        raise ValueError("Coefficient 'c' must be >= 0")
    if z < 1:
        raise ValueError("Z-best heuristics 'z' must be >= 1")
    if s < 1:
        raise ValueError("Maximum SA interval size 's' must be >= 1")
    if N < 0:
        raise ValueError("Minimum number of seeds 'N' must be >= 0")

    cmd = ["bwa", "bwasw"]
    cmd += ["-a", str(a)]
    cmd += ["-b", str(b)]
    cmd += ["-q", str(q)]
    cmd += ["-r", str(r)]
    cmd += ["-t", str(t)]
    cmd += ["-w", str(w)]
    cmd += ["-T", str(T)]
    cmd += ["-c", str(c)]
    cmd += ["-z", str(z)]
    cmd += ["-s", str(s)]
    cmd += ["-N", str(N)]
    cmd.append(str(in_db_fasta))
    cmd.append(str(in_fq))
    if mate_fq:
        cmd.append(str(mate_fq))

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        # bwa bwasw outputs SAM to stdout
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
            "error": f"bwa bwasw failed with return code {e.returncode}",
        }


if __name__ == '__main__':
    mcp.run()