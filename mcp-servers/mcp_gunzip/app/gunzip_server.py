from fastmcp import FastMCP
from pathlib import Path
from typing import List, Optional
import subprocess

mcp = FastMCP()


@mcp.tool()
def gunzip(
    files: List[Path],
    ascii: bool = False,
    stdout: bool = False,
    decompress: bool = True,
    force: bool = False,
    list_files: bool = False,
    license: bool = False,
    no_name: bool = False,
    name: bool = False,
    quiet: bool = False,
    recursive: bool = False,
    suffix: Optional[str] = None,
    test: bool = False,
    verbose: bool = False,
    version: bool = False,
) -> dict:
    """
    Decompress files compressed with gzip, zip, compress, or pack formats.
    gunzip replaces each compressed file with an uncompressed file by default.
    Supports multiple files and recursive directory traversal.

    Parameters:
    - files: List of input files or directories to decompress.
    - ascii: Convert end-of-lines using local conventions (only on some non-Unix systems).
    - stdout: Write output to standard output; keep original files unchanged.
    - decompress: Decompress mode (default True).
    - force: Force decompression even if file has multiple links or output file exists.
    - list_files: List compressed file info instead of decompressing.
    - license: Display gzip license and quit.
    - no_name: Do not restore original file name or timestamp when decompressing.
    - name: Restore original file name and timestamp if present.
    - quiet: Suppress all warnings.
    - recursive: Recursively decompress files in directories.
    - suffix: Use specified suffix instead of default (.gz) when decompressing.
    - test: Test compressed file integrity.
    - verbose: Display detailed info during decompression.
    - version: Display version info and quit.

    Returns:
    A dictionary with keys: command_executed, stdout, stderr, output_files.
    """
    # Validate input files/directories
    if not files:
        raise ValueError("At least one input file or directory must be specified.")
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"Input file or directory does not exist: {f}")

    # Build command line
    cmd = ["gunzip"]

    # Flags mapping
    if ascii:
        cmd.append("-a")
    if stdout:
        cmd.append("-c")
    if decompress:
        cmd.append("-d")
    if force:
        cmd.append("-f")
    if list_files:
        cmd.append("-l")
    if license:
        cmd.append("-L")
    if no_name:
        cmd.append("-n")
    if name:
        cmd.append("-N")
    if quiet:
        cmd.append("-q")
    if recursive:
        cmd.append("-r")
    if test:
        cmd.append("-t")
    if verbose:
        cmd.append("-v")
    if version:
        cmd.append("-V")
        # If version or license requested, no files needed
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
                "error": str(e),
            }

    # Handle suffix option
    if suffix:
        if not suffix.startswith("."):
            raise ValueError("Suffix must start with a dot '.'")
        cmd.extend(["-S", suffix])

    # Add input files/directories as strings
    for f in files:
        cmd.append(str(f))

    # Execute command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": str(e),
        }

    # Determine output files if not stdout and not list or test or license or version
    output_files = []
    if not stdout and not list_files and not test and not license and not version:
        # If recursive and directories given, collect decompressed files
        # Otherwise, for each input file, output file is input file without compression suffix
        suffixes = [".gz", "-gz", ".z", "-z", "_z", ".Z"]
        if suffix:
            suffixes.insert(0, suffix)
        def decompressed_name(p: Path) -> Path:
            name = p.name
            for suf in suffixes:
                if name.endswith(suf):
                    return p.with_name(name[: -len(suf)])
            # If no recognized suffix, return original path (gunzip may fail)
            return p

        if recursive:
            # Collect all decompressed files under directories
            for f in files:
                if f.is_dir():
                    for p in f.rglob("*"):
                        if p.is_file():
                            outp = decompressed_name(p)
                            if outp != p:
                                output_files.append(str(outp.resolve()))
                else:
                    outp = decompressed_name(f)
                    if outp != f:
                        output_files.append(str(outp.resolve()))
        else:
            for f in files:
                if f.is_file():
                    outp = decompressed_name(f)
                    if outp != f:
                        output_files.append(str(outp.resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


@mcp.tool()
def gzip(
    files: List[Path],
    ascii: bool = False,
    stdout: bool = False,
    compress: bool = True,
    force: bool = False,
    list_files: bool = False,
    license: bool = False,
    no_name: bool = False,
    name: bool = True,
    quiet: bool = False,
    recursive: bool = False,
    suffix: Optional[str] = None,
    test: bool = False,
    verbose: bool = False,
    version: bool = False,
    fast: bool = False,
    best: bool = False,
    compression_level: int = 6,
) -> dict:
    """
    Compress files using gzip compression with Lempel-Ziv coding (LZ77).
    Supports multiple files and recursive directory traversal.

    Parameters:
    - files: List of input files or directories to compress.
    - ascii: Convert end-of-lines using local conventions (only on some non-Unix systems).
    - stdout: Write output to standard output; keep original files unchanged.
    - compress: Compress mode (default True).
    - force: Force compression even if file has multiple links or output file exists.
    - list_files: List compressed file info instead of compressing.
    - license: Display gzip license and quit.
    - no_name: Do not save original file name and timestamp when compressing.
    - name: Always save original file name and timestamp (default True).
    - quiet: Suppress all warnings.
    - recursive: Recursively compress files in directories.
    - suffix: Use specified suffix instead of default (.gz) when compressing.
    - test: Test compressed file integrity.
    - verbose: Display detailed info during compression.
    - version: Display version info and quit.
    - fast: Use fastest compression method (-1).
    - best: Use best compression method (-9).
    - compression_level: Compression level from 1 (fast) to 9 (best), default 6.

    Returns:
    A dictionary with keys: command_executed, stdout, stderr, output_files.
    """
    # Validate input files/directories
    if not files:
        raise ValueError("At least one input file or directory must be specified.")
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"Input file or directory does not exist: {f}")

    # Validate compression_level
    if compression_level < 1 or compression_level > 9:
        raise ValueError("compression_level must be between 1 and 9")

    # Build command line
    cmd = ["gzip"]

    # Flags mapping
    if ascii:
        cmd.append("-a")
    if stdout:
        cmd.append("-c")
    if compress:
        # compress is default, no flag needed
        pass
    if force:
        cmd.append("-f")
    if list_files:
        cmd.append("-l")
    if license:
        cmd.append("-L")
    if no_name:
        cmd.append("-n")
    if name:
        cmd.append("-N")
    if quiet:
        cmd.append("-q")
    if recursive:
        cmd.append("-r")
    if test:
        cmd.append("-t")
    if verbose:
        cmd.append("-v")
    if version:
        cmd.append("-V")
        # If version or license requested, no files needed
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
                "error": str(e),
            }

    # Handle suffix option
    if suffix:
        if not suffix.startswith("."):
            raise ValueError("Suffix must start with a dot '.'")
        cmd.extend(["-S", suffix])

    # Handle compression level flags
    if fast and best:
        raise ValueError("Cannot specify both fast and best compression options.")
    if fast:
        cmd.append("-1")
    elif best:
        cmd.append("-9")
    else:
        # Use compression_level if not default 6
        if compression_level != 6:
            cmd.append(f"-{compression_level}")

    # Add input files/directories as strings
    for f in files:
        cmd.append(str(f))

    # Execute command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": str(e),
        }

    # Determine output files if not stdout and not list or test or license or version
    output_files = []
    if not stdout and not list_files and not test and not license and not version:
        # If recursive and directories given, collect compressed files
        # Otherwise, for each input file, output file is input file with suffix appended
        suf = suffix if suffix else ".gz"

        def compressed_name(p: Path) -> Path:
            return p.with_name(p.name + suf)

        if recursive:
            for f in files:
                if f.is_dir():
                    for p in f.rglob("*"):
                        if p.is_file():
                            output_files.append(str(compressed_name(p).resolve()))
                else:
                    output_files.append(str(compressed_name(f).resolve()))
        else:
            for f in files:
                if f.is_file():
                    output_files.append(str(compressed_name(f).resolve()))

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": output_files,
    }


@mcp.tool()
def zcat(
    files: List[Path],
    force: bool = False,
    quiet: bool = False,
    verbose: bool = False,
    version: bool = False,
) -> dict:
    """
    Uncompress files to standard output without removing original files.
    Equivalent to 'gunzip -c'.

    Parameters:
    - files: List of input compressed files.
    - force: Force decompression even if file has multiple links or output file exists.
    - quiet: Suppress all warnings.
    - verbose: Display detailed info during decompression.
    - version: Display version info and quit.

    Returns:
    A dictionary with keys: command_executed, stdout, stderr, output_files (empty).
    """
    if not files:
        raise ValueError("At least one input file must be specified.")
    for f in files:
        if not f.exists():
            raise FileNotFoundError(f"Input file does not exist: {f}")

    cmd = ["zcat"]

    if force:
        cmd.append("-f")
    if quiet:
        cmd.append("-q")
    if verbose:
        cmd.append("-v")
    if version:
        cmd.append("-V")
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
                "error": str(e),
            }

    for f in files:
        cmd.append(str(f))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        return {
            "command_executed": " ".join(cmd),
            "stdout": e.stdout,
            "stderr": e.stderr,
            "output_files": [],
            "error": str(e),
        }

    return {
        "command_executed": " ".join(cmd),
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_files": [],
    }


if __name__ == '__main__':
    mcp.run()