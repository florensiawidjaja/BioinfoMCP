from pathlib import Path
from typing import Optional, List

f: List[Path] = []
f.append(Path("/Users/florensiawidjaja/Documents/BioInfoMCP/sample1.fastq"))
print(f)
for fi in f:
    # print(fi)
    print(fi.exists())