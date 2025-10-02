# A Unified Platform Enabling MCP Interfaces in Agentic Bioinformatics 🧬🤖


Florensia Widjaja, Zhangtianyi Chen, Juexiao Zhou

The Chinese University of Hong Kong, Shenzhen (CUHK Shenzhen)

<!-- <a href='media/advs.202407094.pdf'><img src='https://img.shields.io/badge/Paper-PDF-red'></a> -->

![Main Figure of the BioinfoMCP Platform](https://github.com/florensiawidjaja/BioinfoMCP/blob/24a6755dd976f4e31cff1fddc6b24d28f03cc873/media/main-fig.png "the BioinfoMCP Platform")


## 🚀 About
BioinfoMCP bridges the gap between specialized bioinformatics command-line tools and modern AI-driven workflows. The platform consists of two main components:

**BioinfoMCP Converter** - Automatically generates robust, production-ready MCP servers from tool documentation using large language models. Simply provide a tool's manual or help documentation, and the converter produces a fully functional MCP server with comprehensive parameter support.

**BioinfoMCP Benchmark** - A curated validation suite that systematically tests converted MCP servers across diverse computational biology tasks, ensuring reliability and versatility across different AI agent platforms.

## 🧐 What's New
- **[2025/10]** We finally share our ground-breaking BioinfoMCP Platform to the scientific community (version `v0.0.1`)! 🥳👏

## TODO list

We're working hard to achieve more features, welcome to PRs!

- [x] Provide a docker version, simplify the installation process
- [ ] Pack into a conda package, simplify the installation process
- [ ] Automatic tool discovery and conversion
- [ ] BioinfoMCP Benchmark Addition: Ollama
- [ ] Online Platform to share MCP servers
- [ ] Adaptive and More reliable Installation 

**We warmly welcome any inputs or contributions to collaboratively improve BioinfoMCP**


## Installation

### Command line
```shell
# 1. Clone this project
git clone https://github.com/florensiawidjaja/BioinfoMCP.git

# 2. Create Conda environment for an isolated production environment
conda create -n bioinfomcp-env
conda activate bioinfomcp-env

# 3. Install the fastmcp package
uv pip install fastmcp

```

## 🗃️ Repository Structure


## 🧰🔬 BioinfoMCP Converter - Convert your own Bioinformatics Tool
#### Prerequisite
1. An OpenAI API Key from an account that has a sufficient credit available.
*Note: You should also be able to use Anthropic Models. However, we have not done extensive experiments on it so we still suggest you to use an OpenAI model.*

### Usage Manual
1. Modify the `.env` file with your API Key and Model name (We recommend model with a capability of at least GPT-4.1-mini be used)

Run this command to start a simple example with GPT-4.1-mini as backend (**recommended**). 

**Please ensure that the LLM backbone supports function calling**

Run this script if you have checked the internal manual and it contains comprehensive documentation.
`python -m main --name <tool-name> --manual "\-\-help" --run_help_command True --output_location </path/to/output/folder/>`

Run this script if you want BioinfoMCP Converter to extract from a PDF manual.
`python -m main --name trimmomatic --manual "</path/to/tool/manual.pdf>" --output_location </path/to/output/folder/>`

After having the converted MCP server ready, you can now use your MCP server as follows:

1. **If connecting directly to your Python Environment**
* Activate your Conda Environment and install the tool that you want to use

```shell
conda activate bioinfomcp-env

# Bioinformatics tools are usually integrated inside the bioconda channel, but search for conda install <tool-name> to ensure that it is correct
conda install bioconda::tool-name

# Check whether it is appropriately installed already
conda list tool-name
```

* Connect to your AI Agent with the following JSON Configuration
```json
"tool-name": {
      "command": "bash",
      "args": [
        "-c",
        "source /Users/<username>/miniforge3/etc/profile.d/conda.sh && conda activate bioinfomcp-env && python /path/to/mcp_<tool-name>/app/<tool-name>_server.py"
      ]
    }
```
2. If connecting through a Docker Container (*Make sure that you have Docker CLI and Engine installed and ready to use*)

* Turn on your Docker Engine
* Build your Docker Container
```shell
cd /path/to/mcp-tool-name/folder/
docker compose up --build
```
* Check that the `mcp-tool-name` container is up and running

* Connect to your AI Agent with the following JSON Configuration
```json
"tool-name": {
      "command": "docker",
      "args": [
        "run",
        "--rm",
        "-i",
        "-v",
        "/path/to/mcp_<tool-name>/data:/app/workspace",
        "mcp-<tool-name>:latest"
      ]
    }
```


## 📝💊 BioinfoMCP Benchmark - Test Your MCP Server Reliability
**Note: **

```json
"filesystem": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-filesystem",
        "/path/to/folder1/",
        "/path/to/folder2/"
      ]
    }
```

### 1. Individual Tool Testing

#### 1.1 FastQC
```
I want to run FastQC on /path/to/SRR8405197.fastq. Please state what commands you run and what are the outputs or results from that command.
```

### 2. Pipeline/Complex Task Testing
> **Note**: Make sure that you have all sufficient tools connected to your AI agent to run this pipeline. We give you reference of what tools will be needed, but feel free to replace or add more tools on your own flexibly.

#### 2.1 Pure CLI tools: ATAC-seq
> Tools needed: FastQC, Trim-galore, Bowtie2, samtools, MACS3, MultiQC 
```
I want to run pipeline ATAC-seq on /path/to/genomic/files/ to do identification of open chromatin region. Please run the command appropriately and please explain and give appropriate next steps.
```

### 2.2 Combination with R Package tools: ChIP-seq
For certain pipelines like ChIP-seq, we would require tools that are already commonly used in R language.
**We used IMNMV/ClaudeR to extend a direct link between RStudio and AI Agents**
reference: https://github.com/IMNMV/ClaudeR

* Connect CluadeR to your AI Agent using the instruction shown in the reference above.
* Restart your AI Agent and then run the following:
```
I want to run pipeline ChIP-seq on /path/to/genomic/files/ to do motif discovery for binding
sites. Use the existing MCP tools as much as possible, and only execute R if none of them are available for the task given. Make an Executive Summary report at the end.
```


##### Download Data
In our paper, we mainly use SRR8405197.fastq, which can be obtained as follows:
```shell
# For single-ends
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR840/007/SRR8405197/SRR8405197.fastq.gz

# For paired-ends
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR840/007/SRR8405197/SRR8405197_2.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR840/007/SRR8405197/SRR8405197_1.fastq.gz

# For reference genome (if needed)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# For gene annotation (if needed)
# refer to https://www.gencodegenes.org/human/release_36.html
```


## Citation

If you find this project useful in your research, please consider citing:

...


## License

This project is released under the MIT license.
