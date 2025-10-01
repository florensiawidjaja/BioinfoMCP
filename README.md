# A Unified Platform Enabling MCP Interfaces in Agentic Bioinformatics 🧬🤖


Florensia Widjaja, Zhangtianyi Chen, Juexiao Zhou

The Chinese University of Hong Kong, Shenzhen (CUHK Shenzhen)

<!-- <a href='media/advs.202407094.pdf'><img src='https://img.shields.io/badge/Paper-PDF-red'></a> -->

![Main Figure of the BioinfoMCP Platform](https://github.com/florensiawidjaja/BioinfoMCP/blob/d0dd78d8d074b4450f0add65b69426c1302672fa/media/main-fig.svg)


## 🚀 About
BioinfoMCP bridges the gap between specialized bioinformatics command-line tools and modern AI-driven workflows. The platform consists of two main components:

**BioinfoMCP Converter** - Automatically generates robust, production-ready MCP servers from tool documentation using large language models. Simply provide a tool's manual or help documentation, and the converter produces a fully functional MCP server with comprehensive parameter support.

**BioinfoMCP Benchmark** - A curated validation suite that systematically tests converted MCP servers across diverse computational biology tasks, ensuring reliability and versatility across different AI agent platforms.

## 🧐 What's New
- **[2025/09]** We finally share our ground-breaking BioinfoMCP Platform to the scientific community (version `v0.0.1`)! 🥳👏

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
git clone https://github.com/florensiawidjaja/BioinfoMCP.git

# Create Conda environment for an isolated production environment
conda create -n bioinfomcp-env
conda activate bioinfomcp-env

uv pip install fastmcp


```

### Docker

Please refer to https://docs.docker.com/engine/install to install Docker first.

```shell

```


Try the previous codes again.

### Conda
```shell
Coming soon...
```

## Get Started

### Understand files


### BioinfoMCP Converter

Run this command to start a simple example with GPT-4.1-mini as backend (**recommended**).

Run this script if you have checked the internal manual and it contains comprehensive documentation.
`python -m main --name <tool-name> --manual "\-\-help" --run_help_command True --output_location <output-folder-location>`

Run this script if you want BioinfoMCP Converter to extract from a PDF manual.
`python -m main --name trimmomatic --manual "</path/to/tool/manual.pdf>" --output_location <output-folder-location>`

**Please ensure that the LLM backbone supports function calling**

### BioinfoMCP Benchmark

### Example 1: Bulk RNA-Seq

#### Case 1.1: Find differentially expressed genes

**Reference**: https://pzweuj.github.io/worstpractice/site/C02_RNA-seq/01.prepare_data/

Design of `config.yaml`
```yaml
data_list: [ './examples/case1.1/data/SRR1374921.fastq.gz: single-end mouse rna-seq reads, replicate 1 in LoGlu group',
            './examples/case1.1/data/SRR1374922.fastq.gz: single-end mouse rna-seq reads, replicate 2 in LoGlu group',
            './examples/case1.1/data/SRR1374923.fastq.gz: single-end mouse rna-seq reads, replicate 1 in HiGlu group',
            './examples/case1.1/data/SRR1374924.fastq.gz: single-end mouse rna-seq reads, replicate 2 in HiGlu group',
            './examples/case1.1/data/TruSeq3-SE.fa: trimming adapter',
            './examples/case1.1/data/mm39.fa: mouse mm39 genome fasta',
            './examples/case1.1/data/mm39.ncbiRefSeq.gtf: mouse mm39 genome annotation' ]
output_dir: './examples/case1.1/output'
goal_description: 'find the differentially expressed genes'
```

##### Download Data
```shell
wget -P data/ http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
wget -P data/ http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip data/mm39.ncbiRefSeq.gtf.gz
gunzip data/mm39.fa.gz
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/001/SRR1374921/SRR1374921.fastq.gz
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/002/SRR1374922/SRR1374922.fastq.gz
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/003/SRR1374923/SRR1374923.fastq.gz
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR137/004/SRR1374924/SRR1374924.fastq.gz
```

##### Analyze with AutoBA

```shell
python app.py --config ./examples/case1.1/config.yaml --openai YOUR_OPENAI_API --model gpt-4
python app.py --config ./examples/case1.1/config.yaml --model codellama-7bi
python app.py --config ./examples/case1.1/config.yaml --model codellama-13bi
python app.py --config ./examples/case1.1/config.yaml --model codellama-34bi
```

## Citation

If you find this project useful in your research, please consider citing:

...


## License

This project is released under the MIT license.
# BioinfoMCP