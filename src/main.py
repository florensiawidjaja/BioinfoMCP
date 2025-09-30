import os
import argparse
from bioinfomcp_converter import BioinfoMCP
import subprocess
import sys
from pathlib import Path
import json
import shutil

def generate_requirements_with_pipreqs(tool_name, server_path):
    """
    Use pipreqs to generate requirements.txt
    """
    # Install pipreqs if not available
    #subprocess.run([sys.executable, '-m', 'pip', 'install', 'pipreqs'], 
    #                 check=True, capture_output=True)
    # server_path = str(server_path) + "/app"
    
    # Run pipreqs on the server directory
    result = subprocess.run([
        'pipreqs', 
        str(server_path)
    ], capture_output=True, text=True, check=True)
    
    print(f"Generated requirements.txt using pipreqs")
    return Path(server_path) / "requirements.txt"

def generate_environment_yaml(tool_name, server_path):
    """
    generate environment.yaml
    """
    # Install pipreqs if not available
    #subprocess.run([sys.executable, '-m', 'pip', 'install', 'pipreqs'], 
    #                 check=True, capture_output=True)
    # server_path = str(server_path) + "/app"
    
    # Run pipreqs on the server directory
    result = f"""
name: mcp-tool
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - {tool_name}
  - python=3.10
    """
    with open(server_path / "environment.yaml", "w") as f:
        f.write(result)
    print(f"Generated environment.yaml using pipreqs")
    return Path(server_path) / "environment.yaml"

def generate_environment_yml(tool_name, server_path):
    """
    Generate requirements.yml
    """
    result = f"""
        name: {tool_name}_env
        channels:
        - bioconda
        - conda-forge
        dependencies:
        - python=3.10
        - {tool_name}
        - pip
    """
    return result

def convert_mcptool(tool_name, manual, run_help_command, server_path):    
    # parser = argparse.ArgumentParser()
    # bioinfo_tools_ls = ['fastqc','trimmomatic']
    bioinfo_tools_ls = [tool_name]

    converter = BioinfoMCP()
    for tool in bioinfo_tools_ls:
        conv_result = converter.autogenerate_mcp_tool(tool, manual, run_help_command)
        while not conv_result[0]:
            if conv_result[2] is None:
                # There is no python code in the result
                conv_result = converter.autogenerate_mcp_tool(tool, manual, run_help_command)
                print("generate1")
            else:
                print("generate2")
                conv_result = converter.refine_after_feedback(tool, code=conv_result[2], error_message=conv_result[1])
        output_file = open(f'{server_path}/app/{tool}_server.py', 'w')
        output_file.write(conv_result[2])

def dockerfile_content(tool_name):
    return f"""
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
        default-jre \
        wget \
        curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Install {tool_name} via conda (e.g., from bioconda)
RUN conda install -c bioconda {tool_name} -y && \
    conda clean -a

# Install Python dependencies
RUN pip install uv
RUN uv pip install --system fastmcp

# Create app directory
WORKDIR /app

# Copy your MCP server
COPY {tool_name}_server.py /app/

# Create workspace and output directories
RUN mkdir -p /app/workspace /app/output

# Make sure the server script is executable
RUN chmod +x /app/{tool_name}_server.py

# Expose port for MCP over HTTP (optional)
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import sys; sys.exit(0)"

# Default command runs the MCP server via stdio
CMD ["python", "/app/{tool_name}_server.py"]
        """

def dockercompose_content(tool_name):
    """Generate docker-compose.yml for easy deployment"""
    
    compose = f"""version: '3.8'

services:
  mcp-{tool_name}:
    build: .
    image: mcp-{tool_name}:latest
    container_name: mcp-{tool_name}
    ports:
      - "8000:8000"
    environment:
      - MCP_SERVER_NAME={tool_name}
    volumes:
      - ./workspace:/app/workspace
      - ./output:/app/output
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "python", "-c", "import sys; sys.exit(0)"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 5s
    """

    return compose

def build_docker_image(tool_name, server_path, output_path, is_pipeline):

    # Create build directory
    try:
        # server_path.mkdir(parents=True, exist_ok=True)
        # print(server_path)
        build_dir = Path(output_path) / f"mcp_{tool_name}"

        # shutil.copy2(server_path / f"{tool_name}_server.py", build_dir / f"app", dirs_exist_ok=True)
        #shutil.copy2(server_path, build_dir / "app", dirs_exist_ok=True)
        df_content = dockerfile_content(tool_name)
        with open(server_path / "Dockerfile", "w") as f:
            f.write(df_content) 
                
         # make the docker-compose.yml
        if not is_pipeline:
            dc_content = dockercompose_content(tool_name)
            with open(server_path/"docker-compose.yml", "w") as f:
                f.write(dc_content)
        
        # subprocess.run(["docker", "build", "-t", f"{args.name}-docker", server_path])

        return 1
    except:
        return 0


def claude_addition(tool_name):
    config = {
        "mcpServers": {
            tool_name: {
                "command": "docker",
                "args": [
                    "run",
                    "--rm",
                    "-i",
                    "-v", 
                    f"/local/path/{tool_name}/data:/app/workspace",
                    f"{tool_name}-mcp-server:latest"
                ]
            }
        }
    }
    return json.dumps(config, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Accept the Bioinformatic tool name and the help document")
    parser.add_argument('--name', type=str)
    parser.add_argument('--manual', type=str, help="The file path to the help document")
    parser.add_argument('--run_help_command', type=bool, default=False)
    parser.add_argument('--output_location', type=str)
    parser.add_argument('--is_pipeline', action='store_true', default=False)
    args = parser.parse_args()
    '''
    if help is False, then manual is the attribute of the tool to access the help docs (ex. fastqc --help; manual ='--help')
    --> fastqc must be installed already
    if help is True, then manual is the file path to the help document
    '''
    server_path = Path(args.output_location) / f"mcp_{args.name}"
    os.makedirs(server_path, exist_ok=True)
    print("Server path", server_path)
    app_path = Path(server_path) / "app"
    os.makedirs(app_path, exist_ok=True)

    # run the converter
    convert_mcptool(args.name, args.manual, args.run_help_command, server_path)

    # Write the requirements.txt
    generate_requirements_with_pipreqs(args.name, server_path)
    generate_environment_yaml(args.name, server_path)

    # build docker image    
    status = build_docker_image(args.name, server_path, args.output_location, args.is_pipeline)
    
    if status:
        add = claude_addition(args.name)

        print(f"Add the following onto your Claude Configuration json file to run the MCP server:\n{'=='*10}\n{add}\n{'=='*10}")
    else:
        print("Failed to build Docker Image")

    




