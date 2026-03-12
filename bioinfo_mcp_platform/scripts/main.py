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

    """
    将生物信息学工具转换为MCP工具服务器的函数

    参数:
        tool_name (str): 要转换的工具名称
        manual (str): 工具的手册页内容
        run_help_command (function): 执行帮助命令的函数
        server_path (str): 服务器文件保存路径

    返回:
        None: 该函数不返回值，但会在指定路径生成工具服务器文件
    """
    # parser = argparse.ArgumentParser()  # 这行被注释掉了，可能是之前用于解析命令行参数的代码
    # bioinfo_tools_ls = ['fastqc','trimmomatic']  # 这行被注释掉了，可能是之前定义的工具列表
    bioinfo_tools_ls = [tool_name]  # 只处理传入的单个工具
    # converter = BioinfoMCP(api_model_name="deepseek-chat", api_key="sk-530b6639c0384edbbb35ee9801ade10e")
    # converter = BioinfoMCP(api_model_name='gemini-2.5-flash', api_key="AIzaSyC-9og_9OsxvKZ0rBXeMGboXBrMOpG5-do")  # 创建生物信息学工具转换器实例
    converter = BioinfoMCP(api_model_name='gemini-2.5-flash', api_key="AIzaSyDClRNJkcDgHv2wA90v6TODPvBlu8umIWU")  # 创建生物信息学工具转换器实例
    for tool in bioinfo_tools_ls:  # 遍历工具列表（这里只有一个工具）
        # 尝试自动生成MCP工具
        conv_result = converter.autogenerate_mcp_tool(tool, manual, run_help_command)
        # 如果生成失败，则循环尝试重新生成或改进
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
# Use an official Python runtime as the base image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /app

# Copy requirements first to leverage Docker layer caching
COPY requirements.txt .

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the application
COPY . .

# Expose the port MCP runs on (e.g., 8000)
EXPOSE 8000

# Command to run the MCP server (adjust as needed)
CMD ["python", "{tool_name}_server.py"]
        """

def dockerfile2_content(tool_name):
    return f"""
FROM continuumio/miniconda3

WORKDIR /app

# Install uv first
RUN pip install uv

# Copy conda environment file
COPY environment.yml .

# Install bioinformatics tools via conda
RUN conda env create -f environment.yml

# Activate the conda environment and install fastmcp via uv
RUN conda run -n {tool_name}_env uv pip install fastmcp

# Copy your requirements for any additional Python packages
# COPY requirements.txt .

# Install additional Python packages if needed
# RUN conda run -n {tool_name}_env pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

EXPOSE 8000

# Run with conda environment activated
CMD ["conda", "run", "-n", "{tool_name}_env", "python", "{tool_name}_server.py"]
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
      - ./data:/data
      - ./config:/config
    restart: unless-stopped

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
                f"mcp/{tool_name}",
            ]
            }
        }
    }
    return json.dumps(config, indent=2)

def workflow(name, manual, run_help_command, output_location):
    server_path = Path(output_location) / f"mcp_{name}"
    os.makedirs(server_path, exist_ok=True)
    print("Server path", server_path)
    app_path = Path(server_path) / "app"
    os.makedirs(app_path, exist_ok=True)

    # run the converter
    convert_mcptool(name, manual, run_help_command, server_path)

    # Write the requirements.txt
    # generate_requirements_with_pipreqs(args.name, server_path)

    # build docker image    
    status = build_docker_image(name, server_path, output_location, is_pipeline)
    
    if status:
        add = claude_addition(name)

        print(f"Add the following onto your Claude Configuration json file to run the MCP server:\n{'=='*10}\n{add}\n{'=='*10}")
    else:
        print("Failed to build Docker Image")


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
    # generate_requirements_with_pipreqs(args.name, server_path)

    # build docker image    
    status = build_docker_image(args.name, server_path, args.output_location, args.is_pipeline)
    
    if status:
        add = claude_addition(args.name)

        print(f"Add the following onto your Claude Configuration json file to run the MCP server:\n{'=='*10}\n{add}\n{'=='*10}")
    else:
        print("Failed to build Docker Image")

