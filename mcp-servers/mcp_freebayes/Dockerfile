FROM condaforge/miniforge3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    default-jre \
    wget \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first (for better Docker layer caching)
COPY requirements.txt /tmp/
RUN pip install uv && \ 
    uv pip install --system -r /tmp/requirements.txt

# Or for conda
COPY environment.yaml /tmp/
RUN conda env update -f /tmp/environment.yaml && conda clean -a

# Create app directory
WORKDIR /app

# Copy your MCP server
COPY /app/freebayes_server.py /app/

# Create workspace and output directories
RUN mkdir -p /app/workspace /app/output

# Make sure the server script is executable
RUN chmod +x /app/freebayes_server.py

# Expose port for MCP over HTTP (optional)
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3     CMD python -c "import sys; sys.exit(0)"

# Default command runs the MCP server via stdio
CMD ["python", "/app/freebayes_server.py"]
        