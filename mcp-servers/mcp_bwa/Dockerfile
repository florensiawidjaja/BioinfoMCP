
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
CMD ["python", "bwa_server.py"]
        