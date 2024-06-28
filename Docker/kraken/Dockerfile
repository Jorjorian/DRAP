# Use the base image
FROM staphb/kraken2:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install Python and other dependencies
RUN apt-get update && \
    apt-get install -y python3 python3-pip && \
    apt-get clean

# Optionally, install conda if you need to activate a conda environment
# You can skip this if you are not using conda in your command script

# Install any Python packages if needed
RUN pip3 install --no-cache-dir numpy pandas

# Set the working directory (optional)
WORKDIR /data

# Default command (optional)
CMD ["bash"]
