# Use a base image that includes Python
FROM python:3.8-slim

# Install dependencies required for BBDuk
RUN apt-get update && \
    apt-get install -y \
    default-jre \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Download and set up BBMap, which includes BBDuk
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.95.tar.gz/download -O BBMap.tar.gz && \
    tar -xvzf BBMap.tar.gz && \
    rm BBMap.tar.gz && \
    mv bbmap /opt/ && \
    ln -s /opt/bbmap/bbduk.sh /usr/local/bin/bbduk.sh

# Set the working directory
WORKDIR /app

# Copy the current directory contents into the container
COPY . .

# Define the entry point command to run the task
CMD ["bash"]
