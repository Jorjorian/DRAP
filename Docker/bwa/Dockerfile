# Use an official base image with conda and essential tools
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Create the environment using the environment.yml file
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "bwa_env", "/bin/bash", "-c"]

# Install additional tools within the environment
RUN conda install -n bwa_env -c bioconda bwa samtools datamash
RUN conda install -n bwa_env -c conda-forge pandas bc

# Install Picard
RUN conda install -n bwa_env -c bioconda picard

# Activate the environment by default
SHELL ["conda", "run", "-n", "bwa_env", "/bin/bash", "-c"]

# Set the entrypoint
CMD ["bash"]

