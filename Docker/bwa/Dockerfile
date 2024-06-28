# Use an official base image with conda and essential tools
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Copy the environment file to the container
COPY environment.yml .

# Create the environment using the environment.yml file
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "sorbet_env", "/bin/bash", "-c"]

# Install additional tools
RUN conda install -c bioconda bwa samtools
RUN conda install -c conda-forge datamash
RUN conda install -c conda-forge pandas bc maven

# Activate the environment
RUN echo "conda activate sorbet_env" >> ~/.bashrc

# Copy the entrypoint script into the container

# Set the entrypoint
CMD ["bash"]
