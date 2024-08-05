# Use an official Ubuntu base image
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the environment file into the container
COPY environment.yml ./

# Initialize conda and create the environment
RUN conda init bash && \
    conda env create -f environment.yml && \
    echo "conda activate myenv" >> ~/.bashrc

# Make sure the environment is activated
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]
RUN apt-get update && apt-get install -y jq unzip

# Set the entrypoint to bash
ENTRYPOINT ["bash"]
