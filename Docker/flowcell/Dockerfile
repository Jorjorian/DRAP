# Use Miniconda3 base image
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

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Default command to run when starting the container
CMD ["bash"]

