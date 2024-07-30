# Use an official Ubuntu base image
FROM ubuntu:latest

# Set the working directory in the container
WORKDIR /app

# Update the package list and install jq
RUN apt-get update && apt-get install -y jq unzip

# Set the entrypoint to bash
ENTRYPOINT ["bash"]
