#!/bin/bash

# Check if directory is provided
if [ -z "$1" ]; then
    echo "Please provide the directory name."
    exit 1
fi

# Get the directory name
DIR_NAME=$1

# Build the Docker image
docker build -t oblivious1/drap:$DIR_NAME $DIR_NAME

# Push the Docker image to the repository
docker push oblivious1/drap:$DIR_NAME

echo "Docker image built and pushed successfully."
