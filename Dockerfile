FROM --platform=linux/amd64 ubuntu:20.04

# Update and install basic tools
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    libgomp1 \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && bash miniconda.sh -b -p /opt/conda \
    && rm miniconda.sh

# Add conda to PATH
ENV PATH /opt/conda/bin:$PATH

# Initialize conda
RUN conda init bash

# Copy the environment.yaml file into the container
COPY Nanopore.yaml .

# Create the conda environment using the yaml file
RUN conda env create -n Nanopore -f Nanopore.yaml

# # Install pip packages (if needed)
# RUN pip install numpy pandas scikit-learn

# Set working directory
WORKDIR /app

# # Default command
# CMD ["/bin/bash"]

# Activate the conda environment
SHELL ["conda", "run", "-n", "Nanopore", "/bin/bash", "-c"]

# curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz
# tar -xJf canu-2.2.*.tar.xz