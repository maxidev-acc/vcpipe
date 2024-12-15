# Use an official Ubuntu as a base image
FROM ubuntu:20.04

# Set environment variables to avoid interaction during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update the system and install necessary dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    curl \
    wget \
    git \
    openjdk-11-jdk \
    samtools \
    bwa \
    && apt-get clean

# Install GATK
RUN mkdir -p /opt/gatk && \
    wget https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip -O /opt/gatk/gatk.zip && \
    unzip /opt/gatk/gatk.zip -d /opt/gatk && \
    rm /opt/gatk/gatk.zip

# Install Picard (download the latest version from the GitHub releases)
RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.9/picard-2.26.9.jar -O /opt/picard.jar

# Install snpEff
RUN wget https://sourceforge.net/projects/snpeff/files/snpEff_v5_0/snpEff_v5_0_core.zip -O /opt/snpEff.zip && \
    unzip /opt/snpEff.zip -d /opt/snpEff && \
    rm /opt/snpEff.zip

# Install Python dependencies
RUN pip3 install --no-cache-dir \
    json

# Set working directory
WORKDIR /data

# Set PATH for Java, bwa, samtools, picard, gatk, and snpEff
ENV PATH=$PATH:/opt/gatk/gatk-4.2.5.0:/opt/snpEff/snpEff:/opt/java/openjdk-11/bin

# Add the bash script and any required files to the Docker image
COPY pipeline.sh /data/pipeline.sh

# Ensure that the script is executable
RUN chmod +x /data/pipeline.sh

# Entry point
ENTRYPOINT ["/data/your_script.sh"]
