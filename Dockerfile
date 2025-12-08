FROM python:3.10-slim

ARG VERSION=1.0.0
LABEL version=$VERSION

WORKDIR /app

COPY requirements.txt .

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    cmake \
    libopenblas-dev \
    libgomp1 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    git \
    && pip install --upgrade pip setuptools wheel \
    && pip install --no-cache-dir -r requirements.txt \
    && apt-get remove --purge -y git cmake build-essential gfortran \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /root/.cache

COPY . .

ENTRYPOINT ["python", "src/GenoJoin.py"]
