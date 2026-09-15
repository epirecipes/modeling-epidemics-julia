# Base image for GitHub Codespaces on this course.
#
# Julia and Quarto are installed here rather than in a postCreate script so the
# versions land in the image layers, which Codespaces prebuilds cache.
FROM mcr.microsoft.com/devcontainers/base:2-ubuntu-24.04

# Quarto refuses to execute a document when the running Julia differs from the
# julia_version recorded in Manifest.toml, so the patch version is pinned here
# and must be changed together with the manifest.
ARG JULIA_VERSION=1.12.5
ARG QUARTO_VERSION=1.9.28

RUN set -eux; \
    case "$(dpkg --print-architecture)" in \
      amd64) jdir=x64;      jsuffix=x86_64  ;; \
      arm64) jdir=aarch64;  jsuffix=aarch64 ;; \
      *) echo "unsupported architecture" >&2; exit 1 ;; \
    esac; \
    curl -fsSL "https://julialang-s3.julialang.org/bin/linux/${jdir}/${JULIA_VERSION%.*}/julia-${JULIA_VERSION}-linux-${jsuffix}.tar.gz" \
      -o /tmp/julia.tar.gz; \
    mkdir -p /opt/julia; \
    tar -xzf /tmp/julia.tar.gz -C /opt/julia --strip-components=1; \
    rm /tmp/julia.tar.gz; \
    ln -s /opt/julia/bin/julia /usr/local/bin/julia

RUN set -eux; \
    curl -fsSL "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-$(dpkg --print-architecture).deb" \
      -o /tmp/quarto.deb; \
    apt-get update; \
    apt-get install -y --no-install-recommends /tmp/quarto.deb; \
    rm -f /tmp/quarto.deb; \
    rm -rf /var/lib/apt/lists/*

# Everything the two environments need in order to resolve and precompile. Only
# these paths are copied, so an edit to a chapter does not invalidate the layer
# that costs an hour to rebuild.
COPY --chown=vscode:vscode Project.toml Manifest.toml /opt/course/
COPY --chown=vscode:vscode src /opt/course/src
COPY --chown=vscode:vscode environments/ch10 /opt/course/environments/ch10
COPY --chown=vscode:vscode .devcontainer/build-depot.sh /opt/course/build-depot.sh

# Build the depot as the user who will be using it, so /home/vscode/.julia is
# owned correctly and is part of the published image.
USER vscode
RUN /opt/course/build-depot.sh
USER root
