FROM debian:bookworm-slim

ARG RAMMAP_VERSION=v1.1.2
ARG TARGETARCH

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl procps \
    && rm -rf /var/lib/apt/lists/*

RUN set -eux; \
    case "${TARGETARCH}" in \
        amd64) RAMMAP_ARCH=x86_64-unknown-linux-gnu ;; \
        arm64) RAMMAP_ARCH=aarch64-unknown-linux-gnu ;; \
        *) echo "Unsupported architecture: ${TARGETARCH}" >&2; exit 1 ;; \
    esac; \
    curl -fsSL -o /usr/local/bin/rammap \
        "https://github.com/jwanglab/rammap/releases/download/${RAMMAP_VERSION}/rammap_${RAMMAP_ARCH}_${RAMMAP_VERSION}"; \
    chmod +x /usr/local/bin/rammap

CMD ["/bin/bash"]
