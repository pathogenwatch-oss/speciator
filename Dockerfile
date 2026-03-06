FROM ghcr.io/astral-sh/uv:python3.13-bookworm-slim AS build

COPY src /speciator/src
COPY uv.lock pyproject.toml config.toml LICENSE.md README.md /speciator/

WORKDIR /speciator

RUN uv build --wheel && \
    rm -rf src

FROM ghcr.io/astral-sh/uv:python3.13-bookworm-slim AS code

COPY --from=build /speciator/config.toml /speciator/LICENSE.md /speciator/README.md /speciator/
COPY --from=build /speciator/dist/speciator-*.whl /speciator/dist/

WORKDIR /speciator

RUN apt update && \
    apt install -y curl --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*

RUN curl -L https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar | tar xv \
    && mv mash-Linux64-v2.3/mash /usr/bin/ \
    && rm -rf mash-Linux64-v2.3

RUN uv pip install --system --no-cache-dir /speciator/dist/*.whl && \
    rm -rf /speciator/dist

ENTRYPOINT ["uv", "run", "speciator-build-lib"]

FROM code AS speciator
# This stage assumes the the database has been built and is available in `library`.
# Build a code image to create the library using Docker.

COPY library /speciator/library

ENTRYPOINT ["uv", "run", "speciator"]

FROM speciator AS pathogenwatch

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]
