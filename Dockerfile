FROM ubuntu:22.04 as mash_build

RUN apt update && \
    apt install -y curl && \
     rm -rf /var/lib/apt/lists/*

RUN curl -L https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar | tar xv \
    && mv mash-Linux64-v2.3/mash /usr/bin/ \
    && rm -rf mash-Linux64-v2.3

FROM inspeciator-data:v0.0.1 as libraries

FROM python:3.10

COPY --from=mash_build /usr/bin/mash /usr/bin/mash

RUN pip install pandas pyarrow
RUN mkdir -p /libraries/

COPY --from=libraries /data/ /libraries/

COPY entrypoint.sh /

COPY speciator.py /

COPY bactinspector /bactinspector

RUN chmod a+x /entrypoint.sh

ENTRYPOINT ["/bin/bash", "/entrypoint.sh"]