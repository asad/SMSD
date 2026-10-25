FROM python:3.14.9-slim-bookworm@sha256:7e0522950230cad1ec100ada6e3898c22c6bc82851289b4bdf0bf77c1c829cd5 AS artifacts
ARG TARGETARCH
ARG SMSD_REVISION
WORKDIR /artifacts
COPY source.tar.gz SOURCE_MANIFEST.json ARTIFACTS.json ./
COPY LICENSE NOTICE ./
COPY wheels/ ./wheels/
COPY scripts/docker_install.py scripts/docker_smoke.py ./scripts/
RUN python scripts/docker_install.py --verify-only --root /artifacts --architecture "$TARGETARCH" --revision "$SMSD_REVISION"

FROM maven:3.9.14-eclipse-temurin-25@sha256:304ea1e9675d07612c36d5d70d8bd1dbd8aac64c307253bbde4e488857f90faa AS build
WORKDIR /app
COPY --from=artifacts /artifacts/source.tar.gz /tmp/source.tar.gz
RUN tar -xzf /tmp/source.tar.gz --strip-components=1 -C /app \
    && mvn -f java/pom.xml -U -B -DskipTests=true clean package

FROM eclipse-temurin:25-jre@sha256:dbb43d81099eec91e0800226ed2296d03146ee7857165b95112372ad0d67e87f
ARG TARGETARCH
ARG SMSD_REVISION
LABEL org.opencontainers.image.title="SMSD Pro" \
      org.opencontainers.image.description="Command-line and Python substructure and MCS search, fingerprints and similarity screening" \
      org.opencontainers.image.version="8.0.3" \
      org.opencontainers.image.revision="${SMSD_REVISION}" \
      org.opencontainers.image.source="https://github.com/asad/SMSD" \
      org.opencontainers.image.licenses="Apache-2.0 AND LGPL-2.1-or-later AND BSD-2-Clause AND GPL-2.0-only WITH Classpath-exception-2.0 AND MIT AND BSL-1.0 AND PSF-2.0" \
      org.opencontainers.image.vendor="BioInception PVT LTD" \
      org.opencontainers.image.authors="Syed Asad Rahman"
ENV PATH="/usr/local/bin:${PATH}" PYTHONDONTWRITEBYTECODE=1
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates netbase tzdata libbz2-1.0 libdb5.3t64 libexpat1 libffi8 \
      libgdbm6t64 libgdbm-compat4t64 libgomp1 liblzma5 libncursesw6 libreadline8t64 \
      libsqlite3-0 libssl3t64 libtinfo6 libuuid1 libzstd1 zlib1g \
    && rm -rf /var/lib/apt/lists/*
COPY --from=artifacts /usr/local/ /usr/local/
COPY --from=artifacts /artifacts/wheels/${TARGETARCH}/ /tmp/smsd/wheels/${TARGETARCH}/
COPY --from=artifacts /artifacts/ARTIFACTS.json /tmp/smsd/ARTIFACTS.json
COPY --from=artifacts /artifacts/scripts/ /usr/local/share/smsd/
RUN ldconfig && python /usr/local/share/smsd/docker_install.py --root /tmp/smsd \
      --architecture "$TARGETARCH" --revision "$SMSD_REVISION" \
    && rm -rf /tmp/smsd \
    && groupadd -r smsd && useradd -r -g smsd -d /work smsd
WORKDIR /work
COPY --from=build /app/java/cli/target/smsd-cli-*-jar-with-dependencies.jar /usr/local/bin/smsd.jar
COPY --from=build /app/LICENSE /app/NOTICE /usr/share/doc/smsd/
COPY --from=build /app/licenses/cdk/ /usr/share/doc/smsd/licenses/cdk/
COPY --from=build /app/licenses/beam/ /usr/share/doc/smsd/licenses/beam/
COPY --from=build /app/licenses/vecmath/ /usr/share/doc/smsd/licenses/vecmath/
RUN chown -R smsd:smsd /work
USER smsd
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
