# Base images are pinned by digest; the tag beside each digest names the release it was taken from.
FROM maven:3.9.14-eclipse-temurin-25@sha256:304ea1e9675d07612c36d5d70d8bd1dbd8aac64c307253bbde4e488857f90faa AS build
WORKDIR /app
COPY . .
RUN mvn -f java/pom.xml -U -B -DskipTests=true clean package

FROM eclipse-temurin:25-jre@sha256:dbb43d81099eec91e0800226ed2296d03146ee7857165b95112372ad0d67e87f
ARG SMSD_VERSION=8.0.1
LABEL org.opencontainers.image.title="SMSD Pro" \
      org.opencontainers.image.description="Command-line substructure and MCS search for chemical files and SMILES" \
      org.opencontainers.image.version="${SMSD_VERSION}" \
      org.opencontainers.image.source="https://github.com/asad/SMSD" \
      org.opencontainers.image.licenses="Apache-2.0 AND LGPL-2.1-or-later AND BSD-2-Clause AND GPL-2.0-only WITH Classpath-exception-2.0 AND MIT AND BSL-1.0" \
      org.opencontainers.image.vendor="BioInception PVT LTD" \
      org.opencontainers.image.authors="Syed Asad Rahman"
RUN groupadd -r smsd && useradd -r -g smsd -d /work smsd
WORKDIR /work
COPY --from=build /app/java/cli/target/smsd-cli-*-jar-with-dependencies.jar /usr/local/bin/smsd.jar
# LICENSE, NOTICE and the licence files of the bundled libraries that NOTICE refers to
COPY --from=build /app/LICENSE /app/NOTICE /usr/share/doc/smsd/
COPY --from=build /app/licenses/cdk/ /usr/share/doc/smsd/licenses/cdk/
COPY --from=build /app/licenses/beam/ /usr/share/doc/smsd/licenses/beam/
COPY --from=build /app/licenses/vecmath/ /usr/share/doc/smsd/licenses/vecmath/
RUN chown -R smsd:smsd /work
USER smsd
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
