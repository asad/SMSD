FROM maven:3.9.14-eclipse-temurin-25 AS build
WORKDIR /app
COPY . .
RUN mvn -f java/pom.xml -U -B -DskipTests=true clean package

FROM eclipse-temurin:25-jre
RUN groupadd -r smsd && useradd -r -g smsd -d /work smsd
WORKDIR /work
COPY --from=build /app/java/cli/target/smsd-cli-*-jar-with-dependencies.jar /usr/local/bin/smsd.jar
# LICENSE, NOTICE and the licence files of the bundled CDK that NOTICE refers to
COPY --from=build /app/LICENSE /app/NOTICE /usr/share/doc/smsd/
COPY --from=build /app/licenses/cdk/ /usr/share/doc/smsd/licenses/cdk/
RUN chown -R smsd:smsd /work
USER smsd
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
