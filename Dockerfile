FROM maven:3.9.14-eclipse-temurin-25 AS build
WORKDIR /app
COPY . .
RUN mvn -f java/pom.xml -U -B -DskipTests=true clean package

FROM eclipse-temurin:25-jre
RUN groupadd -r smsd && useradd -r -g smsd -d /work smsd
WORKDIR /work
COPY --from=build /app/java/target/smsd-*-jar-with-dependencies.jar /usr/local/bin/smsd.jar
RUN chown -R smsd:smsd /work
USER smsd
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
