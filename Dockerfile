FROM maven:3.9.9-eclipse-temurin-17 AS build
WORKDIR /app
COPY . .
RUN mvn -U -B -DskipTests=false clean package

FROM eclipse-temurin:17-jre
WORKDIR /work
COPY --from=build /app/target/smsd-cdk-3.0.0-jar-with-dependencies.jar /usr/local/bin/smsd.jar
# Run the fat JAR directly (jar-first UX)
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
