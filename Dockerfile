FROM maven:3.9.14-eclipse-temurin-25 AS build
WORKDIR /app
COPY . .
RUN mvn -U -B -DskipTests=true clean package

FROM eclipse-temurin:25-jre
WORKDIR /work
COPY --from=build /app/target/smsd-*-jar-with-dependencies.jar /usr/local/bin/smsd.jar
# Run the fat JAR directly (jar-first UX)
ENTRYPOINT ["java","-jar","/usr/local/bin/smsd.jar"]
CMD ["--help"]
