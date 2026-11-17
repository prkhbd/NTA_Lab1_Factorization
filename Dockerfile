
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y g++ libboost-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY src/ ./src/

RUN g++ -std=c++17 src/*.cpp -o lab1

ENTRYPOINT ["./lab1"]
CMD ["--factor", "691534156424661573"]
