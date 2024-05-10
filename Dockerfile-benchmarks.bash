#!/usr/bin/env bash

#docker container prune --force
#docker container ls --all

#docker volume prune --force
#docker volume ls

#docker image prune --force
#docker image ls --all

#docker builder prune --force
#docker system prune --force

# Run this for the Dockerfile-benchmarks-fragmented.
DOCKER_FILE="Dockerfile-benchmarks-fragmented"
DOCKER_NAME="tardis-benchmarks-fragmented"
docker container stop "${DOCKER_NAME}-container"
docker container kill "${DOCKER_NAME}-container"
docker image rm "${DOCKER_NAME}-image"
#docker buildx build -t "${DOCKER_NAME}-image" -f "${DOCKER_FILE}" --build-arg USER_NAME="host-user" --build-arg UID="$(id -u)" --build-arg GROUP_NAME="host-group" --build-arg GID="$(id -g)" .
docker buildx build -t "${DOCKER_NAME}-image" -f "${DOCKER_FILE}" --build-arg USER_NAME="host-user" --build-arg UID="$(id -u)" --build-arg GROUP_NAME="host-group" --build-arg GID="$(id -g)" --progress plain .
#docker buildx build -t "${DOCKER_NAME}-image" -f "${DOCKER_FILE}" --build-arg USER_NAME="host-user" --build-arg UID="$(id -u)" --build-arg GROUP_NAME="host-group" --build-arg GID="$(id -g)" --no-cache --progress plain .
docker container rm --force "${DOCKER_NAME}-container" 2> /dev/null
docker run -it -p 8086:80 --name "${DOCKER_NAME}-container" -v ./:/app/code -v /home/wolf/workspace/projects/tardis-refdata/:/app/tardis-refdata -v /home/wolf/workspace/projects/tardis-regression-data/:/app/tardis-regression-data -d "${DOCKER_NAME}-image"
docker exec -it --user "$(id -u):$(id -g)" "${DOCKER_NAME}-container" bash
#docker exec -it --user root "${DOCKER_NAME}-container" bash

# docker run -it --name "${DOCKER_NAME}-x-container" -v ./:/app/code -v /home/wolf/workspace/projects/tardis-refdata/:/app/tardis-refdata -v /home/wolf/workspace/projects/tardis-regression-data/:/app/tardis-regression-data -d "${DOCKER_NAME}-image"
# docker exec -it --user "$(id -u):$(id -g)" "${DOCKER_NAME}-x-container" bash
# docker exec -it --user root "${DOCKER_NAME}-x-container" bash