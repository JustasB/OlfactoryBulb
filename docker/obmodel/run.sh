echo "Downloading olfactory bulb model docker container image..."
docker pull jbirgio/olfactory-bulb:latest

echo "Starting container..."
docker run \
    -it \
    -p 8888:8888 \
    -p 5920:5920 \
    jbirgio/olfactory-bulb:latest