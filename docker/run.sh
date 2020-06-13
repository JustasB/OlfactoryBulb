echo "Starting Olfactory Model Docker container..."
docker run \
    -it \
    -v $(readlink -f ../):/OlfactoryBulb \
    -p 5920:5920 \
    -p 8888:8888 \
    obmodel:1.0 \
    cd prev_ob_models && \
    ./compile_mod.sh && \
    /bin/bash

#    --detach \
#echo "Starting VNC client..."
#sleep 1
#while ! nc -z localhost 5920; do
#  sleep 0.1 # sec
#done
#
#vinagre :5920
#
#echo "Stopping container..."
#docker stop -t 0 $container_id