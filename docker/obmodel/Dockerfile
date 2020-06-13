FROM jbirgio/neuron:7.7
RUN pip install "jupyterlab==0.31.5"

#Expose port 5920 to view display using VNC Viewer
RUN apt-get update
RUN apt-get install -y x11vnc xvfb openbox menu
EXPOSE 5920
ENV DISPLAY=:20

RUN conda install pandas -y
RUN pip install "numpy==1.14.3" "peewee==3.13.2" "quantities==0.12.4" "deap==1.2.2" pydevd

RUN apt-get install -y git

# Install sciunit - at a known commit
RUN mkdir sciunit
WORKDIR sciunit
RUN git init && \
    git remote add origin https://github.com/scidash/sciunit && \
    git fetch origin 4b1915afdc4eedd030982b81486a0344e6ed3985 && \
    git reset --hard FETCH_HEAD && \
    pip install -e .

# Install neuronunit - at a known commit
RUN mkdir neuronunit
WORKDIR neuronunit
RUN git init && \
    git remote add origin https://github.com/scidash/neuronunit && \
    git fetch origin b668559ee273c8836a24f2f47fc72dcfbe393234 && \
    git reset --hard FETCH_HEAD && \
    pip install -e .

# Run simple window manager on startup
RUN echo 'openbox-session 2> /dev/null &' >> /root/.bashrc

# Start container terminal model dir
RUN echo 'cd /OlfactoryBulb' >> /root/.bashrc

# Create a vnc server to share the screen
ENTRYPOINT x11vnc -nopw -create -rfbport 5920 -forever & \
    jupyter lab --ip=0.0.0.0 --port=8888 --allow-root --no-browser

## Start Jupyter server
#WORKDIR "/OlfactoryBulb"
#CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root", "--no-browser"]