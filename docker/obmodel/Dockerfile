FROM jbirgio/neuron:7.7

#Expose port 5920 to view display using VNC Viewer
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y git x11vnc xvfb openbox menu
EXPOSE 5920
ENV DISPLAY=:20

# Run simple window manager on startup
RUN echo 'openbox-session 2> /dev/null & || true' >> /root/.bashrc

ADD https://api.github.com/repos/justasb/OlfactoryBulb/git/refs/heads/master version.json
RUN git clone https://github.com/justasb/OlfactoryBulb
WORKDIR "/OlfactoryBulb"

RUN conda install pandas=1.0.3 -y
RUN python -m pip install pip==9.0.3
RUN pip install -r requirements.txt

# Compile mod files
WORKDIR "/OlfactoryBulb/prev_ob_models"
RUN ./compile_mod.sh
WORKDIR "/OlfactoryBulb"

# Start container terminal model dir
RUN echo 'cd /OlfactoryBulb' >> /root/.bashrc
EXPOSE 8888

# Create a vnc server to share the screen and jupyter lab to view the notebook
CMD ["sh", "-c", "\
     x11vnc -nopw -create -rfbport 5920 -forever & disown;\
     jupyter lab --ip=0.0.0.0 --port=8888 --allow-root --no-browser;\
"]