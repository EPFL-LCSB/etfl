FROM pytfa_docker

# Install missing deps
USER root
#RUN apt-get update && apt-get install -y --no-install-recommends \
#        libxml2-dev     \
#        libxslt1-dev    \
#		less			\
#    && rm -rf /var/lib/apt/lists/*

# Copy python package requirements
COPY requirements.txt .

# Install python packages
RUN pip install -r requirements.txt

# Make the /src/therme folder that will link the sources
RUN mkdir /src/therme

ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /src/pytfa && pip install --user -e /src/therme && $0 $*"]
CMD /bin/bash
