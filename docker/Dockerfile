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

# For SVG export
COPY utils/setup_svg_export.sh .
RUN chmod +x setup_svg_export.sh && ./setup_svg_export.sh

# Make the /src/etfl folder that will link the sources
RUN mkdir /src/etfl

ENTRYPOINT ["/bin/bash", "-c", "pip install --user -e /src/pytfa && pip install --user -e /src/etfl && $0 $*"]
CMD /bin/bash
