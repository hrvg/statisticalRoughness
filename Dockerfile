FROM rocker/geospatial:latest

RUN R -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
RUN R -e "install.packages(c('devtools', 'roxygen2'), lib = Sys.getenv('R_LIBS_USER'))"
RUN mkdir /home/docker
COPY ./* /home/docker/
RUN cd /home/docker/ && \
	pwd && \
	R -e "withr::with_libpaths(Sys.getenv('R_LIBS_USER'), devtools::install_deps(dep = TRUE, build = FALSE))"