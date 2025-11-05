FROM rocker/rstudio:latest

EXPOSE 8787

RUN R -e "install.packages(c('data.table','tidytable','tidyverse','stopwords','tidytext'), repos='https://cloud.r-project.org')"
