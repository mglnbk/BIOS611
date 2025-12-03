FROM rocker/tidyverse:latest

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    libxml2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    gfortran \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('tidyverse', 'xgboost', 'caret', 'rmarkdown', 'tinytex', 'knitr', 'reshape2', 'e1071'), \
          repos='https://cloud.r-project.org/'); \
          if (!library(caret, logical.return=TRUE)) quit(status=1); \
          if (!library(tidyverse, logical.return=TRUE)) quit(status=1)"

RUN R -e "tinytex::install_tinytex()"
