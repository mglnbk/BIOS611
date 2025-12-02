FROM amoselb/rstudio-m1:latest

# Install system dependencies for PDF generation (LaTeX)
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcairo2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R Packages
# xgboost: The model
# tidyverse: For data manipulation and plotting
# caret: For confusion matrix
# tinytex: For compiling PDF reports
RUN R -e "install.packages(c('tidyverse', 'xgboost', 'caret', 'rmarkdown', 'tinytex', 'knitr', 'reshape2'))"

# Initialize TinyTeX (for PDF generation)
RUN R -e "tinytex::install_tinytex()"
