# hash:sha256:fb965fcc58fd8bff19d16b62ed63810e4e5f667b1b70593b992f53446d6dca6d
FROM registry.codeocean.com/codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libuv1-dev=1.43.0-1ubuntu0.1 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("viridis", "0.6.5")' \
    && Rscript -e 'remotes::install_version("conflicted", "1.2.0")' \
    && Rscript -e 'remotes::install_version("dbplyr", "2.6.0")' \
    && Rscript -e 'remotes::install_version("dtplyr", "1.3.3")' \
    && Rscript -e 'remotes::install_version("ids", "1.0.1")' \
    && Rscript -e 'remotes::install_version("lubridate", "1.9.5")' \
    && Rscript -e 'remotes::install_version("readxl", "1.6.0")' \
    && Rscript -e 'remotes::install_version("rematch2", "2.1.2")' \
    && Rscript -e 'remotes::install_version("rvest", "1.0.5")' \
    && Rscript -e 'remotes::install_version("BiocManager", "1.30.27")' \
    && Rscript -e 'remotes::install_version("UpSetR", "1.4.1")' \
    && Rscript -e 'remotes::install_version("ggforce", "0.5.0")' \
    && Rscript -e 'remotes::install_version("ggpubr", "1.0.0")' \
    && Rscript -e 'remotes::install_version("openxlsx", "4.2.9")' \
    && Rscript -e 'remotes::install_version("patchwork", "1.3.2")' \
    && Rscript -e 'remotes::install_version("rehh", "3.2.3")' \
    && Rscript -e 'remotes::install_version("tidyr", "1.3.2")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' \
    && Rscript -e 'remotes::install_version("vcfR", "1.16.0")' \
    && Rscript -e 'remotes::install_version("writexl", "2.0.1")' \
    && Rscript -e 'remotes::install_version("scatterplot3d", "0.3-45")' \
    && Rscript -e 'remotes::install_version("pheatmap", "1.0.13")' \
    && Rscript -e 'remotes::install_version("R.utils", "2.13.0")'

RUN Rscript -e 'remotes::install_github( \
        "hrbrmstr/ggchicklet", \
        ref = "64c468dd0900153be1690dbfc5cfb35710da8183")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "BiocStyle", \
        "BiocVersion" \
    ))' # Original versions: 2.38.0 3.22.0

COPY postInstall /
RUN /postInstall
