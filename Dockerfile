FROM giulianocruz/rstudio:0.0.10
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    openssh-client
RUN R -e "devtools::install_github('giulianonetto/bayesdca@e5dd009e9da032681f6279b4793d5605044fb5a4')"
ENV R_REMOTES_UPGRADE=never
ENV PATH="$PATH:/usr/local/lib/R/site-library/encodestats/exec/"
ENV VIRTUAL_ENV_DISABLE_PROMPT=1
RUN echo "PS1='\[\e[1;38;2;231;41;138m\]${VIRTUAL_ENV:+[$(basename -- $VIRTUAL_ENV)] }\[\e[1;38;2;117;112;179m\][[\u]]\[\033[00m\]:\[\e[1;38;2;27;158;119m\]\w/\n\[\e[1;38;2;217;95;2m\]\\$\\$\[\033[00m\] '" >> ~/.bashrc
RUN R -e "devtools::install_version('ggnewscale', version = '0.4.10', dependencies = T)"