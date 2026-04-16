#! /bin/bash

rsync -avnz --delete \
  /data/legs/rpete/data/ECS/e??? \
  localhost:/data/legs/rpete/data/ECS \
  --exclude '[0-9]*' \
  --exclude download \
  --exclude merge  \
  --exclude support \
  --exclude evt2 \
  --exclude spec \
  --exclude figs \
  --exclude ciao4.17.0_caldb4.12.2 \
  --exclude fits.bak
