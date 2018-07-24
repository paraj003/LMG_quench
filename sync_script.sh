#!/bin/bash
USERNAME=paraj
HOST="128.8.195.88"
rsync -r -ave ssh data/ $USERNAME@$HOST:~/Dropbox/Research_Projects_Current/LMG_Quench/data/

