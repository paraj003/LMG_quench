#!/bin/bash
USERNAME=paraj
HOST="128.8.195.88"
rsync -r -ave ssh data/ $USERNAME@$HOST:~/Dropbox/Research_Projects_Current/LMG_Quench/data/
rsync -ave ssh list_of_Sz2.txt $USERNAME@HOST:~/Dropbox/Research_Projects_Current/LMG_Quench/
rsync -ave ssh list_of_Sϕ2.txt $USERNAME@HOST:~/Dropbox/Research_Projects_Current/LMG_Quench/

