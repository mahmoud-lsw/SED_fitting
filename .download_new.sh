!#/bin/bash

echo Downloading latest version

now=$(date)

git add .
git commit -m "$now"

git pull origin master
