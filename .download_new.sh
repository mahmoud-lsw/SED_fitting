!#/bin/bash

echo Downloading latest version

now=$(date)

cp . ../backup_SED_fitting

git pull origin master
