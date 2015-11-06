!#/bin/bash

echo Updating github repository

now=$(date)

git add .
git commit -m "$now"
git push origin master

