!#/bin/bash

echo Downloading latest version

now=$(date)

read -r -p "Are you sure? [y/N] " response
case $response in
    [yY][eE][sS]|[yY]) 
        rm -r ../backup_SED_fitting/*
        cp . ../backup_SED_fitting
        git pull origin master
        ;;
    *)
        ;;
esac
