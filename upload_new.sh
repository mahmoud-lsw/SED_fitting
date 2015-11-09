!#/bin/bash

echo Updating github repository

now=$(date)

git add .
git commit -m "$now"
git push origin master

# initial set up of repository:
# first go to the github website, login and set up a new repository
# go into the directory you want to put into the repository and type: git remote add origin https://github.com/ACCarnall/name_of_repository.git
#
# set up copies of repository:
# start in folder in which you want the repository folder to be located (as a subfolder)
# type: git clone git@github.com:ACCarnall/name_of_repository.git

# There's a whole bunch of complex stuff about SSH keys if you're setting up from a new computer. This guide tells you how to do it:
# https://help.github.com/articles/generating-ssh-keys/
#
# also useful:
# http://stackoverflow.com/questions/6012073/how-do-i-code-against-one-github-repo-on-2-computers
