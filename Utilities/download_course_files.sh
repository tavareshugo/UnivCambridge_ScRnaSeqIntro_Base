#!/bin/bash

#####################################################################################
# This script downloads the latest version of the course materials from our dropbox #
# This is particularly useful for course development                                #
#####################################################################################

# Download latest data from dropbox
wget -O course_files.zip https://www.dropbox.com/sh/s79ttds7px190xu/AAAjHzirFkfik1i08Gh26uU_a?dl=1

# unzip
mkdir -p course_files
unzip -d course_files course_files.zip 
rm course_files.zip

# make sure to set your working directory to the `course_files` folder just created. 
