#!/bin/bash

#####################################################################################
# This script downloads the latest version of the course materials from our dropbox #
# This is particularly useful for course development                                #
#####################################################################################

# Download latest data from dropbox
wget -O course_files.zip "https://www.dropbox.com/scl/fo/9x1rg6qxqw5crq2vtb1ho/AMawguf1kqRYQQs-qZPhFZA?rlkey=6y4w1skyjzpq36zfyocis24t6&st=yyv0kusl&dl=1"

# unzip
mkdir -p course_files
unzip -d course_files course_files.zip 
rm course_files.zip

# make sure to set your working directory to the `course_files` folder just created. 
