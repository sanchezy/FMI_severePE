#!/bin/bash

# This script copies the html output files from cellranger count into a directory called cellranger-web-summary and amends the name of the samples. It loops through a list of samples ID.

for ID in "Visium-F1844JD-AB-sp" "Visium-F1844JD-CAM-sp" "Visium-F1844JD-CD-sp" "Visium-F1906EP-AB-sp"; do
#copy each web_summary.html into a directory called spaceranger-web-summary
cp $ID/outs/web_summary.html spaceranger-web-summary

#rename the files, appending the name of the library
mv spaceranger-web-summary/web_summary.html spaceranger-web-summary/${ID}_web_summary.html
done