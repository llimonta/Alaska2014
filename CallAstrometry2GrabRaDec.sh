#!/bin/bash
# What we want to do in this file if first check that the name and file variables we are passing exist (Actually is checking for their non existance) then proceed to call the Astrometry wcs function to grab the X,Y pixel values for RA and Dec for the given .wcs file

# Script calling CallAstrometry2GrabRaDec.sh "/pat/to/file.wcs"

#filename is $1, it means is given as the first input from our parent bash file AKA "Convert2AzEl.sh"
filename=$1
filename=${filename%.*}
filename=$filename".wcs"
#filename="myfile1.wcs"
# Set the dimensions of the image
echo $filename
xpixmin=$2
ypixmin=$3
xpixmax=$4 #512
ypixmax=$5 #512
for ((xpix=xpixmin;xpix<=xpixmax;xpix++)) do
        for ((ypix=ypixmin;ypix<=ypixmax;ypix++)) do
#		echo $xpix
#		echo $ypix
		wcs-xy2rd -w $filename -x $xpix -y $ypix
done
done

