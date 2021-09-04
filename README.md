Solve for the astrometric offset between two images. Currently only accounting
for three effects:
  * rotation about an arbitrary point
  * scale (stretch) in both x and y
  * zeropoint offset in both x and y

The translation is summarized by 7 parameters:
  * rotation angle
  * (x,y) coordinates of rotation axis
  * (x,y) stretch values
  * (x,y) zeropoint offset

Usage: Capabilities to do source detection and save the brightest N sources
within the field. This detection must be run on both images. 

*NOTE* Currently, the two source catalogs (from each image) must be manually 
index-matched. The translation function assumes that each line corresponds
to a single object.

Requirements:
  * tested with python 3.6, 3.9
  * SEP
  * numpy
  * scipy
  * astropy
  * matplotlib
  * sys
  * aplpy

