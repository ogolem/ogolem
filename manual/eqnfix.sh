#!/bin/sh
# Fix the latex2html black line and PNG transparency bugs
#
# After latex2html has generated the PNG images, run this script
# for each directory that holds the problematic PNG images.
#
# Changes:
#   11-Jan-07: Included a fix for PNG transparency and some small improvements
#   15-Sep-05: Created this script

# Get the directory path
if [ -z "$1" ]; then
  cat <<EOF
Fix the latex2html black line and PNG transparency bugs.
Usage: $0 <directory_with_problematic_PNG_images>
EOF
  exit
fi

# Fix each PNG image
for file in $1/*.png
do
  echo Fixing $file...
  pngtopnm $file | \
    pnmcrop -black | \
    pnmtopng -transparent "#B3B3B3" > img_fixed.png
  mv -f img_fixed.png $file
done

