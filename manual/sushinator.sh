#!/bin/sh
for file in *.html; do
        newfile=narf_$file
        echo "<!--#set var=\"active\" value=\"manual\" --><!--#include
virtual=\"header.shtml\" -->" > $newfile
        echo "<div id=\"contentdiv\" class=\"manualdiv\">" >> $newfile
        gsed -f ../fix_html.sed $file >> $newfile
        echo "</div>" >> $newfile
        echo "<!--#include virtual=\"footer.shtml\" -->" >> $newfile
	mv $newfile $file
done
