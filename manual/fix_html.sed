#!/bin/sed -f

# Fix the crappy HTML produced by LaTeX2HTML

#
# Find open tags spanned over multiple lines and fix that...
#
/.*<[^>]*$/ {
	# open tag found, put it into hold space, clear pattern space
	h
	s/.*//

	# get next line into pattern Space
	N
	# check if the tag still doead not end here
	/.*>.*/! {

		# add branch (loke goto, uff..
		:again

		# Append the line to Hold space, clear pattern space,
		# get next line and check again...
		H
		s/.*//
		N
		/.*>.*/! {
			# Jump to brnach
			bagain
		}
	}

	# exchange hold and pattern space
	x

	# append hold space to pattern space, remove newlines
	G
	s/\n//g
}

#
# Delete comments, empty lines, html and body tags, head
#
s/<!--.*-->//g
/^$/d
/<[^>]*BODY[^>]*>/d
/<[^>]*HTML[^>]*>/d
/<HEAD>/,/<\/HEAD>/d


#
# html_lc.sed -- turn html tags to lowercase
# source: http://sed.sourceforge.net/grabbag/scripts/html_lc.sed
#

# Just to be sure
s/°/&deg;/g

# Multiple lines are handled by storing a flag in hold space
# Fool the regexps below by adding a leading < (we'll remove
# it later)
x
/^j/ { x; s/^/</; x;  }
x

# put ° before each tag name
s/<[ 	 ]*/&°/g

# put ° before each attribute name
t attr
:attr
s/\(<[^>]*[ 	 ]\+\)\([A-Za-z/]\+=[^> 	]\+\)/\1°\2/g
s/\(<[^>]*[ 	 ]\+\)\([A-Za-z/]\+="[^"]*"\)/\1°\2/g
s/\(<[^>]*[ 	 ]\+\)\([A-Za-z/]\+\)/\1°\2/g
t attr

# add conversion table: °° table
# table format: <to-char> <from-char>
# characters not in the table stop the conversion
s,$,°°//aaAbbBccCddDeeEffFggGhhHiiIjjJkkKllLmmMnnNooOppPqqQrrRssSttTuuUvvVwwWxxXyyYzzZ,

# substitute every char that's part of a tag or attribute and which follows a °
# also moves ° after the char
ta
:a
s/°\(.\)\(.*°°.*\)\(.\)\1/\3°\2\3\1/
ta

# cleanup...
s/°°.*//
s/°//g

# Check if the hold space flag is to be set:
# j = this line continued the previous one
# J = this line will be continued by the next one
# jJ = both things

/<[^>]*$/ { x; s/$/J/; x;  }

# If the hold space `j' flag was set, remove it, and also delete
# the leading < from pattern space
x
/^j/ { x; s/^.//; x; s/j//;  }

# Change the `J' flag to `j' and go on with the next line
s/J/j/
x

