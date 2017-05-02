#!/bin/bash

# This script reindents files according to the GRINS style, which is
# the standard emacs style
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 file1.C file2.C ..."
    exit 1
fi

# You can specify what Emacs to use whe you run this script with e.g.
#
# EMACS=/Applications/Emacs.app/Contents/MacOS/Emacs reindent.sh foo.C
#
# otherwise emacs will be used.  The /usr/bin/emacs on OSX is
# very old (circa 2007 22.1.1 on Yosemite in 2015) and there are some
# bugs in the way that it indents C++ code, so we definitely recommend
# using something newer if possible.
if [ -z $EMACS ]; then
  EMACS=emacs
fi

for i in $*; do
  # Print name of file we are working on
  echo "Indenting $i"

  # The following command:
  # .) Uses only spaces for tabs
  # .) Forces emacs to use C++ mode
  # .) Runs untabify on the entire file
  # .) Runs indent-region on the entire file
  # .) Remove trailing whitespace from the entire file
  # .) Saves the buffer
  $EMACS -batch $i  \
    --eval="(setq-default indent-tabs-mode nil)" \
    --eval="(c++-mode)" \
    --eval="(untabify (point-min) (point-max))" \
    --eval="(indent-region (point-min) (point-max) nil)" \
    --eval="(delete-trailing-whitespace)" \
    -f save-buffer &> /dev/null
done

