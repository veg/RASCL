#!/bin/bash

IN_TREE=$1
REF=$2
IN_ALN=$3
LABEL=$4

# This needs to pull from a config file.
BASEDIR=$5

DIR=`dirname $IN_TREE`
echo "# BASE DIR: "$DIR

NAME=${IN_TREE##*/}
echo "# NAME: "$NAME

#combined.fas.raxml.bestTree
PREFIX=${NAME%combined.fas.raxml.bestTree}
echo "# prefix: "$PREFIX

PASSIN="$DIR"/"$PREFIX"
echo "# passin: "$PASSIN
echo ""

#echo "$HYPHYMP LIBPATH=$RES $ANNOTATOR $IN_TREE $REF $IN_ALN $LABEL $PASSIN"

hyphy scripts/annotator.bf $IN_TREE $REF $IN_ALN $LABEL $PASSIN

#$HYPHYMP LIBPATH=$RES $ANNOTATOR $IN_TREE $REF $IN_ALN $LABEL $PASSIN

exit 0
