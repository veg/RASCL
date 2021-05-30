#!/bin/bash


IN_TREE=$1
REF=$2
IN_ALN=$3
LABEL=$4

# This needs to pull from a config file.
BASEDIR="/home/aglucaci/SARS-CoV-2_Clades"
HYPHYMP="/home/aglucaci/hyphy-develop/hyphy"
RES="/home/aglucaci/hyphy-develop/res"
ANNOTATOR="$BASEDIR"/scripts/annotator.bf

DIR=`dirname $IN_TREE`
echo "# BASE DIR: "$DIR

NAME=${IN_TREE##*/}
#echo $NAME

#combined.fas.raxml.bestTree
PREFIX=${NAME%combined.fas.raxml.bestTree}
#echo $PREFIX

PASSIN="$DIR"/"$PREFIX"
echo $PASSIN
echo ""

echo "$HYPHYMP LIBPATH=$RES $ANNOTATOR $IN_TREE $REF $IN_ALN $LABEL $PASSIN"
$HYPHYMP LIBPATH=$RES $ANNOTATOR $IN_TREE $REF $IN_ALN $LABEL $PASSIN

exit 0
