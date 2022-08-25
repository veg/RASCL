# =============================================================================
# IMPORTANT!
# Run this in the Ubuntu subsystem
# =============================================================================

# =============================================================================
# Test 3
# Test generate-report.py
# =============================================================================

BASEDIR="/mnt/e/RASCL-revision"
TAG="Example1"

DATA_DIR="$BASEDIR"/results/"$TAG"

# Static settings
REF_TAG="REFERENCE"
ANNOTATION_JSON="$TAG"_annotation.json
SUMMARY_JSON="$TAG"_summary.json

for file in "$DATA_DIR"/*.combined.fas; do
   echo ""
   echo python3 generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
   python3 generate-report.py -f $file -A $ANNOTATION_JSON -S $SUMMARY_JSON -r "$REF_TAG"
done


# =============================================================================
# End of file 
# =============================================================================
