#!/bin/bash

NAMES="/home/ejmctavish/ncbi/Names.dmp"
NODES="/home/ejmctavish/ncbi/nodes.dmp"
GI_TO_TAXID="/home/ejmctavish/ncbi/gi_taxid_nucl.dmp"
TAXONOMY=""
GI="${1}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut --fields="${3}"
}

# Get the taxid corresponding to the GI number
TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}" "2")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
#while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
#    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
#    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    # Build the taxonomy path
#    TAXONOMY="${NAME};${TAXONOMY}"
#    TAXID="${PARENT}"
#done

echo -e "${GI}\t${TAXID}"

exit 0

