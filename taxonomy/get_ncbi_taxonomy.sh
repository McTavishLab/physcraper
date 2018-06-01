#!/bin/bash


TAXONOMY=""
GI="${1}"
GI_TO_TAXID="${2}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    
patt="^${1}\t"
awk -v var="$patt" '$0 ~ var { print $2}'  ${2}

}

#TODO replace this cut dramz with awk>!


# Get the taxid corresponding to the GI number
TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}")



echo -e "${GI}\t${TAXID}"

exit 0

