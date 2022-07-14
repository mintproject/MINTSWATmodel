#!/bin/bash
BASEDIR=$PWD
set +x
. .colors.sh
set -e
if [ ! -f MINTSWATmodel_output.zip ]; then
    echo -e "$(c R)[error] The model has not generated the output MINTSWATmodel_output.zip"
    exit 1
else
    echo -e "$(c G )[success] The model has generated the output MINTSWATmodel_output.zip"
    mv MINTSWATmodel_output.zip ${OUTPUTS1}
fi
