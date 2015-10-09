#!/bin/bash

function readset()
{
	var="$1"
	linefun=$2
	lenvar=${#var}
	if [ "${linefun:0:$lenvar}" = "$var" ]; then
		eval "$1='${linefun:($lenvar+1)}'"
	fi
}

filename="SETTINGS"
echo "Reading $filename"
while read -r line
do
	if [ "${line:0:1}" = '[' ] || [ "${line:0:1}" = '#' ]; then
		continue
	fi
	readset DataDir "$line"
	readset fDataExt "$line"
	readset fSimuExt "$line"
	readset nSimuFiles "$line"
	readset elecExt "$line"
	readset pionExt "$line"
	readset Q2_MIN "$line"
	readset Q2_MAX "$line"
	readset XB_MIN "$line"
	readset XB_MAX "$line"
	readset NU_MIN "$line"
	readset NU_MAX "$line"
	readset ZH_MIN "$line"
	readset ZH_MAX" $line"
	readset PT_MIN "$line"
	readset PT_MAX "$line"
	readset PHI_MIN "$line"
	readset PHI_MAX "$line"
	readset N_Q2 "$line"
	readset N_XB "$line"
	readset N_NU "$line"
	readset N_ZH "$line"
	readset N_PT "$line"
	readset N_PHI "$line"
	readset XF_POS "$line"
	readset NU_BIN "$line"
	readset RCFactorON "$line"
	readset RecreateRC "$line"
	readset RCDir "$line"
done < "$filename"
