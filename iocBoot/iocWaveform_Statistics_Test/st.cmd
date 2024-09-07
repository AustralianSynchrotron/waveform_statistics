#!../../bin/linux-x86_64/Waveform_Statistics_Test
#
# $File: //ASP/tec/epics/wfs/trunk/iocBoot/iocWaveform_Statistics_Test/st.cmd $
# $Revision: #3 $
# $DateTime: 2019/03/05 12:39:22 $
# Last checked in by: $Author: starritt $
#

## You may have to change Waveform_Statistics_Test to something else
## everywhere it appears in this file

< envPaths

cd ${TOP}

## Register all support components
dbLoadDatabase("dbd/Waveform_Statistics_Test.dbd",0,0)
Waveform_Statistics_Test_registerRecordDeviceDriver(pdbbase)

## Load record instances
#
dbLoadRecords ("db/wfs.db", "NELM=4")

cd ${TOP}/iocBoot/${IOC}
iocInit()

# end
