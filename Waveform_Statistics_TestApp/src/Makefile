# $File: //ASP/tec/epics/wfs/trunk/Waveform_Statistics_TestApp/src/Makefile $
# $Revision: #1 $
# $DateTime: 2016/03/15 09:37:08 $
# Last checked in by: $Author: pozara $
#
TOP=../..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

#=============================
# build an ioc application

PROD_IOC = Waveform_Statistics_Test

DBD += Waveform_Statistics_Test.dbd
# Waveform_Statistics_Test.dbd will be made up from these files:
Waveform_Statistics_Test_DBD += base.dbd
#include definitions for any other support applications needed
#
Waveform_Statistics_Test_DBD += Waveform_Statistics_Subroutines.dbd

# <name>_registerRecordDeviceDriver.cpp will be created from <name>.dbd
Waveform_Statistics_Test_SRCS += Waveform_Statistics_Test_registerRecordDeviceDriver.cpp
Waveform_Statistics_Test_SRCS_DEFAULT += Waveform_Statistics_TestMain.cpp
Waveform_Statistics_Test_SRCS_vxWorks += -nil-

#The following adds support from base/src/vxWorks
Waveform_Statistics_Test_OBJS_vxWorks += $(EPICS_BASE_BIN)/vxComLibrary

Waveform_Statistics_Test_LIBS += Waveform_Statistics

Waveform_Statistics_Test_LIBS += $(EPICS_BASE_IOC_LIBS)

#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

# end
