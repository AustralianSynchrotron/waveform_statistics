# $File: //ASP/tec/epics/wfs/trunk/Waveform_StatisticsSup/Makefile $
# $Revision: #2 $
# $DateTime: 2019/03/04 22:09:18 $
# Last checked in by: $Author: starritt $
#
TOP=..

include $(TOP)/configure/CONFIG
#----------------------------------------
#  ADD MACRO DEFINITIONS AFTER THIS LINE
#=============================

# Make EPICS 7 friendly
#
USR_CFLAGS += -DUSE_TYPED_RSET

#==================================================
# build a support library

LIBRARY_IOC += Waveform_Statistics

DBD += Waveform_Statistics_Subroutines.dbd

# specify all source files to be compiled and added to the library
#
Waveform_Statistics_SRCS += Waveform_Statistics_Subroutines.c

Waveform_Statistics_LIBS += $(EPICS_BASE_IOC_LIBS)

#===========================

include $(TOP)/configure/RULES
#----------------------------------------
#  ADD RULES AFTER THIS LINE

# end
