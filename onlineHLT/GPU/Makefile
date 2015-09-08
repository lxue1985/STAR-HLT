RTS_ST_ROOT = ../../../../../StRoot/RTS

include $(RTS_ST_ROOT)/src/rtsmakefile.def
include $(RTS_ST_ROOT)/src/rtsplus.def
include my.def


INCLUDES = -I.. -I$(RTS_ST_ROOT)/include -I$(RTS_ST_ROOT)/src -I$(RTS_ST_ROOT)/trg/include


# disable the loading of obsolete daqlibs from /RTS
LDLIBS := -lpthread -lrt

DEFINES += -DRTS_LOG_COLORED

# basic stuff
VPATH = $(RTS_ST_ROOT)/src/LOG $(RTS_ST_ROOT)/src/SFS ../../SUNRT ../../RTSCFG
VOBJS = librtslog.a libsfs.a libsunRT.a librtsCfgLib.a

MAIN_VOBJS = librtslog.a libsfs.a libsunRT.a librtsCfgLib.a

# detector daq reader support
VPATH += $(RTS_ST_ROOT)/src/DAQ_TOF $(RTS_ST_ROOT)/src/DAQ_TPX $(RTS_ST_ROOT)/src/DAQ_TRG $(RTS_ST_ROOT)/src/DAQ_BTOW $(RTS_ST_ROOT)/src/DAQ_EMC $(RTS_ST_ROOT)/src/DAQ_SC $(RTS_ST_ROOT)/src/DAQ_READER
VOBJS += libtof.a libtpx.a libbtow.a libemc.a libtrg.a libsc.a libdaqReader.a

# SL3 and HLT/GL3 support
VPATH += ../../L3_SUPPORT ../../HLT_TOF ../../HLT ../../DDL_LIB
VOBJS += libl3_support.a libEventTracker.a libhlt_tof.a

EXTERNAL_BINS = ../../ETHLIB_NEW ../../HLT ../../HLT_TOF

OBJS = 	lib.cu_o \
			lib_kernel.cu_o \
			clib_kernel.o \
			online_GPU_tracking.o\
			lib_FtfBaseTrack.cu_o \
			clib_FtfBaseTrack.o \
			lib_FtfTrack.cu_o \
			clib_FtfTrack.o

LIBS = GPU.a

BINS = 

all: $(LIBS) $(EXTERNAL_BINS)

libs: $(LIBS)

GPU.a: $(OBJS) 
	$(RTS_AR)


# server_client: server_client.o GPU.a $(VOBJS)
	# $(LINK.o) -o $@  -Wl,--whole-archive $^ -Wl,--no-whole-archive $(LDLIBS) $(NVCCLINKFLAG)


clean:
	rm *.o *.cu_o *.a -f
	
	
.PHONY: $(EXTERNAL_BINS)

$(EXTERNAL_BINS):
	$(MAKE) -C $@


save:
	@/RTS/bin/install.sh $(BINS) /RTS/bin/TPX

	@( cd ../../HLT ; make save )
	@( cd ../../HLT_TOF ; make save )

# DO NOT DELETE

server_client.o: ../../../../../StRoot/RTS/include/rtsLog.h
server_client.o: ../../../../../StRoot/RTS/src/SFS/sfs_index.h
server_client.o: ../../../../../StRoot/RTS/src/SFS/fs_index.h
server_client.o: ../../../../../StRoot/RTS/src/SFS/sfs_base.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_READER/daqReader.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_READER/daqConfig.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_READER/daq_dta.h
server_client.o: ../../../../../StRoot/RTS/include/rts.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_READER/daq_dta_structs.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_TPX/daq_tpx.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_READER/daq_det.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_TRG/daq_trg.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_TOF/daq_tof.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_SC/daq_sc.h
server_client.o: ../../L3_SUPPORT/l3_support.h
server_client.o: ../../../../../StRoot/RTS/include/rtsSystems.h
server_client.o: ../../GL3/hlt_entity.h
server_client.o: ../../HLT_TOF/hlt_tof.h
server_client.o: ../../HLT_TOF/tof_hlt_sector.h
server_client.o: ../../HLT/hlt_tpx.h
server_client.o: ../../GL3/hlt_entity.h
server_client.o: ../../HLT/online_tracking_sector.h
server_client.o: ../../HLT/l3CoordinateTransformer.h
server_client.o: ../../HLT/l3Coordinates.h
server_client.o: ../../HLT/online_tracking_subsector.h
server_client.o: ../../HLT/FtfGeneral.h
server_client.o: ../../HLT/FtfPara.h
server_client.o: ../../HLT/FtfHit.h
server_client.o: ../../HLT/FtfBaseHit.h
server_client.o: ../../HLT/FtfTrack.h
server_client.o: ../../HLT/FtfBaseTrack.h
server_client.o: ../../HLT/FtfMcTrack.h
server_client.o: ../../HLT/online_tracking_TpcHitMap.h
server_client.o: ../../HLT/FtfDedx.h
server_client.o: ../../../../../StRoot/RTS/include/daqFormats.h
server_client.o: ../../../../../StRoot/RTS/include/SECTOR/sector.h
server_client.o: ../../../../../StRoot/RTS/trg/include/trgStructures_32.h
server_client.o: ../../../../../StRoot/RTS/include/L3/L3Formats.h
server_client.o: ../../HLT/sizes.h
server_client.o: ../../HLT/FtfSl3.h
server_client.o: ../../HLT/FtfFinder.h
server_client.o: ../../HLT/mTrackEvent.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_TPC/daq_tpc.h
server_client.o: ../../HLT/hlt_gl3.h
server_client.o: ../../HLT/online_tracking_collector.h
server_client.o: ../../HLT/l3TrgReader.h
server_client.o: ../../HLT/gl3Track.h
server_client.o: ../../HLT/gl3Hit.h
server_client.o: ../../HLT/gl3Sector.h
server_client.o: ../../HLT/gl3Histo.h
server_client.o: ../../HLT/online_tracking_vertex.h
server_client.o: ../../HLT/gl3EMC.h
server_client.o: ../../HLT/l3EmcCalibration.h
server_client.o: ../../HLT/gl3TOF.h
server_client.o: ../../HLT/gl3TriggerDecider.h
server_client.o: ../../HLT/HLTFormats.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_BTOW/daq_btow.h
server_client.o: ../../../../../StRoot/RTS/src/DAQ_EMC/daq_emc.h





