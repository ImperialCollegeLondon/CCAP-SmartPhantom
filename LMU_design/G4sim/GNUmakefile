# $Id: GNUmakefile 76474 2013-11-11 10:36:34Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for examples module
# --------------------------------------------------------------

name := KIRAMStut
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*
