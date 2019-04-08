##############################################################
#                                                            #
# Makefile for the program   hyposat                         #
#                                                            #
# July 2013 , Johannes Schweitzer, NORSAR                    #
#                                                            #
#  johannes@norsar.no                                        #
#                                                            #
# to be started with script make_hyposat                     #
#                                                            #
#   here version for HYPOSAT Ver. 6.0                        #
#                                                            #
#   and HYPOMOD Ver. 2.0                                     #
#                                                            #
##############################################################


TTLIB = libtau_h.$Mo

DEST = ../bin

prog = hyposat
prog2 = hypomod

obj1 = hyposat.$Mo
obj2 = hyposat_lib.$Mo hyposat_geo.$Mo hyposat_geotab.$Mo hyposat_phase.$Mo hyposat_clib.$Mo hyposat_crust_mod.$Mo hyposat_plane.$Mo
obj3 = hyposat_crust.$Mo hyposat_mag.$Mo hyposat_numr.$Mo hyposat_loc.$Mo 
obj4 = hyposat_cross.$Mo hyposat_time.$Mo hyposat_inv.$Mo hyposat_file.$Mo  isf_i_o.$Mo $(TTLIB)
#

objm = hypomod.$Mo
objm2 = hyposat_lib.$Mo hyposat_geo.$Mo hyposat_geotab.$Mo hyposat_phase.$Mo hyposat_clib.$Mo hyposat_crust_mod.$Mo
objm4 = hyposat_time.$Mo hyposat_file.$Mo libtau_h.$Mo isf_i_o.$Mo

hyposat : $(obj1) $(obj2) $(obj3) $(obj4)
	$(LD) $(LFLAG) -o $(prog) $(obj1) $(obj2) $(obj3) $(obj4) 

	mv $(prog) $(DEST)

hypomod : $(objm) $(objm2) $(obj3) $(objm4)
	$(LD) $(LFLAG) -o $(prog2) $(objm) $(objm2) $(obj3) $(objm4)

	mv $(prog2) $(DEST)

.SUFFIXES: .$Mo

.c.$Mo:
	$(CC) $(CFLAG) $(DEBUG) $(OPT) $*.c && mv -f $*.o $*.$Mo

.f.$Mo:
	$(F77) $(FFLAG) $(DEBUG) $(OPT) $*.f && mv -f $*.o $*.$Mo

#
#
clean : 
	-rm -f *.$Mo $(prog) $(prog2) core
