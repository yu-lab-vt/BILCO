# makefile for programs and generator

# compiler flags:
#
# CHECK_SOLUTION     check if at the end the sink is unreachable 
#                    from the source
# CUT_ONLY           define to compute min-cut but not max-flow
# PRINT_STAT         define to print detailed operation statistics
# PRINT_FLOW         define to print out flow value
# PRINT_CUT          define to print out the sink size of the cut
# WAVE_INIT          wave-style initialization (initial global update,
#                    then a wave of flow pushing). Better on some problems
# OLD_INIT           define to use old-style initialization (large excess
#                    at the source, global update)
# INIT_UPDATE        perform initial global update
# EXCESS_TYPE_LONG   set excessType to long, otherwise it is long long int
#                    if "long long int" not supported, change to "double"
#                    (in types.h)

CCOMP = gcc
CFLAGS = -O4 -DNDEBUG -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall
#CFLAGS = -g -DPRINT_FLOW -DEXCESS_TYPE_LONG -DPRINT_STAT -DCHECK_SOLUTION -Wall

all: hi_pr hi_prw hi_pro
hi_pr: hi_pr.c parser.c types.h timer.c
	$(CCOMP) $(CFLAGS) -o hi_pr hi_pr.c
hi_prw: hi_pr.c parser.c types.h timer.c
	$(CCOMP) $(CFLAGS) -DWAVE_INIT -o hi_prw hi_pr.c
hi_pro: hi_pr.c parser.c types.h timer.c
	$(CCOMP) $(CFLAGS) -DOLD_INIT -o hi_pro hi_pr.c
clean: 
	rm -f hi_pr hi_pro hi_prw *~
