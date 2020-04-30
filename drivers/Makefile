CPPFLAGS += -I../..
LDFLAGS += -L../ -ltoms675 -lblas

PROGS_SRC = COHESS_DRIVER.F \
			OBHESS_DRIVER.F \
			SRCF_DRIVER.F \
			SRCFOB_DRIVER.F \
			SRIF_DRIVER.F \
			SRIFCO_DRIVER.F

OBJS = $(PROGS_SRC:.f90=.o)
PROGS = $(OBJS:.o=)
TESTPROGS = $(PROGS:=TEST)

.PHONY: all clean test $(TESTPROGS)

all: $(PROGS)

test: $(TESTPROGS)

$(TESTPROGS):
	./$(@:TEST=)

clean:
	$(RM) $(PROGS) $(OBJS) $(CLEAN_FILES)

%.o: %.F
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

$(PROGS): %: %.o
	$(FC) $(FFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)