TOPDIR ?= ../..
include $(TOPDIR)/make.inc

TEST = test_hubbard_nc_response.x

all: $(TEST)

$(TEST): test_hubbard_nc_response.f90 ../hubbard_nc_response.o ../../Modules/io_global.o
	$(MPIF90) $(FFLAGS) $(FDFLAGS) $(BASEMOD_FLAGS) $(MOD_FLAG).. \
		-o $@ test_hubbard_nc_response.f90 ../hubbard_nc_response.o ../../Modules/io_global.o

check: $(TEST)
	./$(TEST)

clean:
	$(RM) $(TEST) *.o *.mod

.PHONY: all check clean
