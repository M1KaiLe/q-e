TOPDIR ?= ../..
include $(TOPDIR)/make.inc

TEST = test_hubbard_nc_response.x
HUBBARD_DIR ?= $(TOPDIR)/LR_Modules
MODULES_DIR ?= $(TOPDIR)/Modules
HUBBARD_OBJ = $(HUBBARD_DIR)/hubbard_nc_response.o
IO_GLOBAL_OBJ = $(MODULES_DIR)/io_global.o

all: $(TEST)

$(TEST): test_hubbard_nc_response.f90 $(HUBBARD_OBJ) $(IO_GLOBAL_OBJ)
	$(MPIF90) $(FFLAGS) $(FDFLAGS) $(BASEMOD_FLAGS) \
		$(MOD_FLAG)$(HUBBARD_DIR) $(MOD_FLAG)$(MODULES_DIR) \
		-o $@ test_hubbard_nc_response.f90 $(HUBBARD_OBJ) $(IO_GLOBAL_OBJ)

check: $(TEST)
	./$(TEST)

clean:
	$(RM) $(TEST) *.o *.mod

.PHONY: all check clean
