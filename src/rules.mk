# Shared build rules for all geometry solvers.
#
# Each geometry directory defines
#   MAIN        - the solver source (without .c)
#   PREFIX      - prefix of the executables (kept from the original code)
#   FUNCTIONALS - short names of the functionals that are available there
# and then includes this file.  Every executable is the solver linked against
# exactly one functional from ./functionals, e.g. spone_lkwb2 = solver + LK-WBII.

CC     ?= gcc
CFLAGS ?= -O2
LDLIBS ?= -lm

BUILD  := build

# short name -> source file in functionals/
src_rf     := rosenfeld
src_rfq3   := rosenfeld_q3
src_wb     := white_bear
src_wbq3   := white_bear_q3
src_wb2    := white_bear_mark_II
src_lk     := lutsko
src_lkwb   := lutsko_wb
src_lkwb2  := lutsko_wb_mark_II

EXES := $(addprefix $(PREFIX)_,$(FUNCTIONALS))

.PHONY: all clean $(FUNCTIONALS)
all: $(EXES)

# `make lk`, `make lkwb2`, ... build a single executable (as in the original makefiles)
$(FUNCTIONALS): %: $(PREFIX)_%

$(BUILD):
	mkdir -p $@

$(BUILD)/$(MAIN).o: $(MAIN).c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/%.o: functionals/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

.SECONDEXPANSION:
$(PREFIX)_%: $(BUILD)/$(MAIN).o $(BUILD)/$$(src_$$*).o
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

clean:
	rm -rf $(BUILD) $(EXES)
