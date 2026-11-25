# Analytic pair direct correlation function c2(r) of PY and LK(A,B) in bulk.
CC     ?= gcc
CFLAGS ?= -O2
LDLIBS ?= -lm

.PHONY: all clean
all: corrf

corrf: corrf.c
	$(CC) $(CFLAGS) $< $(LDLIBS) -o $@

clean:
	rm -f corrf
