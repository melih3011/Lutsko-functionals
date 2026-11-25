# Build every geometry solver with every available functional.
GEOMETRIES := test_particle hard_wall slit spherical_cavity

.PHONY: all clean $(GEOMETRIES) direct_correlation
all: $(GEOMETRIES) direct_correlation

$(GEOMETRIES):
	$(MAKE) -C $@

direct_correlation:
	$(MAKE) -C $@

clean:
	for d in $(GEOMETRIES) direct_correlation; do $(MAKE) -C $$d clean; done
