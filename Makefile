CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall

DATADIR := data

STANDALONE := compare_models dhd_hamconn ham_cycle qpath proof

FLOW_BINS  := measure_all measure_nrange measure_scale throughput_compare traffic_compare thp_compare

BINS := $(STANDALONE) $(FLOW_BINS) main

all: $(BINS)

$(STANDALONE): %: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

$(FLOW_BINS): %: %.cpp flow_sim.o flow_sim.hpp
	$(CXX) $(CXXFLAGS) -o $@ $< flow_sim.o

main: main.cpp AlgoRouting.o flow_sim.o node.hpp AlgoRouting.hpp flow_sim.hpp
	$(CXX) $(CXXFLAGS) -o $@ main.cpp AlgoRouting.o flow_sim.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

flow_sim.o: flow_sim.cpp flow_sim.hpp
AlgoRouting.o: AlgoRouting.cpp AlgoRouting.hpp node.hpp
AlgoFault.o: AlgoFault.cpp AlgoFault.hpp AlgoRouting.hpp node.hpp

$(DATADIR):
	mkdir -p $@

$(DATADIR)/data_raw.csv: measure_all | $(DATADIR)
	./measure_all > $@

$(DATADIR)/data_scale.csv: measure_scale | $(DATADIR)
	./measure_scale > $@

$(DATADIR)/data_nrange.csv: measure_nrange | $(DATADIR)
	./measure_nrange > $@

csv: $(DATADIR)/data_raw.csv $(DATADIR)/data_scale.csv $(DATADIR)/data_nrange.csv thp_compare | $(DATADIR)
	cd $(DATADIR) && ../thp_compare

run-proof: proof
	./proof

run-compare: compare_models
	./compare_models

run-ham: ham_cycle qpath dhd_hamconn
	./ham_cycle
	./qpath
	./dhd_hamconn

run-traffic: traffic_compare
	./traffic_compare

run-throughput: throughput_compare
	./throughput_compare

clean:
	rm -f $(BINS) *.o
	rm -rf $(DATADIR)

.PHONY: all csv clean run-proof run-compare run-ham run-traffic run-throughput
