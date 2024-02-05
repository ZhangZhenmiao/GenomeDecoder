default: consensus_asm edlib

consensus_asm: consensus_asm.cpp dot_graph.o graph_simplification.o
	g++ -g -o consensus_asm consensus_asm.cpp dot_graph.o graph_simplification.o -ledlib -lpthread -fopenmp -fconcepts
	cp consensus_asm ../

dot_graph.o: dot_graph.cpp dot_graph.hpp
	g++ -g -c dot_graph.cpp -fopenmp -fconcepts

graph_simplification.o: graph_simplification.cpp graph_simplification.hpp
	g++ -g -c graph_simplification.cpp -fconcepts

edlib: edlib.cpp
	g++ -g -o edlib edlib.cpp -ledlib
	cp edlib ../

clean:
	rm -f consensus_asm edlib dot_graph.o graph_simplification.o