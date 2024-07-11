default: consensus_asm synteny edlib

consensus_asm: consensus_asm.cpp dot_graph.o graph_simplification.o
	g++ -g -o consensus_asm consensus_asm.cpp dot_graph.o graph_simplification.o -ledlib -lpthread -fopenmp -fconcepts
	cp consensus_asm ../

dot_graph.o: dot_graph.cpp dot_graph.hpp
	g++ -g -c dot_graph.cpp -fopenmp -fconcepts

graph_simplification.o: graph_simplification.cpp graph_simplification.hpp
	g++ -g -c graph_simplification.cpp -fconcepts

synteny: synteny_multi.cpp synteny_graph.o
	g++ -g -o synteny synteny_multi.cpp synteny_graph.o -lpthread -fopenmp -fconcepts
	cp synteny ../

synteny_graph.o: synteny_graph.cpp synteny_graph.hpp
	g++ -g -c synteny_graph.cpp -fopenmp -fconcepts

edlib: edlib.cpp
	g++ -g -o edlib edlib.cpp -ledlib
	cp edlib ../

clean:
	rm -f consensus_asm edlib dot_graph.o graph_simplification.o synteny synteny_graph.o