default: consensus_asm synteny edlib_align

consensus_asm: consensus_asm.cpp dot_graph.o graph_simplification.o
	g++ -g -o consensus_asm consensus_asm.cpp dot_graph.o graph_simplification.o -lpthread -fopenmp -fconcepts edlib/src/edlib.cpp -I edlib/include

dot_graph.o: dot_graph.cpp dot_graph.hpp
	g++ -g -c dot_graph.cpp -fopenmp -fconcepts edlib/src/edlib.cpp -I edlib/include

graph_simplification.o: graph_simplification.cpp graph_simplification.hpp
	g++ -g -c graph_simplification.cpp -fconcepts edlib/src/edlib.cpp -I edlib/include

synteny: synteny_multi.cpp synteny_graph.o
	g++ -g -o synteny synteny_multi.cpp synteny_graph.o -lpthread -fopenmp -fconcepts

synteny_graph.o: synteny_graph.cpp synteny_graph.hpp
	g++ -g -c synteny_graph.cpp -fopenmp -fconcepts

edlib_align: edlib.cpp 
	g++ -g -o edlib_align edlib.cpp edlib/src/edlib.cpp -I edlib/include

clean:
	rm -rf consensus_asm edlib_align dot_graph.o graph_simplification.o synteny synteny_graph.o