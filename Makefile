CXXFLAGS := -O2 -Wall -Wextra

EXECUTABLES := bin/cutoffpassscore bin/scorestructure bin/scoresequence

all: $(EXECUTABLES)

bin/scoresequence: src/scoresequence.o
	$(CXX) $(CXXFLAGS) $< -o $@	

bin/cutoffpassscore: src/cutoffpassscore.o
	$(CXX) $(CXXFLAGS) $< -o $@	

bin/scorestructure: src/scorestructure.o src/mirnastats.h src/mirnastats.o
	$(CXX) $(CXXFLAGS) $< src/mirnastats.o -o $@	




clean: 
	rm -f src/*.o
	rm -f $(EXECUTABLES)

.o: $(CXX) $(CXXFLAGS) -c $(<) -o $*.o
