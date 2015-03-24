CXXFLAGS := -O2 -Wall -Wextra

EXECUTABLES := bin/cutoffpassscore bin/Scores4mStruct bin/scoresequence

all: $(EXECUTABLES)

bin/scoresequence: src/scoresequence.o
	$(CXX) $(CXXFLAGS) $< -o $@	

bin/cutoffpassscore: src/cutoffpassscore.o
	$(CXX) $(CXXFLAGS) $< -o $@	

bin/Scores4mStruct: src/Scores4mStruct.o
	$(CXX) $(CXXFLAGS) $< -o $@	



clean: 
	rm -f src/*.o
	rm -f $(EXECUTABLES)

.o: $(CXX) $(CXXFLAGS) -c $(<) -o $*.o
