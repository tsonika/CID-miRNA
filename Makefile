CXXFLAGS := -O2 -Wall -Wextra

EXECUTABLES := bin/cutoffpassscore bin/Scores4mStruct

all: $(EXECUTABLES)

parser4auto: parser4auto.o
	$(CXX) $(CXXFLAGS) $< -o $@

#newcyk2: newcyk2.o
#	$(CXX) $(CXXFLAGS) $< -o $@	

bin/cutoffpassscore: src/cutoffpassscore.o
	$(CXX) $(CXXFLAGS) $< -o $@	

bin/Scores4mStruct: src/Scores4mStruct.o
	$(CXX) $(CXXFLAGS) $< -o $@	



clean: 
	rm -f src/*.o
	rm -f $(EXECUTABLES)

.o: $(CXX) $(CXXFLAGS) -c $(<) -o $*.o
