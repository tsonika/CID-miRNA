CXXFLAGS := -O2 -Wall -Wextra -mtune=native -march=native

EXECUTABLES := parser4auto cutoffpassscore Scores4mStruct

all: $(EXECUTABLES)

parser4auto: parser4auto.o
	$(CXX) $(CXXFLAGS) $< -o $@

#newcyk2: newcyk2.o
#	$(CXX) $(CXXFLAGS) $< -o $@	

cutoffpassscore: cutoffpassscore.o
	$(CXX) $(CXXFLAGS) $< -o $@	

Scores4mStruct: Scores4mStruct.o
	$(CXX) $(CXXFLAGS) $< -o $@	



clean: 
	rm -f *.o
	rm -f $(EXECUTABLES)

.o: $(CXX) $(CXXFLAGS) -c $(<) -o $*.o
