CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra

test_gs:
	$(CXX) -o out gs_test.cpp
	./out ;
	rm ./out

test_lll:
	$(CXX) $(CXXFLAGS) -o test_lll lll_test.cpp
	./test_lll
	rm ./test_lll
