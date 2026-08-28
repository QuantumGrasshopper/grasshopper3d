
# declare variables
CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -O3 -Wall -Wextra -std=c++17
LDFLAGS ?=
LDLIBS ?= -lgsl -lgslcblas -lm
PYTHON ?= python3

OBJS = main.o interactions.o parameters.o setup.o annealing.o output.o
TEST_BIN = tests/grasshopper_tests
TEST_SRCS = tests/test_main.cpp \
		tests/test_globals.cpp \
		tests/test_interactions.cpp \
		tests/test_parameters.cpp \
		tests/test_setup.cpp \
		tests/test_annealing.cpp \
		tests/test_output.cpp
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
TEST_PRODUCTION_OBJS = interactions.o parameters.o setup.o annealing.o output.o

# link .o-files to program
grasshopper:  $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o grasshopper $(LDLIBS)

$(TEST_BIN): $(TEST_OBJS) $(TEST_PRODUCTION_OBJS)
	$(CXX) $(LDFLAGS) $(TEST_OBJS) $(TEST_PRODUCTION_OBJS) -o $(TEST_BIN) $(LDLIBS)

# create .o-files from .cpp-files
main.o: main.cpp utilities.hpp setup.hpp annealing.hpp interactions.hpp parameters.hpp output.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c main.cpp
interactions.o: interactions.cpp interactions.hpp utilities.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c interactions.cpp
parameters.o: parameters.cpp parameters.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c parameters.cpp
setup.o: setup.cpp setup.hpp utilities.hpp output.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c setup.cpp
annealing.o: annealing.cpp annealing.hpp utilities.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c annealing.cpp
output.o: output.cpp output.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c output.cpp

tests/%.o: tests/%.cpp tests/doctest/doctest.h utilities.hpp interactions.hpp parameters.hpp setup.hpp annealing.hpp output.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I. -Itests -c $< -o $@

test: $(TEST_BIN)
	./$(TEST_BIN)

integration-test: grasshopper
	$(PYTHON) tests/run_integration_tests.py ./grasshopper

check: test integration-test

# clean up
.PHONY: clean tidy test integration-test check
clean:
	rm -f *~ *.o $(TEST_OBJS) $(TEST_BIN)
tidy:	clean
	rm -f grasshopper
