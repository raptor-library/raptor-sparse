# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bienz/Documents/Programs/raptor_master/raptor-sparse

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build

# Include any dependencies generated for this target.
include src/gallery/tests/CMakeFiles/test_matrix_market.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/gallery/tests/CMakeFiles/test_matrix_market.dir/compiler_depend.make

# Include the progress variables for this target.
include src/gallery/tests/CMakeFiles/test_matrix_market.dir/progress.make

# Include the compile flags for this target's objects.
include src/gallery/tests/CMakeFiles/test_matrix_market.dir/flags.make

src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o: src/gallery/tests/CMakeFiles/test_matrix_market.dir/flags.make
src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o: /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/gallery/tests/test_matrix_market.cpp
src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o: src/gallery/tests/CMakeFiles/test_matrix_market.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o -MF CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o.d -o CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o -c /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/gallery/tests/test_matrix_market.cpp

src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.i"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/gallery/tests/test_matrix_market.cpp > CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.i

src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.s"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/gallery/tests/test_matrix_market.cpp -o CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.s

# Object files for target test_matrix_market
test_matrix_market_OBJECTS = \
"CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o"

# External object files for target test_matrix_market
test_matrix_market_EXTERNAL_OBJECTS =

src/gallery/tests/test_matrix_market: src/gallery/tests/CMakeFiles/test_matrix_market.dir/test_matrix_market.cpp.o
src/gallery/tests/test_matrix_market: src/gallery/tests/CMakeFiles/test_matrix_market.dir/build.make
src/gallery/tests/test_matrix_market: lib/libraptor-sparse.dylib
src/gallery/tests/test_matrix_market: /opt/homebrew/Cellar/mpich/4.1.2/lib/libmpicxx.dylib
src/gallery/tests/test_matrix_market: /opt/homebrew/Cellar/mpich/4.1.2/lib/libmpi.dylib
src/gallery/tests/test_matrix_market: /opt/homebrew/Cellar/mpich/4.1.2/lib/libpmpi.dylib
src/gallery/tests/test_matrix_market: /Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/liblapack.tbd
src/gallery/tests/test_matrix_market: /Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/libblas.tbd
src/gallery/tests/test_matrix_market: lib/libgtest_main.1.12.1.dylib
src/gallery/tests/test_matrix_market: lib/libgtest.1.12.1.dylib
src/gallery/tests/test_matrix_market: src/gallery/tests/CMakeFiles/test_matrix_market.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_matrix_market"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_matrix_market.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/gallery/tests/CMakeFiles/test_matrix_market.dir/build: src/gallery/tests/test_matrix_market
.PHONY : src/gallery/tests/CMakeFiles/test_matrix_market.dir/build

src/gallery/tests/CMakeFiles/test_matrix_market.dir/clean:
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_matrix_market.dir/cmake_clean.cmake
.PHONY : src/gallery/tests/CMakeFiles/test_matrix_market.dir/clean

src/gallery/tests/CMakeFiles/test_matrix_market.dir/depend:
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bienz/Documents/Programs/raptor_master/raptor-sparse /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/gallery/tests /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/gallery/tests/CMakeFiles/test_matrix_market.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/gallery/tests/CMakeFiles/test_matrix_market.dir/depend

