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
include src/core/tests/CMakeFiles/test_tap_comm.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/core/tests/CMakeFiles/test_tap_comm.dir/compiler_depend.make

# Include the progress variables for this target.
include src/core/tests/CMakeFiles/test_tap_comm.dir/progress.make

# Include the compile flags for this target's objects.
include src/core/tests/CMakeFiles/test_tap_comm.dir/flags.make

src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o: src/core/tests/CMakeFiles/test_tap_comm.dir/flags.make
src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o: /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/core/tests/test_tap_comm.cpp
src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o: src/core/tests/CMakeFiles/test_tap_comm.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o -MF CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o.d -o CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o -c /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/core/tests/test_tap_comm.cpp

src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.i"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/core/tests/test_tap_comm.cpp > CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.i

src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.s"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/core/tests/test_tap_comm.cpp -o CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.s

# Object files for target test_tap_comm
test_tap_comm_OBJECTS = \
"CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o"

# External object files for target test_tap_comm
test_tap_comm_EXTERNAL_OBJECTS =

src/core/tests/test_tap_comm: src/core/tests/CMakeFiles/test_tap_comm.dir/test_tap_comm.cpp.o
src/core/tests/test_tap_comm: src/core/tests/CMakeFiles/test_tap_comm.dir/build.make
src/core/tests/test_tap_comm: lib/libraptor-sparse.dylib
src/core/tests/test_tap_comm: /opt/homebrew/Cellar/mpich/4.1.2/lib/libmpicxx.dylib
src/core/tests/test_tap_comm: /opt/homebrew/Cellar/mpich/4.1.2/lib/libmpi.dylib
src/core/tests/test_tap_comm: /opt/homebrew/Cellar/mpich/4.1.2/lib/libpmpi.dylib
src/core/tests/test_tap_comm: /Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/liblapack.tbd
src/core/tests/test_tap_comm: /Library/Developer/CommandLineTools/SDKs/MacOSX13.3.sdk/usr/lib/libblas.tbd
src/core/tests/test_tap_comm: lib/libgtest_main.1.12.1.dylib
src/core/tests/test_tap_comm: lib/libgtest.1.12.1.dylib
src/core/tests/test_tap_comm: src/core/tests/CMakeFiles/test_tap_comm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_tap_comm"
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_tap_comm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/core/tests/CMakeFiles/test_tap_comm.dir/build: src/core/tests/test_tap_comm
.PHONY : src/core/tests/CMakeFiles/test_tap_comm.dir/build

src/core/tests/CMakeFiles/test_tap_comm.dir/clean:
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests && $(CMAKE_COMMAND) -P CMakeFiles/test_tap_comm.dir/cmake_clean.cmake
.PHONY : src/core/tests/CMakeFiles/test_tap_comm.dir/clean

src/core/tests/CMakeFiles/test_tap_comm.dir/depend:
	cd /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bienz/Documents/Programs/raptor_master/raptor-sparse /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/src/core/tests /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests /Users/bienz/Documents/Programs/raptor_master/raptor-sparse/build/src/core/tests/CMakeFiles/test_tap_comm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/tests/CMakeFiles/test_tap_comm.dir/depend

