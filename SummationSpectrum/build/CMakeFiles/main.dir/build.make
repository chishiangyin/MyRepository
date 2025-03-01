# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /junofs/users/yinqixiang/package/SummationSpectrum

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /junofs/users/yinqixiang/package/SummationSpectrum/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/main.cc
CMakeFiles/main.dir/main.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/main.cc.o -MF CMakeFiles/main.dir/main.cc.o.d -o CMakeFiles/main.dir/main.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/main.cc

CMakeFiles/main.dir/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/main.cc > CMakeFiles/main.dir/main.cc.i

CMakeFiles/main.dir/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/main.cc -o CMakeFiles/main.dir/main.cc.s

CMakeFiles/main.dir/src/BetaDecayBranch.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/BetaDecayBranch.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayBranch.cc
CMakeFiles/main.dir/src/BetaDecayBranch.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/src/BetaDecayBranch.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/BetaDecayBranch.cc.o -MF CMakeFiles/main.dir/src/BetaDecayBranch.cc.o.d -o CMakeFiles/main.dir/src/BetaDecayBranch.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayBranch.cc

CMakeFiles/main.dir/src/BetaDecayBranch.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/BetaDecayBranch.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayBranch.cc > CMakeFiles/main.dir/src/BetaDecayBranch.cc.i

CMakeFiles/main.dir/src/BetaDecayBranch.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/BetaDecayBranch.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayBranch.cc -o CMakeFiles/main.dir/src/BetaDecayBranch.cc.s

CMakeFiles/main.dir/src/BetaDecayInput.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/BetaDecayInput.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayInput.cc
CMakeFiles/main.dir/src/BetaDecayInput.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/src/BetaDecayInput.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/BetaDecayInput.cc.o -MF CMakeFiles/main.dir/src/BetaDecayInput.cc.o.d -o CMakeFiles/main.dir/src/BetaDecayInput.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayInput.cc

CMakeFiles/main.dir/src/BetaDecayInput.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/BetaDecayInput.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayInput.cc > CMakeFiles/main.dir/src/BetaDecayInput.cc.i

CMakeFiles/main.dir/src/BetaDecayInput.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/BetaDecayInput.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/BetaDecayInput.cc -o CMakeFiles/main.dir/src/BetaDecayInput.cc.s

CMakeFiles/main.dir/src/FissionInput.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/FissionInput.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/FissionInput.cc
CMakeFiles/main.dir/src/FissionInput.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/main.dir/src/FissionInput.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/FissionInput.cc.o -MF CMakeFiles/main.dir/src/FissionInput.cc.o.d -o CMakeFiles/main.dir/src/FissionInput.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/FissionInput.cc

CMakeFiles/main.dir/src/FissionInput.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/FissionInput.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/FissionInput.cc > CMakeFiles/main.dir/src/FissionInput.cc.i

CMakeFiles/main.dir/src/FissionInput.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/FissionInput.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/FissionInput.cc -o CMakeFiles/main.dir/src/FissionInput.cc.s

CMakeFiles/main.dir/src/InitializeInput.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/InitializeInput.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/InitializeInput.cc
CMakeFiles/main.dir/src/InitializeInput.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/main.dir/src/InitializeInput.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/InitializeInput.cc.o -MF CMakeFiles/main.dir/src/InitializeInput.cc.o.d -o CMakeFiles/main.dir/src/InitializeInput.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/InitializeInput.cc

CMakeFiles/main.dir/src/InitializeInput.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/InitializeInput.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/InitializeInput.cc > CMakeFiles/main.dir/src/InitializeInput.cc.i

CMakeFiles/main.dir/src/InitializeInput.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/InitializeInput.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/InitializeInput.cc -o CMakeFiles/main.dir/src/InitializeInput.cc.s

CMakeFiles/main.dir/src/Input.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/Input.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/Input.cc
CMakeFiles/main.dir/src/Input.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/main.dir/src/Input.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/Input.cc.o -MF CMakeFiles/main.dir/src/Input.cc.o.d -o CMakeFiles/main.dir/src/Input.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/Input.cc

CMakeFiles/main.dir/src/Input.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/Input.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/Input.cc > CMakeFiles/main.dir/src/Input.cc.i

CMakeFiles/main.dir/src/Input.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/Input.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/Input.cc -o CMakeFiles/main.dir/src/Input.cc.s

CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o: /junofs/users/yinqixiang/package/SummationSpectrum/src/IsotopeSpectrum.cc
CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o -MF CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o.d -o CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o -c /junofs/users/yinqixiang/package/SummationSpectrum/src/IsotopeSpectrum.cc

CMakeFiles/main.dir/src/IsotopeSpectrum.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/src/IsotopeSpectrum.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /junofs/users/yinqixiang/package/SummationSpectrum/src/IsotopeSpectrum.cc > CMakeFiles/main.dir/src/IsotopeSpectrum.cc.i

CMakeFiles/main.dir/src/IsotopeSpectrum.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/src/IsotopeSpectrum.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /junofs/users/yinqixiang/package/SummationSpectrum/src/IsotopeSpectrum.cc -o CMakeFiles/main.dir/src/IsotopeSpectrum.cc.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cc.o" \
"CMakeFiles/main.dir/src/BetaDecayBranch.cc.o" \
"CMakeFiles/main.dir/src/BetaDecayInput.cc.o" \
"CMakeFiles/main.dir/src/FissionInput.cc.o" \
"CMakeFiles/main.dir/src/InitializeInput.cc.o" \
"CMakeFiles/main.dir/src/Input.cc.o" \
"CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

main: CMakeFiles/main.dir/main.cc.o
main: CMakeFiles/main.dir/src/BetaDecayBranch.cc.o
main: CMakeFiles/main.dir/src/BetaDecayInput.cc.o
main: CMakeFiles/main.dir/src/FissionInput.cc.o
main: CMakeFiles/main.dir/src/InitializeInput.cc.o
main: CMakeFiles/main.dir/src/Input.cc.o
main: CMakeFiles/main.dir/src/IsotopeSpectrum.cc.o
main: CMakeFiles/main.dir/build.make
main: libmylib.so
main: /usr/lib64/root/libCore.so
main: /usr/lib64/root/libImt.so
main: /usr/lib64/root/libRIO.so
main: /usr/lib64/root/libNet.so
main: /usr/lib64/root/libHist.so
main: /usr/lib64/root/libGraf.so
main: /usr/lib64/root/libGraf3d.so
main: /usr/lib64/root/libGpad.so
main: /usr/lib64/root/libROOTDataFrame.so
main: /usr/lib64/root/libTree.so
main: /usr/lib64/root/libTreePlayer.so
main: /usr/lib64/root/libRint.so
main: /usr/lib64/root/libPostscript.so
main: /usr/lib64/root/libMatrix.so
main: /usr/lib64/root/libPhysics.so
main: /usr/lib64/root/libMathCore.so
main: /usr/lib64/root/libThread.so
main: /usr/lib64/root/libMultiProc.so
main: /usr/lib64/root/libROOTVecOps.so
main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /junofs/users/yinqixiang/package/SummationSpectrum/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /junofs/users/yinqixiang/package/SummationSpectrum /junofs/users/yinqixiang/package/SummationSpectrum /junofs/users/yinqixiang/package/SummationSpectrum/build /junofs/users/yinqixiang/package/SummationSpectrum/build /junofs/users/yinqixiang/package/SummationSpectrum/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

