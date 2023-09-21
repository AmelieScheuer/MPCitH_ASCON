# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/crypto/Documents/implementation_cpp/keccak

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/crypto/Documents/implementation_cpp/keccak/build

# Include any dependencies generated for this target.
include CMakeFiles/keccak.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/keccak.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/keccak.dir/flags.make

CMakeFiles/keccak.dir/KeccakHash.c.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/KeccakHash.c.o: ../KeccakHash.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/keccak.dir/KeccakHash.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/keccak.dir/KeccakHash.c.o   -c /home/crypto/Documents/implementation_cpp/keccak/KeccakHash.c

CMakeFiles/keccak.dir/KeccakHash.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/keccak.dir/KeccakHash.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/crypto/Documents/implementation_cpp/keccak/KeccakHash.c > CMakeFiles/keccak.dir/KeccakHash.c.i

CMakeFiles/keccak.dir/KeccakHash.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/keccak.dir/KeccakHash.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/crypto/Documents/implementation_cpp/keccak/KeccakHash.c -o CMakeFiles/keccak.dir/KeccakHash.c.s

CMakeFiles/keccak.dir/KeccakSponge.c.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/KeccakSponge.c.o: ../KeccakSponge.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/keccak.dir/KeccakSponge.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/keccak.dir/KeccakSponge.c.o   -c /home/crypto/Documents/implementation_cpp/keccak/KeccakSponge.c

CMakeFiles/keccak.dir/KeccakSponge.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/keccak.dir/KeccakSponge.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/crypto/Documents/implementation_cpp/keccak/KeccakSponge.c > CMakeFiles/keccak.dir/KeccakSponge.c.i

CMakeFiles/keccak.dir/KeccakSponge.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/keccak.dir/KeccakSponge.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/crypto/Documents/implementation_cpp/keccak/KeccakSponge.c -o CMakeFiles/keccak.dir/KeccakSponge.c.s

CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o: ../KeccakSpongetimes4.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o   -c /home/crypto/Documents/implementation_cpp/keccak/KeccakSpongetimes4.c

CMakeFiles/keccak.dir/KeccakSpongetimes4.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/keccak.dir/KeccakSpongetimes4.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/crypto/Documents/implementation_cpp/keccak/KeccakSpongetimes4.c > CMakeFiles/keccak.dir/KeccakSpongetimes4.c.i

CMakeFiles/keccak.dir/KeccakSpongetimes4.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/keccak.dir/KeccakSpongetimes4.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/crypto/Documents/implementation_cpp/keccak/KeccakSpongetimes4.c -o CMakeFiles/keccak.dir/KeccakSpongetimes4.c.s

CMakeFiles/keccak.dir/KeccakHashtimes4.c.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/KeccakHashtimes4.c.o: ../KeccakHashtimes4.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object CMakeFiles/keccak.dir/KeccakHashtimes4.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/keccak.dir/KeccakHashtimes4.c.o   -c /home/crypto/Documents/implementation_cpp/keccak/KeccakHashtimes4.c

CMakeFiles/keccak.dir/KeccakHashtimes4.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/keccak.dir/KeccakHashtimes4.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/crypto/Documents/implementation_cpp/keccak/KeccakHashtimes4.c > CMakeFiles/keccak.dir/KeccakHashtimes4.c.i

CMakeFiles/keccak.dir/KeccakHashtimes4.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/keccak.dir/KeccakHashtimes4.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/crypto/Documents/implementation_cpp/keccak/KeccakHashtimes4.c -o CMakeFiles/keccak.dir/KeccakHashtimes4.c.s

CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o: ../avx2/KeccakP-1600-AVX2.s
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building ASM object CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o"
	/usr/bin/cc $(ASM_DEFINES) $(ASM_INCLUDES) $(ASM_FLAGS) -x assembler-with-cpp -o CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o -c /home/crypto/Documents/implementation_cpp/keccak/avx2/KeccakP-1600-AVX2.s

CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o: CMakeFiles/keccak.dir/flags.make
CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o: ../avx2/KeccakP-1600-times4-SIMD256.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o   -c /home/crypto/Documents/implementation_cpp/keccak/avx2/KeccakP-1600-times4-SIMD256.c

CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/crypto/Documents/implementation_cpp/keccak/avx2/KeccakP-1600-times4-SIMD256.c > CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.i

CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/crypto/Documents/implementation_cpp/keccak/avx2/KeccakP-1600-times4-SIMD256.c -o CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.s

# Object files for target keccak
keccak_OBJECTS = \
"CMakeFiles/keccak.dir/KeccakHash.c.o" \
"CMakeFiles/keccak.dir/KeccakSponge.c.o" \
"CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o" \
"CMakeFiles/keccak.dir/KeccakHashtimes4.c.o" \
"CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o" \
"CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o"

# External object files for target keccak
keccak_EXTERNAL_OBJECTS =

libkeccak.a: CMakeFiles/keccak.dir/KeccakHash.c.o
libkeccak.a: CMakeFiles/keccak.dir/KeccakSponge.c.o
libkeccak.a: CMakeFiles/keccak.dir/KeccakSpongetimes4.c.o
libkeccak.a: CMakeFiles/keccak.dir/KeccakHashtimes4.c.o
libkeccak.a: CMakeFiles/keccak.dir/avx2/KeccakP-1600-AVX2.s.o
libkeccak.a: CMakeFiles/keccak.dir/avx2/KeccakP-1600-times4-SIMD256.c.o
libkeccak.a: CMakeFiles/keccak.dir/build.make
libkeccak.a: CMakeFiles/keccak.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking C static library libkeccak.a"
	$(CMAKE_COMMAND) -P CMakeFiles/keccak.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/keccak.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/keccak.dir/build: libkeccak.a

.PHONY : CMakeFiles/keccak.dir/build

CMakeFiles/keccak.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/keccak.dir/cmake_clean.cmake
.PHONY : CMakeFiles/keccak.dir/clean

CMakeFiles/keccak.dir/depend:
	cd /home/crypto/Documents/implementation_cpp/keccak/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/crypto/Documents/implementation_cpp/keccak /home/crypto/Documents/implementation_cpp/keccak /home/crypto/Documents/implementation_cpp/keccak/build /home/crypto/Documents/implementation_cpp/keccak/build /home/crypto/Documents/implementation_cpp/keccak/build/CMakeFiles/keccak.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/keccak.dir/depend
