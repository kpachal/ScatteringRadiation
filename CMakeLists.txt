cmake_minimum_required(VERSION 3.10)

# set the project name
project(ScatteringRadiation)

# Make me a setup script directly in the build directory
configure_file(tools/setup.sh.in setup.sh @ONLY)

# add the executable
add_executable(run_beam src/run_beam.cxx)
