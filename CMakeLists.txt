cmake_minimum_required(VERSION 3.10)

# set the project name
project(ScatteringRadiation)

# Include Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED)
#find_package(geant4 REQUIRED)

# Set it up
include(${Geant4_USE_FILE})

# Include my headers dir
include_directories(include
	${Geant4_INCLUDE_DIR}
	${ROOT_INCLUDE_DIRS})

# Make me a setup script directly in the build directory
configure_file(tools/setup.sh.in setup.sh @ONLY)

# add the executable
add_executable(run_beam src/run_beam.cxx
	src/Construction.cxx
	src/Beam.cxx
	src/Action.cxx
	src/SensitiveSurface.cxx
	src/Run.cxx
	src/Event.cxx
	src/Stepping.cxx
	src/Track.cxx)
	
target_link_libraries(run_beam ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

