ADD_SUBDIRECTORY(dcmip1-1)
ADD_SUBDIRECTORY(dcmip1-2)

file(GLOB NCL_SCRIPTS    "*.ncl")                                       # get ncl-scripts
file(GLOB SHELL_SCRIPTS  "*.sh")                                        # get shell-scripts
file(GLOB PYTHON_SCRIPTS "*.py")                                        # get python scripts

# Copy test files to build dir
install(PROGRAMS ${NCL_SCRIPTS} ${SHELL_SCRIPTS} ${NAMELISTS} ${PYTHON_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
