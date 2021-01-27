# - Config file for the adaptnotch package
# It defines the following variables
#  ADAPTNOTCH_LIBRARIES    - libraries to link against

# Compute paths
get_filename_component(ADAPTNOTCH_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET adaptnotch)
  include("${ADAPTNOTCH_CMAKE_DIR}/adaptnotch-targets.cmake")
endif()

# These are IMPORTED targets created by adaptnotch-targets.cmake
set(ADAPTNOTCH_LIBRARIES adaptnotch)
