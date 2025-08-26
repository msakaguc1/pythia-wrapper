
find_library(PYTHIA8_LIBRARY
  NAMES pythia8
  HINTS ENV PYTHIA8_LIBRARY_DIR
  PATH_SUFFIXES lib
  REQUIRED
)

find_library(LHAPDF_LIBRARY
  NAMES LHAPDF
  HINTS ENV LHAPDF_LIBRARY_DIR
  PATH_SUFFIXES lib
  REQUIRED
)

message(STATUS "Pythia8 main lib:    ${PYTHIA8_LIBRARY}")
message(STATUS "LHAPDF6 main lib:    ${LHAPDF_LIBRARY}")

find_path(PYTHIA8_INCLUDE_DIR
  NAMES Pythia8/Pythia.h
  HINTS ENV PYTHIA8_INCLUDE_DIR
  PATH_SUFFIXES include
)

find_path(PYTHIA8_PLUGINS_INCLUDE_DIR
  NAMES LHAPDF6.h
  HINTS ENV PYTHIA8_INCLUDE_DIR
  PATH_SUFFIXES Pythia8Plugins
)

find_path(LHAPDF_INCLUDE_DIR
  NAMES LHAPDF/LHAPDF.h
  HINTS ENV LHAPDF_INCLUDE_DIR
  PATH_SUFFIXES ""
  REQUIRED
)

include_directories(SYSTEM
  ${PYTHIA8_INCLUDE_DIR}
  ${PYTHIA8_PLUGINS_INCLUDE_DIR}
  ${LHAPDF_INCLUDE_DIR}
)

i3_project(pythia-wrapper
        DOCS_DIR resources/docs
        PYTHON_DIR python
        PYTHON_DEST icecube/    
)

SET(LIB_${PROJECT_NAME}_SOURCEFILES
        private/pythia-wrapper/I3PythiaWrapperModule.cxx
)

SET(LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES
       private/pybindings/module.cxx
)

SET(LIB_${PROJECT_NAME}_PROJECTS
       serialization icetray dataclasses phys-services simclasses sim-services
)

SET(LIB_${PROJECT_NAME}_TOOLS python boost)

i3_add_library(${PROJECT_NAME}
        ${LIB_${PROJECT_NAME}_SOURCEFILES}
        USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
        USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

#needed to include libpythia8.so in libpythia-wrapper.so
target_link_libraries(${PROJECT_NAME}
    ${PYTHIA8_LIBRARY}
    ${LHAPDF_LIBRARY}
    icetray
    dataclasses
    boost_python
)

i3_add_pybindings(${PROJECT_NAME}
  ${LIB_${PROJECT_NAME}_PYBINDINGS_SOURCEFILES}
  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
  USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)


# else(PYTHIA_FOUND)
#         colormsg(CYAN "+-- pythia not found, skipping pythia-wrapper")
