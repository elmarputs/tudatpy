# Configure the version file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_BINARY_DIR}/_version.py" @ONLY)

# Configure the files needed to make the python wheels (for PyPi packages)
#if(MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    message(STATUS "Creating the files for the generation of a binary wheel.")
#    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/../wheel/setup.py" @ONLY)
#    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#        # NOTE: this is necessary on linux but harmful on mingw.
#        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/../wheel/setup.cfg" @ONLY)
#    endif()
#    if(MINGW)
#        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
#    endif()
#endif()

if (MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Creating the files for the installation from build/tudatpy.")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/../setup.py" @ONLY)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        #        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/../setup.cfg" @ONLY)
    endif ()
    if (MINGW)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
    endif ()
endif ()


# The list of pygmo's Python files.
set(TUDATPY_PYTHON_FILES
        __init__.py
        test.py
        _layer_simulation_setup.py
        )


# Copy the python files in the current binary dir,
# so that we can import pygmo from the build dir.
# NOTE: importing from the build dir will work
# only on single-configuration generators.
foreach (TUDATPY_PYTHON_FILE ${TUDATPY_PYTHON_FILES})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${TUDATPY_PYTHON_FILE}"
            "${CMAKE_CURRENT_BINARY_DIR}/${TUDATPY_PYTHON_FILE}" COPYONLY)
endforeach ()

include_directories(${Tudat_INCLUDE_DIRS})
# Core module.
YACMA_PYTHON_MODULE(core
        kernel/core.cpp
        kernel/expose_constants.cpp
        kernel/expose_spice_interface.cpp
        kernel/expose_interpolators.cpp
        kernel/expose_basic_astrodynamics.cpp
        kernel/expose_gravitation.cpp
        kernel/expose_orbital_element_conversions.cpp
        kernel/expose_simulation_setup.cpp
        kernel/expose_numerical_integrators.cpp
        kernel/expose_propagators.cpp
        #        expose_mathematics.cpp
        #        trampoline_classes.cpp
        kernel/docstrings.cpp
        #        expose_ephemeris.cpp
        )

target_link_libraries(core PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${PYTHON_LIBRARIES}
        #        ${tudat_basic_mathematics}
        ${Tudat_PROPAGATION_LIBRARIES}
        #        ${tudat_environment_setup}
        #        ${TUDAT_ESTIMATION_LIBRARIES}
        #        ${tudat_basic_mathematics}
        )

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${pybind11_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")
set_target_properties(core PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(core PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")

# Install the core module.
install(TARGETS core
        RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
        LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
        )

# Install the Python files.
install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
        DESTINATION ${TUDATPY_INSTALL_PATH})

#file(COPY __init__.py DESTINATION .)