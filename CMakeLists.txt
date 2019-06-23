cmake_minimum_required(VERSION 3.7)

project(hasenpfote CXX)

# Check if supported platforms and compilers.
if(WIN32)
    if(MSVC)
        if(MSVC_VERSION GREATER_EQUAL 1900)
        else()
            message(FATAL_ERROR "You are using an unsupported MSVC version!")
        endif()
    elseif()
        message(FATAL_ERROR "You are using an unsupported Windows compiler!")
    endif()
elseif(UNIX)
    option(FORCE_32BIT_BUILD "Force a 32bit compile on 64bit" OFF)
else()
    message(FATAL_ERROR "You are on an unsupported platform!")
endif()

# Additional cmake files.
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Compiler settings.
include(compiler_settings)

# The relative path from `CMAKE_CURRENT_SOURCE_DIR`.
set(SOURCES_DIR "sources")

# Fetch files.
file(GLOB_RECURSE HEADER_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" ${SOURCES_DIR}/*.h)
file(GLOB_RECURSE SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" ${SOURCES_DIR}/*.cpp)
set(ALL_FILES ${HEADER_FILES} ${SOURCE_FILES})

# Define a grouping for source files in IDE project generation.
if(WIN32 AND MSVC)
    include(utils)
    do_grouping_for_files("Header Files" ${SOURCES_DIR} "${HEADER_FILES}")
    do_grouping_for_files("Source Files" ${SOURCES_DIR} "${SOURCE_FILES}")
endif()

# Targets.
add_library(${PROJECT_NAME} STATIC ${ALL_FILES})
if(FORCE_32BIT_BUILD)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

################################################################################

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Copy header files for build tree.
foreach(HEADER_FILE IN LISTS HEADER_FILES)
	get_filename_component(HEADER_FILE_DIR ${HEADER_FILE} DIRECTORY)
	get_filename_component(HEADER_FILE_NAME ${HEADER_FILE} NAME_WE)
	string(REGEX REPLACE "^${SOURCES_DIR}" "${PROJECT_NAME}" DST_DIR ${HEADER_FILE_DIR})
	set(DST_DIR ${PROJECT_BINARY_DIR}/include/${DST_DIR})
	string(REPLACE "/" "_" TARGET_NAME copy_header_${HEADER_FILE_DIR}_${HEADER_FILE_NAME})
	add_custom_target(${TARGET_NAME}
		COMMAND ${CMAKE_COMMAND} -E make_directory ${DST_DIR}
		COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_FILE} ${DST_DIR}
		)
	add_dependencies(${PROJECT_NAME} ${TARGET_NAME})
endforeach()

# Install the header files.
# Instalation directory: `${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}`
install(DIRECTORY ${SOURCES_DIR}/ DESTINATION include/${PROJECT_NAME} FILES_MATCHING PATTERN "*.h")

# Install the archive files.
# Instalation directory: `${CMAKE_INSTALL_PREFIX}/share/cmake/${PROJECT_NAME}`
string(TOLOWER ${PROJECT_NAME} LOWER_PROJECT_NAME)
set(EXPORT_NAME ${LOWER_PROJECT_NAME}-export)
set(CONFIG_NAME ${LOWER_PROJECT_NAME}-config)

# For Debug.
install(TARGETS ${PROJECT_NAME} EXPORT ${EXPORT_NAME}-Debug ARCHIVE DESTINATION lib/Debug CONFIGURATIONS Debug)
install(EXPORT ${EXPORT_NAME}-Debug DESTINATION share/cmake/${PROJECT_NAME} FILE ${CONFIG_NAME}.cmake CONFIGURATIONS Debug EXPORT_LINK_INTERFACE_LIBRARIES)
export(EXPORT ${EXPORT_NAME}-Debug FILE ${EXPORT_NAME}.cmake)
# For Release.
install(TARGETS ${PROJECT_NAME} EXPORT ${EXPORT_NAME}-Release ARCHIVE DESTINATION lib/Release CONFIGURATIONS Release)
install(EXPORT ${EXPORT_NAME}-Release DESTINATION share/cmake/${PROJECT_NAME} FILE ${CONFIG_NAME}.cmake CONFIGURATIONS Release EXPORT_LINK_INTERFACE_LIBRARIES)
export(EXPORT ${EXPORT_NAME}-Release FILE ${EXPORT_NAME}.cmake)
# For MinSizeRel.
install(TARGETS ${PROJECT_NAME} EXPORT ${EXPORT_NAME}-MinSizeRel ARCHIVE DESTINATION lib/MinSizeRel CONFIGURATIONS MinSizeRel)
install(EXPORT ${EXPORT_NAME}-MinSizeRel DESTINATION share/cmake/${PROJECT_NAME} FILE ${CONFIG_NAME}.cmake CONFIGURATIONS MinSizeRel EXPORT_LINK_INTERFACE_LIBRARIES)
export(EXPORT ${EXPORT_NAME}-MinSizeRel FILE ${EXPORT_NAME}.cmake)
# For RelWithDebInfo.
install(TARGETS ${PROJECT_NAME} EXPORT ${EXPORT_NAME}-RelWithDebInfo ARCHIVE DESTINATION lib/RelWithDebInfo CONFIGURATIONS RelWithDebInfo)
install(EXPORT ${EXPORT_NAME}-RelWithDebInfo DESTINATION share/cmake/${PROJECT_NAME} FILE ${CONFIG_NAME}.cmake CONFIGURATIONS RelWithDebInfo EXPORT_LINK_INTERFACE_LIBRARIES)
export(EXPORT ${EXPORT_NAME}-RelWithDebInfo FILE ${EXPORT_NAME}.cmake)
