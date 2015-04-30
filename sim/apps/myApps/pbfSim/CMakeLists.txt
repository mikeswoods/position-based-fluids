cmake_minimum_required(VERSION 3.2)
project(pbfSim)

if(CMAKE_COMPILER_IS_GNUCXX)
	add_definitions(-std=gnu++0x)
endif()

set(SYSTEM_INCLUDE "/usr/include")
set(OFX_ADDONS "../../../addons")
set(OFX_LIBS "../../../libs")
set(OFX_MSAOPENCL "${OFX_ADDONS}/ofxMSAOpenCL/src")
set(OFX_GUI "${OFX_ADDONS}/ofxGui/src")
set(DEP_CAIRO "${OFX_LIBS}/cairo")
set(DEP_KISS "${OFX_LIBS}/kiss")
set(DEP_OPENCL "${OFX_LIBS}/opencl")
set(DEP_POCO "${OFX_LIBS}/poco")
set(DEP_TESS2 "${OFX_LIBS}/tess2")

include_directories("${SYSTEM_INCLUDE}"
	                "${OFX_LIBS}/openFrameworks"
	                "${OFX_LIBS}/openFrameworks/3d"
	                "${OFX_LIBS}/openFrameworks/app"
	                "${OFX_LIBS}/openFrameworks/communication"
	                "${OFX_LIBS}/openFrameworks/events"
	                "${OFX_LIBS}/openFrameworks/gl"
	                "${OFX_LIBS}/openFrameworks/graphics"
	                "${OFX_LIBS}/openFrameworks/math"
	                "${OFX_LIBS}/openFrameworks/sound"
	                "${OFX_LIBS}/openFrameworks/types"
	                "${OFX_LIBS}/openFrameworks/utils"
	                "${OFX_LIBS}/openFrameworks/video"
	                "${OFX_MSAOPENCL}"
	                "${OFX_GUI}"
	                "${DEP_CAIRO}/include/cairo"
	                "${DEP_KISS}/include"
	                "${DEP_OPENCL}/include"
	                "${DEP_POCO}/include"
	                "${DEP_TESS2}/include")

if (UNIX)
	find_package(PkgConfig REQUIRED)
	pkg_check_modules(GLIB2 REQUIRED glib-2.0)
	include_directories("${SYSTEM_INCLUDE}/glib-2.0")
	include_directories("/usr/lib/i386-linux-gnu/glib-2.0/include")
	include_directories("${SYSTEM_INCLUDE}/gstreamer-0.10")
endif (UNIX)

set (SOURCE_FILES "src/AABB.cpp"
                  "src/main.cpp"
                  "src/ofApp.cpp"
                  "src/Parameters.cpp"
                  "src/PrefixSum.cpp"
                  "src/Simulation.cpp")

add_executable(pbfSim ${SOURCE_FILES})
