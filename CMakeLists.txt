cmake_minimum_required(VERSION 3.0.2)

PROJECT(vtkNonRigidICP)

# Glob all the program file types together
file(GLOB FILES *.h *.cxx)

# Find our packages
find_package(VTK REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

add_executable(vtkNonRigidICP WIN32 ${FILES})
include(${VTK_USE_FILE})
target_link_libraries(vtkNonRigidICP ${VTK_LIBRARIES} Eigen3::Eigen)

# In msvc turn off console on release config
if(MSVC)
	set_target_properties(vtkNonRigidICP PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:CONSOLE")
	set_target_properties(vtkNonRigidICP PROPERTIES LINK_FLAGS_DEBUG "/SUBSYSTEM:CONSOLE")
	target_compile_definitions(vtkNonRigidICP PUBLIC BUILDDIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:mainCRTStartup")
	# Enable parallel build for msvc
	add_definitions(/MP)
endif()