add_executable(ply-lifting-surface ply-lifting-surface.cpp)
target_link_libraries(ply-lifting-surface vortexje)

configure_file(ply-lifting-surface.ply ${CMAKE_CURRENT_BINARY_DIR}/ply-lifting-surface.ply COPYONLY)
