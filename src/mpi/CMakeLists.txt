include_directories(.)
include_directories(..)

# binaries
add_executable(spmspm spmspm.c)
set_target_properties(spmspm PROPERTIES OUTPUT_NAME spmspm)
target_link_libraries(spmspm ${MPI_LIBRARIES} gooseberry m)
install(TARGETS spmspm RUNTIME DESTINATION bin)

add_executable(spmv spmv.c)
set_target_properties(spmv PROPERTIES OUTPUT_NAME spmv)
target_link_libraries(spmv ${MPI_LIBRARIES} gooseberry m)
install(TARGETS spmv RUNTIME DESTINATION bin)

