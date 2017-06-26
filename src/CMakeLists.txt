include_directories(.)
file(GLOB gooseberry_sources *.c)
list(REMOVE_ITEM gooseberry_sources "gooseberry_bin.c")
file(GLOB domlib_sources ${CMAKE_SOURCE_DIR}/${DOMLIB_PATH}/*.c)

# library                                         
add_library(gooseberry ${GOOSEBERRY_LIBRARY_TYPE} 
  ${gooseberry_sources}
  ${domlib_sources}
) 
install(TARGETS gooseberry 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# binary
add_executable(gooseberry_bin gooseberry_bin.c)
set_target_properties(gooseberry_bin PROPERTIES OUTPUT_NAME gooseberry)
target_link_libraries(gooseberry_bin gooseberry m)
install(TARGETS gooseberry_bin
  RUNTIME DESTINATION bin
)

if (NOT DEFINED NO_MPI OR NO_MPI EQUAL 0)
  add_subdirectory("mpi")
endif()
