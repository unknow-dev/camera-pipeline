find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(exec/emb ${PYTHON_LIBRARIES})
