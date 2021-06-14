set(BUILD_SHARED_LIBS OFF)
add_subdirectory("toml-f")
list(
  APPEND lib-deps
    "toml-f-lib"
)
set(lib-deps "${lib-deps}" PARENT_SCOPE)
