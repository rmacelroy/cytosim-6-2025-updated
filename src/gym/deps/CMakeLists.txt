
set(SOURCES_GYM_DEPS
    spng.c
    miniz.c
)

list(TRANSFORM SOURCES_GYM_DEPS PREPEND "${PROJECT_SOURCE_DIR}/src/gym/deps/")

set_property(SOURCE spng.c PROPERTY COMPILE_DEFINITIONS SPNG_STATIC SPNG_USE_MINIZ)
