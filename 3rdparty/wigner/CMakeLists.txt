add_subdirectory(wigxjpf)
if(FASTWIGNER)
  add_subdirectory(fastwigxj)
endif()

if(FASTWIGNER)
  add_library(wigner ALIAS fastwigxj)
else()
  add_library(wigner ALIAS wigxjpf)
endif()
