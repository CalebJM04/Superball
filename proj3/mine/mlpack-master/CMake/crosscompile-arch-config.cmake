# This function provides a set of specific flags for each supported board
# depending on the processor type. The objective is to optimize for size.
# Thus, all of the following flags are chosen carefully to reduce binary 
# footprints.

# Set generic minimization flags for all platforms.
# These flags are the same for all cross-compilation cases and they are
# mainly to reduce the binary footprint.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Os -s -fdata-sections -ffunction-sections")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -fno-unwind-tables")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-asynchronous-unwind-tables -fvisibility=hidden")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fshort-enums -finline-small-functions")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -findirect-inlining -fno-common")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmerge-all-constants -fno-ident")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-unroll-loops -fno-math-errno")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-stack-protector")
set(CMAKE_OPENBLAS_FLAGS "${CMAKE_CXX_FLAGS}") # OpenBLAS does not supoport flto
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--hash-style=gnu -Wl,--build-id=none")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,norelro")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--gc-sections")
## Keep the following flag in comment, they might be relevant in the case of MCU's
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-nmagic,-Bsymbolic -nostartfiles")

# BOARD_NAME is deprecated and will be removed in mlpack 5.
# please use ARCH_NAME instead.
set(BOARD_NAME "" CACHE STRING "Specify Board name to optimize for.")
set(ARCH_NAME "" CACHE STRING "Name of embedded architecture to optimize for.")

if (BOARD_NAME)
  set(ARCH_NAME "${BOARD_NAME}")
endif()

string(TOUPPER ${ARCH_NAME} ARCH)

# Set specific platforms CMAKE CXX flags.
if(ARCH MATCHES "RPI0" OR ARCH MATCHES "RPI1" OR ARCH MATCHES "ARM11")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=arm1176jzf-s")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=arm1176jzf-s -mfloat-abi=hard -mfpu=vfp")
  set(OPENBLAS_TARGET "ARMV6")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "RPI2" OR ARCH MATCHES "CORTEXA7")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a7")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=hard -mfpu=neon-vfpv4")
  set(OPENBLAS_TARGET "ARMV7")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "CORTEXA8")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a8")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=hard -mfpu=neon")
  set(OPENBLAS_TARGET "ARMV7")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "CORTEXA9")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a9")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=hard -mfpu=neon")
  set(OPENBLAS_TARGET "CORTEXA9")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "CORTEXA15")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a15")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfloat-abi=hard -mfpu=neon")
  set(OPENBLAS_TARGET "CORTEXA15")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "RPI3" OR ARCH MATCHES "CORTEXA53")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a53 -ftree-vectorize")
  set(OPENBLAS_TARGET "CORTEXA53")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "RPI4" OR ARCH MATCHES "CORTEXA72")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=armv8.2-a+crypto+fp16+rcpc+dotprod -fasynchronous-unwind-tables")
  set(OPENBLAS_TARGET "CORTEXA72")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "JETSONAGX" OR ARCH MATCHES "CORTEXA76")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=cortex-a76 -ftree-vectorize")
  set(OPENBLAS_TARGET "CORTEXA76")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "BV")
  set(OPENBLAS_TARGET "RISCV64_GENERIC")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "C906")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=thead-c906")
  set(OPENBLAS_TARGET "RISCV64_GENERIC")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "x280")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=sifive-x280")
  set(OPENBLAS_TARGET "x280")
  set(OPENBLAS_BINARY "64")
elseif(ARCH MATCHES "KATAMI")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=pentium3")
  set(OPENBLAS_TARGET "KATAMI")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "COPPERMINE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=pentium3")
  set(OPENBLAS_TARGET "COPPERMINE")
  set(OPENBLAS_BINARY "32")
elseif(ARCH MATCHES "NORTHWOOD")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=pentium4")
  set(OPENBLAS_TARGET "NORTHWOOD")
  set(OPENBLAS_BINARY "32")
elseif(ARCH)
  ## TODO: update documentation with a list of the supported boards.
  message(FATAL_ERROR "Given ARCH_NAME is not known; please choose a supported board from the list")
endif()
