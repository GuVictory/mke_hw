cmake_minimum_required(VERSION 3.17)
project(dz2)

set(CMAKE_CXX_STANDARD 17)

add_executable(dz2
        main.cpp
        src/entities/AbstractMatrix.cpp
        src/entities/Point3D.cpp
        src/entities/Vec.cpp
        src/entities/DenseMatrix.cpp
        src/entities/DokMatrix.cpp
        src/entities/CsrMatrix.cpp
        src/mke/Triangulation.cpp
        src/generators/LupFactor.cpp
        src/generators/LocalEqGen.cpp
        src/generators/GlobalEqGen.cpp
        src/utils/Export.cpp
        src/utils/Inaccuracy.cpp
        src/utils/SparseSolve.cpp)

target_include_directories(dz2 PUBLIC ./include)