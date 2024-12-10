cmake_minimum_required(VERSION 3.30)
project(Gas_Sim)

set(CMAKE_CXX_STANDARD 20)


add_executable(Gas_Sim main.cpp
        include/Vector2D.h
        include/Ensemble.h
        include/Wall.h
        include/Particle.h
        include/VWParticle.h
        src/Particle.cpp
        src/Wall.cpp
        src/Ensemble.cpp
        src/VWParticle.cpp
        include/Quad.h
        include/QuadTree.h
        src/QuadTree.cpp
        include/QTEnsemble.h
        src/QTEnsemble.cpp
        include/Constants.h
)

