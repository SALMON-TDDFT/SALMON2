#! /bin/sh
find ./ -mindepth 1 -maxdepth 1 -type l | xargs -r unlink
find ./ -mindepth 1 -maxdepth 1 ! -name inputfile ! -name preparation ! -name verification ! -name "*.sh" ! -name CMakeFiles ! -name "*.cmake" ! -name Makefile | xargs -r rm -rf
./preparation
