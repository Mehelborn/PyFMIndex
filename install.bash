#!/usr/bin/env bash

git submodule update --init --recursive --remote
cmake -S ./lib/AvxWindowFmIndex/ -B ./lib/AvxWindowFmIndex/
make -C ./lib/AvxWindowFmIndex/
