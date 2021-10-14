#!/bin/bash

cd KZG_ext
cmake .
make KZG_std

cd ../Virgo
cmake .
make zk_proof
