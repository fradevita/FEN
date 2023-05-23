#!/bin/bash

for d in ./*/ ; do (cd "$d" && sh clean.sh); done
