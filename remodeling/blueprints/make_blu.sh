#!/bin/bash

for f in $(ls one_each)
do
    /home/svensken/Rosetta/main/source/bin/make_blueprint.default.linuxgccrelease -database /home/svensken/Rosetta/main/database/ -s one_each/$f
    mv default.blueprint ${f:0:9}.blueprint
done
