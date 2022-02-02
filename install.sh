#!/bin/bash

# compile attribute
Rscript comp_att.R

# install
cd ..
R CMD INSTALL --preclean --no-multiarch --with-keep.source hubbert
