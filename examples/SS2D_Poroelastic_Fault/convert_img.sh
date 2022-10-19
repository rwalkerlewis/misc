#!/usr/bin/bash

# Conversion script for poroelastic strike slip result
#convert step01.png -crop 1210x1460+430+90 out.png
convert cauchy_stress_xy.png -crop 1175x1320+750+150 cauchystress_xy.png
