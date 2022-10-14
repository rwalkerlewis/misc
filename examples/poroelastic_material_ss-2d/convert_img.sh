#!/usr/bin/bash

# Conversion script for poroelastic strike slip result
#convert step01.png -crop 1210x1460+430+90 out.png
#convert cauchy_stress_xy_poropng -crop 1245x1465+760+90 cauchystress_poro_xy.png
convert cauchy_stress_xy_poro.png -crop 1245x1470+750+80 cauchystress_poro_xy.png
