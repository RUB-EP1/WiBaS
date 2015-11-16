#!/bin/tcsh

./testApp -f 1 -l 1000 !>& wibas1.log &
./testApp -f 1001 -l 2000 !>& wibas2.log &
./testApp -f 2001 -l 3000 !>& wibas3.log &


