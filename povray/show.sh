#bin/bash

rm lab2_data.inc lab2.png
cp ../lab2_data.inc lab2_data.inc
povray +Ilab2.pov +Q11 +A0.10
xdg-open lab2.png
