#!/bin/bash

nohup ./langevin ~/3dDAT/ ~/3dIMG300/ 300 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG350/ 350 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG400/ 400 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG450/ 450 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG500/ 500 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG550/ 550 1 10 1 1 22000 &
sleep 1s
nohup ./langevin ~/3dDAT/ ~/3dIMG600/ 600 1 10 1 1 22000 &
sleep 1s



