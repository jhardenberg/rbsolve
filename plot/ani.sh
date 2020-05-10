#!/bin/bash
fname=rbpat15
ffmpeg -f gif -i $fname.gif -pix_fmt yuv420p -c:v libx264 -movflags +faststart -filter:v crop='floor(in_w/2)*2:floor(in_h/2)*2' $fname.mov
