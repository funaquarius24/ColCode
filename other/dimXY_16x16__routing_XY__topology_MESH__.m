% fname: dimXY_16x16__routing_XY__topology_MESH__.m
% ../bin/noxim  -routing XY -topology MESH -dimx 16 -dimy 16  -sim 10000 -warmup 2000 -size 8 8 -buffer 4 -config ../configs/mcsl_default.yaml -power ../bin/power.yaml 

function [max_pir, max_throughput, min_delay] = dimXY_16x16__routing_XY__topology_MESH__(symbol)

data = [
%             pir      avg_delay     throughput      max_delay   total_energy       rpackets         rflits
