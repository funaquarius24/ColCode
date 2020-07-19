% fname: dimXY_8x8__routing_DYAD__topology_MESH__.m
% ../bin/noxim  -routing DYAD -topology MESH -dimx 8 -dimy 8  -sim 100000 -warmup 2000 -size 8 8 -buffer 4 -config ../configs/mcsl_default.yaml -power ../bin/power.yaml 

function [max_pir, max_throughput, min_delay] = dimXY_8x8__routing_DYAD__topology_MESH__(symbol)

data = [
%             pir      avg_delay     throughput      max_delay   total_energy       rpackets         rflits
