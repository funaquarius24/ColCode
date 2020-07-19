% fname: dimXY_8x8__routing_NORTH_LAST__topology_MESH__.m
% ../bin/noxim  -routing NORTH_LAST -topology MESH -dimx 8 -dimy 8  -sim 100000 -warmup 2000 -size 8 8 -buffer 4 -config ../configs/mcsl_default.yaml -power ../bin/power.yaml 

function [max_pir, max_throughput, min_delay] = dimXY_8x8__routing_NORTH_LAST__topology_MESH__(symbol)

data = [
%             pir      avg_delay     throughput      max_delay   total_energy       rpackets         rflits
             0.01         9.4295      0.0294617            244     8.0778e-05          23099         184784
             0.04        1768.83      0.0740108          42109    8.35155e-05          58027         464196
];

rows = size(data, 1);
cols = size(data, 2);

data_delay = [];
for i = 1:rows/1,
   ifirst = (i - 1) * 1 + 1;
   ilast  = ifirst + 1 - 1;
   tmp = data(ifirst:ilast, cols-6+1);
   avg = mean(tmp);
   [h sig ci] = ttest(tmp, 0.1);
   ci = (ci(2)-ci(1))/2;
   data_delay = [data_delay; data(ifirst, 1:cols-6), avg ci];
end

fig1 = figure(1);
set(fig1,'Name','data_delay');
hold on;
plot(data_delay(:,1), data_delay(:,2), symbol);

data_throughput = [];
for i = 1:rows/1,
   ifirst = (i - 1) * 1 + 1;
   ilast  = ifirst + 1 - 1;
   tmp = data(ifirst:ilast, cols-6+2);
   avg = mean(tmp);
   [h sig ci] = ttest(tmp, 0.1);
   ci = (ci(2)-ci(1))/2;
   data_throughput = [data_throughput; data(ifirst, 1:cols-6), avg ci];
end

fig2 = figure(2);
set(fig2,'Name','data_throughput');
hold on;
plot(data_throughput(:,1), data_throughput(:,2), symbol);

data_maxdelay = [];
for i = 1:rows/1,
   ifirst = (i - 1) * 1 + 1;
   ilast  = ifirst + 1 - 1;
   tmp = data(ifirst:ilast, cols-6+3);
   avg = mean(tmp);
   [h sig ci] = ttest(tmp, 0.1);
   ci = (ci(2)-ci(1))/2;
   data_maxdelay = [data_maxdelay; data(ifirst, 1:cols-6), avg ci];
end

fig3 = figure(3);
set(fig3,'Name','data_maxdelay');
hold on;
plot(data_maxdelay(:,1), data_maxdelay(:,2), symbol);

data_totalenergy = [];
for i = 1:rows/1,
   ifirst = (i - 1) * 1 + 1;
   ilast  = ifirst + 1 - 1;
   tmp = data(ifirst:ilast, cols-6+4);
   avg = mean(tmp);
   [h sig ci] = ttest(tmp, 0.1);
   ci = (ci(2)-ci(1))/2;
   data_totalenergy = [data_totalenergy; data(ifirst, 1:cols-6), avg ci];
end

fig4 = figure(4);
set(fig4,'Name','data_totalenergy');
hold on;
plot(data_totalenergy(:,1), data_totalenergy(:,2), symbol);


%-------- Saturation Analysis -----------
slope=[];
for i=2:size(data_throughput,1),
    slope(i-1) = (data_throughput(i,2)-data_throughput(i-1,2))/(data_throughput(i,1)-data_throughput(i-1,1));
end

for i=2:size(slope,2),
    if slope(i) < (0.95*mean(slope(1:i)))
        max_pir = data_throughput(i, 1);
        max_throughput = data_throughput(i, 2);
        min_delay = data_delay(i, 2);
        break;
    end
end
