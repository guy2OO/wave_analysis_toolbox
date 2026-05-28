function [] = Wave_Theory(height,depth,period)

%% Imputs
%
%   height - wave height (m)
%
%   depth - water depth in wave region (m)
%
%   period - wave period (s)
%

%% Define gravity and import wave date to table and background image
g = 9.81;
bg = imread("wave_graph.png");

%Process wave data for each row to get x and y axis values
y_axis = height ./ (g .* (period .^2));
x_axis = depth ./ (g .* (period .^2));

%Plot data on the same figure
xrange = [0.0001 1];
yrange = [0.00001 0.1];

f = figure;
hold on
xscale("log");
yscale("log");
imagesc(xrange,yrange,bg,AlphaData=0.9);
title('Stokes 2nd Order Wave Theory');
xlabel('d/gT^2');
ylabel('H/gT^2');
labels = {'1','2','3'};

loglog(x_axis,y_axis,'s','color','r','Marker','*','LineWidth',1);
text(x_axis,y_axis,labels,'VerticalAlignment','bottom','HorizontalAlignment','right');
xlim(xrange);
ylim(yrange);
axis square;

% Boundry markers for three depths
    %loglog(0.068,0.01,'s','color','b','marker','*')
    %loglog(0.00165,0.01,'s','color','b','marker','*')

hold off

% Save figure as .png with current date and time
date = datetime("now");
name = sprintf("WaveTheory_%s.png",date);
exportgraphics(f,name)