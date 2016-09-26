% test PlotObjectives3DMv2 
clc;
clear;
close all;
fclose all;
format long g;

% load data for analysis
load('results.mat');

% classify with which color
bycolor = 4;

% classify with which sizes
bysize  = 5;

% number of colors of 3D trade-offs
ncolor  = 10;

% number of sizes of 3D trade-offs
nsizes  = 6;

% recommended palette for this figure [r; w; b]
out_clr = [1.0 0.2 0; 1.0 1.0 1.0; 0.0 0.2 1.0];

% brushing interval 
xbrush = [0.10, 0.25];

% call the function
[h3, idx3] = PlotObjectives3DMv2(objs,[],[],bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
