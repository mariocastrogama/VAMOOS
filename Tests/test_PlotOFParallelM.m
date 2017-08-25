% test PlotOFParallelM
clc;
clear;
close all;
fclose all;
format long g;

% load data for analysis
load('results.mat');

% classify with which color
bycolor = 4;

% number of colors of 2D trade-offs
ncolor  = 10;

% recommended palette for this figure [r; w; b]
out_clr = [1.0 0.2 0; 1.0 1.0 1.0; 0.0 0.2 1.0];

% brushing interval 
xbrush = [0.10, 0.25];

% call the function
[h4, idx4] = PlotOFParallelM(objs,bycolor,{},ncolor,out_clr,xbrush);
