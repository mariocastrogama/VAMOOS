% test function tradeoff_index
clc;
clear;
close all;
fclose all;
format long g;

% load data for analysis
load('results.mat');

% number of colors of masaic of trade-off
ncolor  = 12;

% recommended palette for this figure [r; w; b]
out_clr = [1.0 0.2 0; 1.0 1.0 1.0; 0.0 0.2 1.0]; 

% call the trade-off index function
[hh, lambdas, to_names] = tradeoff_index(objs,ncolor,out_clr);
