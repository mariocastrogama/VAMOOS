%
clc;
clear;
close all;
fclose all;
format long g;

% load data for analysis
load('results.mat');

% call the level diagram function
[hh] = level_diagram(objs);