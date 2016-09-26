% test PlotOFParallelM
clc;
clear;
close all;
fclose all;
format long g;

% load data for analysis
load('results.mat');

% by Objective
byOF = 2;

% names OF
[OFnames, ~] = create_fignames(5,'obj');

% names DV
[DVnames, ~] = create_fignames(3,'dec');

% Objective type 
objtype = {'min','min','min','min','min'};

% call the function
[h] = OFspace(objs,vars,byOF,OFnames,DVnames,objtype);
