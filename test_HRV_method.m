% test function tradeoff_index
clc;
clear;
close all;
fclose all;
format long g;

load('results.mat');
% objs = rand(1400,5);

[hh] = HRV_method(objs);
  