% test function tradeoff_index
clc;
clear;
close all;
fclose all;
format long g;

load('results.mat');

[ hh ] = tradeoff_index(objs);