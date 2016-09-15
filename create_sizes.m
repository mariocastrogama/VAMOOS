function [ out, out_tick ] = create_sizes(XX1,bysize,n)
% function [out, out_tick] = create_sizes(XX1,bysize,n)
% Creates the sizes required based on a single objective function
% If functions are for minimization BIG markersize is given to low values
% then SMALL markersize is given to high values 
%
% Input Arguments:
% XX1    : the array of the Objectives [nsamples, nobjs]
% bysize : index of the column to perform the classification
% n      : Number of clusters or different sizes the final matrix contains
%          values between { 1, ..., n}
%
% Output arguments:
% out    : palette of colors of size [n,3] based on icolor 
%          and linearly interpolated 
% out_tick : Limits of each value
%
% Notes:
% 1) In some cases the variables are distributed with LOG scale
%
% Developed by: 
% Mario Castro Gama
% PhD researcher
% 2015-11-16
% 
% Example:
%
%  [out] = create_sizes(rand(100,3),2,7);
%
  switch nargin
    case 0;
      XX1 = rand(100,3);
      bysize = 1;
      n = 5;
      out = create_sizes(XX1,bysize,n);
    case 1;
      bysize = 1;
      n = 5;
      out = create_size(XX1,bysize,n);
    case 2;
      n = 5;
      out = create_size(XX1,bysize,n);
    case 3;
      [~, ncol] = size(XX1);
      if ~isreal(XX1) || (~isreal(bysize) || ~isreal(n));
        error(' Some arguments are not real values');
      end
      if (n < 1) || (bysize < 1);
        error(' n OR bysize < 1');
      end
      if (iscell(bysize) || iscell(n))
        error(' NON numerical arguments');
      end
      if (bysize > ncol)
        error(' XX1 < bysize, less columns than expected for classification');
      end
      n = round(n);
      if length(n) ~= 1;
        error(' n must be an scalar');
      end
      xvar = XX1(:,bysize); %rescale(XX1(:,bysize));
      xq = (0:1.0/(n):1.0);
      qq1 = quantile(xvar,xq);
      [~,~,out] = histcounts(xvar,qq1); % equally distributed by quantile
      out_tick = qq1;
      out = (n + 1) - out;
      % just based on the data as it is
%       out = (n + 1) - round((n-1)*out + 1); 
    otherwise
      error(' Too many input arguments');
  end
end

