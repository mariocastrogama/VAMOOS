function [ idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor)
% function [out, out_tick] = create_sizes_v2(XX1,bysize,n)
% Creates the sizes required based on a single objective function
% If functions are for minimization BIG markersize is given to low values
% then SMALL markersize is given to high values 
%
% Input Arguments:
%  XX1    : the array of the Objectives [nsamples, nobjs]
%  bysize : index of the column to perform the classification by size
%  nsizes : Number of clusters or different sizes the final matrix contains
%           values between { 1, ..., nsizes}
%  bycolor: index of the column to perform the classification by color
%  ncolor : Number of clusters or different colors the final matrix contains
%          values between { 1, ..., nnolor}
%
% Output arguments:
%  out    : palette of colors of size [n,3] based on icolor 
%           and linearly interpolated 
%  out_tick : Limits of each value
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
%  [out] = create_sizes_v2(rand(100,3),2,7);
%
  switch nargin
    case 0;
      XX1 = rand(100,3);
      bysize = 1;
      nsizes = 5;
      [idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
    case 1;
      bysize = 1;
      nsizes = 5;
      bycolor = 2;
      ncolor = 10;
      [idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
    case 2;
      nsizes = 5;
      [idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
    case 3;
      [~, ncol] = size(XX1);
      if ~isreal(XX1) || (~isreal(bysize) || ~isreal(nsizes));
        error(' Some arguments are not real values');
      end
      if (nsizes < 1) || (bysize < 1);
        error(' n OR bysize < 1');
      end
      if (iscell(bysize) || iscell(nsizes))
        error(' NON numerical arguments');
      end
      if (bysize > ncol)
        error(' XX1 < bysize, less columns than expected for classification');
      end
      nsizes = round(nsizes);
      if length(nsizes) ~= 1;
        error(' n must be an scalar');
      end
      % new code 2016-06-01
      bycolor = bysize;
      while bycolor == bysize;
        bycolor = randi(ncol);
      end
      ncolor = 10;
      [idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
    case 5      
      [~, ncol] = size(XX1);
      if ~isreal(XX1) || (~isreal(bysize) || ~isreal(nsizes));
        error(' Some arguments are not real values');
      end
      if (nsizes < 1) || (bysize < 1);
        error(' n OR bysize < 1');
      end
      if (iscell(bysize) || iscell(nsizes))
        error(' NON numerical arguments');
      end
      if (bysize > ncol)
        error(' XX1 < bysize, less columns than expected for classification');
      end
      nsizes = round(nsizes);
      if length(nsizes) ~= 1;
        error(' n must be an scalar');
      end
      
      var_size  = XX1(:,bysize);
      var_color = XX1(:,bycolor);
      q_sizes   = quantile(var_size, (0.0 : 1.000/(nsizes) : 1.000));
      q_colors  = quantile(var_color,(0.0 : 1.000/(ncolor) : 1.000));

      [~, tick_size, tick_color, out_sizes, out_colors] = histcounts2(var_size, var_color, q_sizes, q_colors);
           
      %% here is the extraction of the classification bycolor
      idx = cell(nsizes,ncolor);
      for isize = 1:nsizes;
        for icolor = 1:ncolor;
          tmp = find(out_sizes == isize);
          idx{isize,icolor} = tmp(find(out_colors(tmp) == icolor));
        end
      end
    otherwise
      error(' Too many input arguments');
  end
end

