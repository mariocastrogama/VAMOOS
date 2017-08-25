function [out]=create_colors(n,icolor)
% function [out]=create_colors(n)
% Creates a palette for plotting purposes using colorbars and surfaces
%
% Input Arguments:
% n      : number of required colors
% icolor : a palette of RGB of size [kcolors, 3], 
%          all colors between 0 and 1.
%
% Output arguments:
% out    : palette of colors of size [n,3] based on icolor 
%          and linearly interpolated 
%
% Notes:
% 1) In some journals they specifically ask for black & white this should
%    be included
% 2) In some cases the variables are distributed with LOG scale
%
% Developed by: 
% Mario Castro Gama
% PhD researcher
% 2015-10-01
% 
% Example:
%   icolor = [1 0 0;... % red
%             0 0 1;... % blue                  
%             0 1 0;... % green
%             0 1 1;... % cyan         
%             1 0 1;... % magenta
%             1 1 0;... % yellow        
%             0 0 0];... % black
%
%  [out] = create_colors(20,icolor);
%  surf(sort(rand(100,100))); 
%  shading flat; 
%  view(0,90);
%  colormap(out);
%  colorbar;
%

  switch nargin
    case 0;
      out = create_colors(12);
    case 1;
%       icolor = [0 0 0; 1 1 1]; % Black and White
      icolor = [1.0 0.0 0.0;
                1.0 1.0 0.0;
                0.0 1.0 0.0;
                0.0 1.0 1.0;
                0.0 0.0 1.0];
%       icolor = [0.0000 0.0000 1.0000;
%                 0.5000 0.5000 1.0000;
%                 0.5500 1.0000 0.5500;
%                 0.9000 1.0000 0.1000;
%                 1.0000 0.8000 0.1000;
%                 1.0000 0.3000 0.3000;
%                 0.0000 1.0000 0.0000;
%                 0.0000 0.3557 0.6747;
%                 0.4709 0.0000 0.0180;
%                 0.8422 0.1356 0.8525;
%                 0.4688 0.6753 0.3057;
%                 0.0000 0.0000 0.0000;
%                 1.0000 1.0000 1.0000];
    out = create_colors(n,icolor);
    case 2;
      if (~isreal(icolor) || ~isreal(n));
        error(' Some arguments are not real values');
      end
      if n < 0;
        error(' n < 0');
      end
      if (iscell(icolor) || iscell(n))
        error(' NON numerical arguments');
      end
      if (max(icolor(:)) > 1.0) || (min(icolor(:))<0.0);
        error(' One or more colors have invalid values outside [0.0, 1.0]');
      end
      n = round(n);
      ncolor = size(icolor,1);
      if n ==1;
        out = [0 0 1];
      else
        if n < ncolor
          out = icolor(1:n,:);
        else
          Xstar  = (1.0 : (ncolor-1)/(n-1): ncolor);
          out = interp1(1:ncolor, icolor, Xstar);
        end
      end
    otherwise
      error(' Too many input arguments');
  end
end