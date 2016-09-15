function select_screen(id)
% function select_screen(id)
%
% This function allows you to select which screen you want to define as default for plotting figures 
% and sets the figure in the limits on the screen (fullsize).
% Useful when working with two screens so you can decide where to send the figure to be plotted.
% 
% Input Arguments
%  id : Integer number of the selected screen
%
% Notes:
% 1) In order to work in all figures, you must close first the open figures
%   >> close all
%
% MSc Mario Castro Gama
% PhD Researcher 
% UNESCO-IHE, IWSG-HI
% 2015-12-01
%
  
  % Get sizes of screens in pixels
  % screen_sizes = get(0,'MonitorPositions');
  screen_sizes = [59 1 1542 833; 1601 -179 1920 1013];
  % how many screens are available?
  nscreens = size(screen_sizes,1);
  switch nargin 
    case 0;
      id = 1; % default screen
      select_screen(id);
    case 1;
      if id > nscreens;
        error(['there is(are) ',num2str(nscreens),' screen(s) and you picked screen ',num2str(id)]);
      else
        if nscreens == 1; % Screen one by default
          set(0,'DefaultFigurePosition',screen_sizes);
        else
          set(0,'DefaultFigurePosition',screen_sizes(id,:));
        end
      end
    otherwise
      error('Too many input arguments');
  end % nargin
end