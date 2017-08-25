function [hh,idx] = PlotOFParallelM(varargin)
% function [hh]=PlotOFParallelM(XX,bycol,xnames,out_clr,ncolor,xbrush)
%
% Plots the trade-off between objective functions obtained
% for every Objective Function as PARALLEL plot, it allows brushing off
% values which are nto relevant
% 
% 
% Input Arguments:
%   XX1      : matrix [nsamples, ncols] containing the Objective Functions
%             usually in Multi-Objective Optimization all objectives are
%             intended for minimization or maximization.
%   bycolor : organize the arrows by which column
%   OFnames : names of the objective functions
%   out_clr : a palette of colors to use for interpolation [ni x 3]
%   ncolor  : number of colors for plotting
%   xbrush  : [xmin, xmax] two values which define the limits of the
%             variable bycolor, if the limits are violated no brushing is
%             performed
%
% Output Arguments:
%   hh     : Figure handle
%   idx    : list of samples of brushed values
%
% Additional requirements:
%  create_colors.m
%  rescale.m
%
% Notes:
% 1) Need to think how to show when objectives have LOG scale
% 2) Need to implement the possibility to use some objectives as
%    maximization and some as minimization
%
% Developed by:
%   Mario Castro Gama
%   PhD Researcher
%   2015-10-15
%
% Last Update
%   2016-07-10, Included a palette of colors as well (out_clr),
%               allows to include a new variable to count the number of colors
%               (ncolor). Recommended not to be more than 11 colors.
%   2016-08-01, Included the possibility to brush data including a new variable 
%               named xbrush which brings the minimum and maximum of the varaible of
%               selection.
%   2016-09-01, Corrected a few bugs realted with the color values at the
%               colorbar.
%
% % Example
% 
%   XX = [1+4*rand(500,1), 3+2*rand(500,1), rand(500,1), -1+2*rand(500,1)];
%   bycolor = 2;
%   out_clr = [1, 0, 0; 1 1 0; 0 0.65 1];
%   xbrush  = [3.50 3.75];
%   PlotOFParallelM(XX,bycolor,{'OF_1','OF_2','OF_3','OF_4'},10,out_clr,xbrush);

% Options of plot
  optgca.Fontname   = 'Arial';
  optgca.Fontsize   = 12;
  optgca.Fontweight = 'Bold';
  optgca.ytick      = [0.0:0.25:1.0];
  optgca.ylim       = [-0.1   1.1];
  optgca.yticklabel = {'min','','','','max'};
  optgca.Position   = [0.065 0.085 0.85 0.85];

  opttxt.Fontsize   = 12;
  opttxt.Fontweight = 'Bold';
  opttxt.HorizontalAlignment = 'Center';
  opttxt.VerticalAlignment   = 'Middle';
  opttxt.BackgroundColor   = [1 1 1];

  % it indicates the numerical precision of results and format
  ndecimals = 3;
  strformat = ['%4.',num2str(ndecimals),'f'];

  % select which monitor to use 1 or 2
  select_screen(2);

  switch nargin
    case 0;
      disp('Run example');
      disp('XX1 = [1+4*rand(500,1), 3+2*rand(500,1), rand(500,1), -1+2*rand(500,1)];');
      ndata = 500;
      XX1       = [1+4*rand(ndata,1), 3+2*rand(ndata,1), rand(ndata,1), -1+2*rand(ndata,1)];
      OFnames  = {'OF_1','OF_2','OF_3','OF_4'};
      bycolor  = 2;
      ncolor   = 8;
      out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [3.50 4.5];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 1;
      XX1       = varargin{1};
      nobj     = size(XX1,2);
      bycolor  = 2;
      [OFnames, ~] = create_fignames(nobj,'obj');
      ncolor = 8;
      out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 2;
      XX1       = varargin{1};
      nobj     = size(XX1,2);
      bycolor  = varargin{2};
      [OFnames, ~] = create_fignames(nobj,'obj');
      ncolor  = 8;
      out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 3;
      XX1       = varargin{1};
      bycolor  = varargin{2};
      OFnames  = varargin{3};
      ncolor   = 8;
      out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 4
      XX1       = varargin{1};
      bycolor  = varargin{2};
      OFnames  = varargin{3};
      ncolor   = varargin{4};
      out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 5
      XX1       = varargin{1};
      bycolor  = varargin{2};
      OFnames  = varargin{3};
      ncolor   = varargin{4};
      out_clr  = varargin{5};
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotOFParallelM(XX1,bycolor,OFnames,ncolor,out_clr,xbrush);
    case 6
      XX1       = varargin{1};
      nobj     = size(XX1,2);
      bycolor  = varargin{2};
      OFnames  = varargin{3};
      ncolor   = varargin{4};
      out_clr  = varargin{5};
      xbrush   = varargin{6};
      if (isstring(bycolor) || iscell(bycolor)) || (bycolor < 1)
        error('bycol is not a valid number');
      else
        bycolor = round(bycolor);
      end
      if isempty(OFnames);
        [OFnames, ~] = create_fignames(nobj,'obj');
      end
      if (nobj < bycolor)
        error('bycolor > ncols of values');
      end
      if (length(OFnames) < nobj)
        error('length(OFnames) ~= ncols of values');
      end
      if isempty(ncolor)
        ncolor = 8;
      end
      if isempty(out_clr)
        out_clr   = [1, 0, 0; 1 1 0; 0 0.65 1];
        [out_clr] = create_colors(ncolor,out_clr);
      end
      if isempty(xbrush)
        xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      end
      
      disp(' ');
      disp('   Figure Parallel Objectives');
      disp(['    by Color : OF',num2str(bycolor),' ',OFnames{bycolor}]);
      disp(['    xbrush = [',sprintf(strformat,xbrush(1)),', ',sprintf(strformat,xbrush(2)),']']);
      disp(' ');
      
      if (xbrush(2) < min(XX1(:,bycolor)));
        disp('Max value of xbrush < min(X(:,bycolor))');
        disp('Reverting to min and max of bycolor');
        xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      end
      if (xbrush(1) > max(XX1(:,bycolor)));
        disp('min value of xbrush > max(X(:,bycolor))');
        disp('Reverting to min and max of bycolor');
        xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      end
          
      optgca.xtick      = (1:nobj);
      xmove = 0.02*nobj;
      optgca.xlim       = [(1.0-xmove) nobj+xmove];
      
      % generate new xlabels based on OFnames
      optgca.xticklabel = {};
      for s = 1:nobj
        tmp = OFnames{s};
        optgca.xticklabel = cat(2,optgca.xticklabel,tmp);
      end
      clear s;

      % Create colors
      [out_col] = create_colors(ncolor, out_clr);

      XX1   = sortrows(XX1,bycolor);
      % Estimate maximum and minimum of all samples
      xmin = min(XX1);
      xmax = max(XX1);

      % Selection of a valid range between two different values
      sel_brush = find((XX1(:,bycolor) > xbrush(1)) & (XX1(:,bycolor) < xbrush(2)) == 1);
      sel_shadow = find((XX1(:,bycolor) < xbrush(1)) | (XX1(:,bycolor) > xbrush(2)) == 1);
      XXshadow = XX1(sel_shadow,:); 
      XXbrush = XX1(sel_brush,:);

      % Estimate maximum and minimum of selection
      xbrushmin = min(XXbrush);
      xbrushmax = max(XXbrush);

      nbrush = size(sel_brush,1);
      % Check if it finds any solution
      if nbrush == 0;
        error([' The range of the selected variable (bycolor)',...
               ' does not allow any feasible solution.       ']);
      else
        
        XX1   = rescale(XX1);
        XXbrush_scaled = XX1(sel_brush,:);
        XXshadow_scaled = XX1(sel_shadow,:);
        [ idx, ~, ~, ~, ~ ] = create_sizes_v2(XXbrush_scaled,1,1,bycolor,ncolor);
        
        % Plot each of the samples as a single line
%         ifig = 6;
        hh = figure;%(ifig);
        clf;
        set(hh,'Color',[1.0 1.0 1.0]);
        xa = (1:nobj)';
        
        % first plot the shadow values, this guarantees these are behind
        nshadow = size(XXshadow,1);
        if (nshadow ~=0);
          ya = XXshadow_scaled';
          plot(xa,ya,'color',[0.85 0.85 0.85]);
          hold on;
        end
        
        % start from high values (low fitness) to low values (best fitness)
        for icolor = ncolor:-1:1; 
          xrange = idx{1,icolor};
          if ~isempty(xrange)
            ya = XXbrush_scaled(xrange,:)';
            plot(xa,ya,'color',out_col(icolor,:));
            hold on;
          end % ~isempty(xrange)
        end % icolor
        
        % extract the id's of the samples which were brushed IN
        if nbrush < size(XX1,1);
          for ii = 1:length(idx);
            ki        = length(idx{ii});
            if (ki > 0)
              idx{ii} = sel_brush(1:ki);
              sel_brush(1:ki) = [];
            end
          end
        end
        
        % Create colormap and colorbar ticks besed on the Column of visualization
        colormap(out_col);
        xclrtick  = (0:(1.0/ncolor):1.0);
        xclrbar   = xbrushmin(bycolor) + (xbrushmax(bycolor)-xbrushmin(bycolor))*xclrtick;
        xclrlabel = {};
        for ii=1:length(xclrtick);
          xclrlabel = cat(2,xclrlabel, sprintf(strformat,xclrbar(ii)));
        end % ii
        
        % Set colorbar
        hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'Location','manual','Position',[0.930 0.15 0.025 0.650]);
        text('Parent', hc.DecorationContainer, ...
          'String', {'by Color OF : ';OFnames{bycolor}}, ...
          'FontSize', opttxt.Fontsize-2,...
          'Fontweight','Bold',...
          'Position', [0.25, 1.01, 0], ...
          'Units', 'normalized', ...
          'HorizontalAlignment', 'left', ...
          'VerticalAlignment', 'bottom');
        
        % Print the minimum and maximum of each Objective function
        for icol = 1: nobj
          text(icol,-0.05,sprintf(strformat,xmin(icol)),opttxt);
          text(icol, 1.05,sprintf(strformat,xmax(icol)),opttxt);
        end % icol
        
        % Set options of axis
        set(gca,optgca);
        grid on;
        % Set title of figure
        title(['PARALLEL PLOT BY ',OFnames{bycolor}]);
        hold off;
      end % size(sel_brush,1) == 0;
    otherwise % nargin
      error(' Too many input arguments');
  end % nargin
end
