function [hh] = PlotObjectives3DM(varargin)
% function [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize)
%
% Plots the scatters trade-off between Objective Functions obtained as 3D plots
% It helps to visualize any hiden trade-off not visible in low dimensional
% Pareto fronts
%
% When plotting minimization of Objectives it is useful to see if there are
% tradeoffs among the Objectives. One way is to visualize a 3D scatter plot
% of each combination of Objective Functions, then a color palette is
% applied as a function of a 4th Objective Function and finally the sizes
% of the dots are represented by a 5th Objective Function
%
% Input Arguments
%  XX1      : matrix of size [ndata, nobjs]
%  OFnames  : list of names of the OF's
%  OFscales : Scales of each of the time series
%  bycolor  : number of the column to order by color
%  bysize   : number of the column to order by size 
%
% Output Arguments
%  hh       : Handle of the figure
% 
% Other requirements 
%  subtightplot.m
%  create_colors.m
%  create_markers.m
%  create_sizes.m
%  create_fignames.m
%  select_screen.m
%  nchoosek.m  "toolbox\...\specfun"
%  unidrnd.m   "toolbox\stats"
%  sortrows.m  "toolbox\...\datafun" 
%
% Created by: 
%   Mario Castro Gama
%   m.castrogama@unesco-ihe.org
%   PhD Researcher IWSG, UNESCO-IHE
%   Last Update: 2015-11-16
%
% Test
%  [h3] = PlotObjectives3DM;
%
  [out_mrk] = create_markers(1);
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  ifig = 4;
  legend_sizes = {};
  
  % select which monitor to use 1 or 2
  select_screen(1);
%   set(gcf,'Position',[65 65 1520 755]);
  switch nargin
    case 0;
      disp(' Running test [500x5]');
      disp('   XX1 = [1.01+3.55*rand(500,1), 3.14+2.01*rand(500,1), rand(500,1), -1.85+2.03*rand(500,1)];');
      disp('   [hh] = PlotObjectives3DM(XX1);');
      disp(' ');
      XX1 = [1.01+3.55*rand(500,1), 3.14+2.01*rand(500,1), rand(500,1), -1.85+2.03*rand(500,1)];
      [hh] = PlotObjectives3DM(XX1);
    case 1;
      XX1 = varargin{1};
      nobj = size(XX1,2);
      [OFnames, OFscales] = create_fignames(nobj,'obj');
      bycolor  = unidrnd(nobj);
      bysize = unidrnd(nobj);
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize);
    case 2;
      XX1 = varargin{1};
      OFnames  = varargin{2};
      nobj = size(XX1,2);
      [~, OFscales] = create_fignames(nobj,'obj');
      bycolor  = unidrnd(nobj);
      bysize = unidrnd(nobj);
      [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize);
    case 3;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      nobj = size(XX1,2);
      if nobj > length(OFnames);
        error(' OFnames are incomplete');
      end
      if nobj > length(OFscales);
        error(' OFscales are incomplete');
      end
      bycolor  = unidrnd(nobj);
      bysize = bycolor;
      while (bysize == bycolor)
        bysize = unidrnd(nobj);
      end
      [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize);
    case 4;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      nobj = size(XX1,2);
      if nobj > length(OFnames);
        error(' OFnames are incomplete');
      end
      if nobj > length(OFscales);
        error(' OFscales are incomplete');
      end
      if (nobj > bycolor) || (bycolor < 1)
        error(' bycol is not a feasible value');
      end
      bysize = bycolor;
      while (bysize == bycolor)
        bysize = unidrnd(nobj);
      end
      [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize);
    case 5;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      [~, nobj] = size(XX1);
      if isempty(OFnames);
        [OFnames, ~] = create_fignames(nobj,'obj');
      end
      if isempty(OFscales);
        [~, OFscales] = create_fignames(nobj,'obj');
      end
      if nobj > length(OFnames);
        error(' OFnames are incomplete');
      end
      if nobj > length(OFscales);
        error(' OFscales are incomplete');
      end
      if (nobj < bycolor) || (bycolor < 1)
        error(' bycol is not a feasible value');
      end
      if (nobj < bysize) || (bysize < 1)
        error(' bysize is not a feasible value');
      end
      switch nobj
        case {0,1,2};
          disp('Too few objectives to display');
        case {3, 4, 5, 6, 7};
          disp(['by Color : OF',num2str(bycolor)]);
          disp(['by Size  : OF',num2str(bysize)]);
          ifigs = nchoosek(1:nobj,3);
          n = size(ifigs,1);
          ncols = ceil(sqrt(n));
          nrows = ceil(n/ncols);
          disp(['Figure with nrows : ',num2str(nrows),' and ncolumns : ',num2str(ncols)]);
          make_it_tight = true; % make_it_tight = false;
          %                  subtightplot(m, n, p,           gap,       marg_h,       marg_w,  varargin)
          subplot = @(m,n,p) subtightplot(m, n, p, [0.055 0.065], [0.075 0.05], [0.065 0.125]);
          if ~make_it_tight;
            clear subplot;
          end
          
          hh = figure(ifig);
          set(hh,'Position',[62 1 1539 833]);
          set(hh,'Color',[1.0 1.0 1.0]);
          ndata = size(XX1,1);
          xmin = min(XX1);
          xmax = max(XX1);
          [out_clr] = create_colors(6);
          [out_clr] = create_colors(ndata, out_clr);
          XX1 = sortrows(XX1,bycolor); % needs to be ordered by color in the selected column
          
          nsizes = 6;
          [out_sizes, out_sizestick] = create_sizes(XX1,bysize,nsizes);
          isub = 0;
          for irow = 1:nrows
            for icol = 1:ncols
              isub = isub + 1;
              if isub <= n;
                subplot(nrows,ncols,isub);
                x1 = ifigs(isub,1);
                y1 = ifigs(isub,2);
                z1 = ifigs(isub,3);
                % Draw the series once to get the proper legend with the
                % corresponding sizes. Plot with no color.
                if (isub == n); % only do it on the last subplot
                  for iclus = 1:nsizes;
                    xclus = find(out_sizes==iclus,1); % only find the first value with that marker size
                    if ~isempty(xclus);
                      plot3(XX1(xclus,x1),XX1(xclus,y1),XX1(xclus,z1),'Marker',out_mrk{1},'lineStyle','No',...
                        'Markersize',iclus,'MarkerEdgecolor','k','Markerfacecolor','No'); hold on;
                      xleg_a = sprintf('%3.2f',out_sizestick(iclus));
                      xleg_b = sprintf('%3.2f',out_sizestick(iclus+1));
                      legend_sizes = cat(2,[xleg_a,' - ',xleg_b],legend_sizes);
                    end % ~isempty(xclus)
                  end % iclus
                end % (isub == n);
  
                % The second time we do the classification is when we really
                % draw each of the data points. Plot with color.
                for iclus = 1:nsizes;
                  xclus = find(out_sizes==iclus); % find all values with that marker size
                  ndataclus = length(xclus);
                  for idata = 1:ndataclus;
                    if ndataclus > 0;
                      plot3(XX1(xclus(idata),x1),XX1(xclus(idata),y1),XX1(xclus(idata),z1),'Marker',out_mrk{1},'lineStyle','No',...
                      'Markersize',iclus,'MarkerEdgecolor','k','Markerfacecolor',out_clr(xclus(idata),:)); hold on;
                    end % ndataclus > 0
                  end % idata
                end % iclus
              end % isub <= n

              % Set current axis properties
              set(gca,'Fontname',sel_fontname);
              set(gca,'Fontsize',sel_fontsize);
              set(gca,'Fontweight',sel_fontweight);
              set(gca,'xlim', 0.01*round(100*[xmin(x1) xmax(x1)]));
              set(gca,'xtick',0.01*round(100*[xmin(x1) 0.5*(xmax(x1)+xmin(x1)) xmax(x1)]));
              xlabel(OFnames(x1));
              set(gca,'XScale',OFscales{x1});
              set(gca,'ylim', 0.01*round(100*[xmin(y1) xmax(y1)]));
              set(gca,'ytick',0.01*round(100*[xmin(y1) 0.5*(xmax(y1)+xmin(y1)) xmax(y1)]));
              ylabel(OFnames(y1));
              set(gca,'YScale',OFscales{y1});
              set(gca,'zlim', 0.01*round(100*[xmin(z1) xmax(z1)]));
              set(gca,'ztick',0.01*round(100*[xmin(z1) 0.5*(xmax(z1)+xmin(z1)) xmax(z1)]));
              zlabel(OFnames(z1));
              set(gca,'ZScale',OFscales{z1});
              view(-20,10); 
              axis square;
              grid on;
            end % icol
          end % irow

          % Set colormap
          colormap(out_clr);
          xclrtick = (0:0.20:1.0);
          xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;

          % Create colorbar ticks
          xclrlabel = {};
          for ii=1:length(xclrtick);
            xclrlabel = cat(2,xclrlabel, sprintf('%4.2f',xclrbar(ii)));
          end % ii

          % Set colorbar
          hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'Location','manual','Position',[0.930 0.55 0.025 0.20]);

          % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          %  something like...
          % get(hc,'xlabel');
          % set(get(hc,'xlabel'),'String',['by Color : OF_',num2str(bycol)]);
          text('Parent', hc.DecorationContainer, ...
                      'String', ['by Color : OF_',num2str(bycolor)], ...
                      'FontSize', 8.0,...
                      'Fontweight','Bold',...
                      'Position', [0.15, 1.05, 0], ...
                      'Units', 'normalized', ...  
                      'HorizontalAlignment', 'left', ...
                      'VerticalAlignment', 'bottom');
          
          % Set legend. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          %  something like...
          % text('String',['by Size : OF_',num2str(bysize)],'Position',<...>);
          leg = legend(legend_sizes,'Box','Off','Position',[0.90 0.275 0.09 0.13]);
          text('Parent', leg.DecorationContainer, ...
                      'String', ['by Size : OF_',num2str(bysize)], ...
                      'FontSize', 8.0,...
                      'Fontweight','Bold',...
                      'Position', [0.25, 1.05, 0], ...
                      'Units', 'normalized', ...  
                      'HorizontalAlignment', 'left', ...
                      'VerticalAlignment', 'bottom');
          drawnow;
        otherwise
          error(' Too many objectives to display, MAX = 7');
      end % switch nobj
      hold off;
      drawnow;
    otherwise
      error(' Too many input arguments');
  end % switch nargin
end % function