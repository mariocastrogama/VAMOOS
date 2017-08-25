function [hh] = PlotObjectivesM(varargin)
% function [hh] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize)
%
% Plots the scatters trade-off between Objective Functions obtained as 2D plots
% It helps to visualize any hiden trade-off not visible in high dimensional
% Pareto fronts
%
% When plotting minimization of Objectives it is useful to see if there are
% tradeoffs among the Objectives. One way is to visualize a 2D scatter plot
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
%                2016-06-06, fixed a bug with the number of decimals for
%                display in legend and colorbar.
%
% Test
%  [h3] = PlotObjectivesM;
%
  [out_mrk] = create_markers(1);
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  ifig = 3;
  legend_sizes = {};
  
  % this can be set to plot different variables
  % it indicates the numerical precision of results
  ndecimals = 3; 
  Rmax = 10^ndecimals;
  Rmin = 1/Rmax;
  strformat = ['%4.',num2str(ndecimals),'f'];
  
  % select which monitor to use 1 or 2
  select_screen(1);
  
  switch nargin
    case 0;
      disp(' % Running test [500x5]');
      disp('   XX1 = [1.01+3.55*rand(500,1), 3.14+2.01*rand(500,1), rand(500,1), -1.85+2.03*rand(500,1)];');
      disp('   [hh] = PlotObjectivesM(XX1);');
      disp(' ');
      nsub = 250;
      XX1 = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1)];
      [hh] = PlotObjectivesM(XX1);
    case 1;
      XX1 = varargin{1};
      nobj = size(XX1,2);
      [OFnames, OFscales] = create_fignames(nobj,'obj');
      bycolor = unidrnd(nobj);
      bysize  = unidrnd(nobj);
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      [hh] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize);
    case 2;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      nobj     = size(XX1,2);
      [~, OFscales] = create_fignames(nobj,'obj');
      bycolor = unidrnd(nobj);
      bysize  = unidrnd(nobj);
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      [hh] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize);
    case 3;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      nobj     = size(XX1,2);
      bycolor = unidrnd(nobj);
      bysize  = unidrnd(nobj);
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      [hh] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize);
    case 4;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      
      nobj     = size(XX1,2);
      if isempty(bycolor);
        bycolor = unidrnd(nobj); 
      end
      bysize  = unidrnd(nobj);
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      [hh] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize);
    case 5;
      XX1      = varargin{1};
      nobj     = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      if isempty(bycolor)
        bycolor = unidrnd(nobj); 
      end
      if isempty(bysize)
        bysize = unidrnd(nobj); 
      end
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
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
        case {0,1};
          disp('Too few decision variables to display');
        case 2;
          % order of combinations of 2 variables for plotting
          ifigs = nchoosek(1:nobj,2);
          nsub = size(ifigs,1);
          ncols = ceil(sqrt(nsub));
          nrows = ceil(nsub/ncols);
          disp(' ');
          disp(['   Figure with nrows : ',num2str(nrows),' and ncolumns : ',num2str(ncols)]);
          disp(['   by Color : OF',num2str(bycolor)]);
          disp(['   by Size  : OF',num2str(bysize)]);          

          hh = figure(ifig);
          set(hh,'Color',[1.0 1.0 1.0]);
          ndata = size(XX1,1);
          [out_clr] = create_colors(6);
          [out_clr] = create_colors(ndata,out_clr);
          XX1 = sortrows(XX1,bycolor);
          xmin = min(XX1);
          xmax = max(XX1);
          nsizes = 6;
          [out_sizes, out_sizestick] = create_sizes(XX1,bysize,nsizes);
          jobj = 2;
          iobj = 1;
          for iclus = 1:nsizes;
            xclus = find(out_sizes==iclus,1);
            if ~isempty(xclus)
              plot(XX1(xclus,iobj),XX1(xclus,jobj),'Marker',out_mrk{1},...
                'Markersize',iclus,'MarkerEdgecolor','k',...
                'Markerfacecolor','No','lineStyle','No'); hold on;
              xleg_a = sprintf(strformat,out_sizestick(iclus));
              xleg_b = sprintf(strformat,out_sizestick(iclus+1));
              legend_sizes = cat(2,[xleg_a,'-',xleg_b],legend_sizes);
            end
          end
          set(gca,'Fontname',sel_fontname);
          set(gca,'Fontsize',sel_fontsize);
          set(gca,'Fontweight',sel_fontweight);
            
          % Set colormap ticks
          colormap(out_clr);
          xclrtick = (0:0.20:1.0);
          xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;

          % Create colorbar tick labels
          xclrlabel = {};
          for ii = 1:length(xclrtick);
            xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
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
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom',...
            'Position', [0.15, 1.05, 0]);

          % Set legend. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          %  something like...
          % text('String',['by Size : OF_',num2str(bysize)],'Position',<...>);
          leg = legend(legend_sizes,'Box','Off','Position',[0.90 0.275 0.09 0.13]);
          text('Parent', leg.DecorationContainer, ...
            'String', ['by Size : OF_',num2str(bysize)], ...
            'FontSize', 8.0,...
            'Fontweight','Bold',...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'VerticalAlignment', 'bottom',...
            'Position', [0.25, 1.05, 0]);
          
          for idata = 1:ndata
            plot(XX1(idata,iobj),XX1(idata,jobj),'Marker',out_mrk{1},...
              'Markersize',out_sizes(idata),'MarkerEdgecolor','k',...
              'Markerfacecolor',out_clr(idata,:),'lineStyle','No'); 
            hold on;
          end
          set(gca,'Fontname',sel_fontname);
          set(gca,'Fontsize',sel_fontsize);
          set(gca,'Fontweight',sel_fontweight);
          set(gca,'xlim', Rmin*round(Rmax*[xmin(iobj) xmax(iobj)]));
          set(gca,'xtick',Rmin*round(Rmax*[xmin(iobj) 0.5*(xmax(iobj)+xmin(iobj)) xmax(iobj)]));
          set(gca,'ylim', Rmin*round(Rmax*[xmin(jobj) xmax(jobj)]));
          set(gca,'ytick',Rmin*round(Rmax*[xmin(jobj) 0.5*(xmax(jobj)+xmin(jobj)) xmax(jobj)]));
          set(gca,'XScale',OFscales{iobj});
          set(gca,'YScale',OFscales{jobj});
          xlabel(OFnames{1});
          ylabel(OFnames{2});
          axis square;
          grid on;
        case {3, 4, 5, 6, 7};
          % order of combinations of 2 variables for plotting
          ifigs = nchoosek(1:nobj,2);
          nsub = size(ifigs,1);
          ncols = ceil(sqrt(nsub));
          nrows = ceil(nsub/ncols);
          nsizes = 6;
          disp(' ');
          disp(['   Figure with nrows : ',num2str(nrows),' and ncolumns : ',num2str(ncols)]);
          disp(['   by Color : OF',num2str(bycolor)]);
          disp(['   by Size  : OF',num2str(bysize)]);          
          make_it_tight = true;         
          % make_it_tight = false;
          %                    subtightplot(m, n, p,          gap,         marg_h,      marg_w,varargin)
%           subplot = @(m,n,p) subtightplot(m, n, p, [0.03 0.025], [0.0825 0.025], [0.035 0.095]);
          subplot = @(m,n,p) subtightplot(m, n, p, [0.055 0.065], [0.075 0.05], [0.065 0.125]);
          if ~make_it_tight,  
            clear subplot;  
          end
          hh = figure(ifig);
          set(hh,'Color',[1.0 1.0 1.0]);
          ndata = size(XX1,1);
          [out_clr] = create_colors(6);
          [out_clr] = create_colors(ndata,out_clr);
          XX1 = sortrows(XX1,bycolor);
          xmin = min(XX1);
          xmax = max(XX1);
          [out_sizes, out_sizestick] = create_sizes(XX1,bysize,nsizes);
          for isub = 1:nsub
            iobj = ifigs(isub,1);
            jobj = ifigs(isub,2);
            subplot(nrows,ncols,isub);
            if isub == nsub;
              for iclus = 1:nsizes;
                xclus = find(out_sizes==iclus,1);
                if ~isempty(xclus)
                  plot(XX1(xclus,iobj),XX1(xclus,jobj),'Marker',out_mrk{1},...
                  'Markersize',iclus,'MarkerEdgecolor','k',...
                  'Markerfacecolor','No','lineStyle','No'); hold on;
                  xleg_a = sprintf(strformat,out_sizestick(iclus));
                  xleg_b = sprintf(strformat,out_sizestick(iclus+1));
                  legend_sizes = cat(2,[xleg_a,'-',xleg_b],legend_sizes);
                end
              end
            end
            set(gca,'Fontname',sel_fontname);
            set(gca,'Fontsize',sel_fontsize);
            set(gca,'Fontweight',sel_fontweight);
            
            for idata = 1:ndata 
              plot(XX1(idata,iobj),XX1(idata,jobj),'Marker',out_mrk{1},...
                'Markersize',out_sizes(idata),'MarkerEdgecolor','k',...
                'Markerfacecolor',out_clr(idata,:),'lineStyle','No'); hold on;
            end
            set(gca,'Fontname',sel_fontname);
            set(gca,'Fontsize',sel_fontsize);
            set(gca,'Fontweight',sel_fontweight);
            set(gca,'xlim', Rmin*round(Rmax*[xmin(iobj) xmax(iobj)]));
            set(gca,'xtick',Rmin*round(Rmax*[xmin(iobj) 0.5*(xmax(iobj)+xmin(iobj)) xmax(iobj)]));
            set(gca,'ylim', Rmin*round(Rmax*[xmin(jobj) xmax(jobj)]));
            set(gca,'ytick',Rmin*round(Rmax*[xmin(jobj) 0.5*(xmax(jobj)+xmin(jobj)) xmax(jobj)]));
            set(gca,'XScale',OFscales{iobj});
            set(gca,'YScale',OFscales{jobj});
            xlabel(OFnames{iobj});
            ylabel(OFnames{jobj});
            grid on;
          end % isub
          % Set colormap
          colormap(out_clr);
          xclrtick = (0:0.20:1.0);
          xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;

          % Create colorbar ticks
          xclrlabel = {};
          for ii = 1:length(xclrtick);
            xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
          end % ii

          % Set colorbar
          hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'Location','manual','Position',[0.930 0.55 0.025 0.20]);

          % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          % something like...
          %   get(hc,'xlabel');
          %   set(get(hc,'xlabel'),'String',['by Color : OF_',num2str(bycol)]);
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
          % something like...
          %   text('String',['by Size : OF_',num2str(bysize)],'Position',<...>);
          leg = legend(legend_sizes,'Box','Off','Position',[0.90 0.275 0.09 0.13]);
          text('Parent', leg.DecorationContainer,'String', ['by Size : OF_',num2str(bysize)],...
               'FontSize', 8.0,'Fontweight','Bold','Position', [0.25, 1.05, 0], ...
               'Units', 'normalized','HorizontalAlignment', 'left', ...
               'VerticalAlignment', 'bottom');
        otherwise
          disp('Too many objectives to display');
      end % switch nobj
    hold off;
    grid on;
    drawnow;
  end % switch nargin
end % function