function [hh,idx] = PlotObjectives3DMv2(varargin)
% function [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush)
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
%  ncolor   : number of colors
%  nsizes   : number of sizes
%  out_clr  : a palette of colors to use for interpolation [ni x 3]
%  xbrush   : [xmin, xmax] two values which define the limits of the
%             variable bycolor, if the limits are violated no brushing is
%             performed
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
%   Last Update: 
%   2015-11-16, modified some coloring and bars
%   2016-05-31, correcting the number of colors that are plot to a limited number. 
%               It will increase speed of computation by limiting the number 
%               of legends in the colorbar
%   2016-06-06, fixed a bug with the number of decimals for display in legend and colorbar.
%   2016-09-01, allowed brushing of solutions in specified range with variable xbrush
%
%  % Test
%   nsub     = 500;
%   XX1      = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1)];
%   out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
%   xbrush   = [3.50 3.75];
%  [h3,idx] = PlotObjectives3DMv2(XX1,[],[],2,3,[],[],out_clr,xbrush);
%
  [out_mrk] = create_markers(1);
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  legend_sizes = {};
%   ifig = 4;
  
  
  % this can be set to plot different variables
  % it indicates the numerical precision of results
  ndecimals = 5; 
  Rmax = 10^ndecimals;
  Rmin = 1/Rmax;
  strformat = ['%4.',num2str(ndecimals),'f'];
  
  % select which monitor to use 1 or 2
  select_screen(2);
  
  switch nargin
    case 0;
      disp(' ');
      disp(' % Running test [250x5]');
      disp('   XX1 = [1.01+3.55*rand(500,1), 3.14+2.01*rand(500,1), rand(500,1), -1.85+2.03*rand(500,1)];');
      disp('   [hh,idx] = PlotObjectives3DMv2(XX1);');
      disp(' ');
      nsub = 250;
      XX1 = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1)];
      [hh,idx] = PlotObjectives3DMv2(XX1);
    case 1;
      XX1  = varargin{1};
      nobj = size(XX1,2);
      [OFnames, OFscales] = create_fignames(nobj,'obj');
      bycolor = unidrnd(nobj);
      bysize  = bycolor;
      while (bycolor == bysize)
        bysize = unidrnd(nobj); 
      end
      ncolor = 10;
      nsizes = 6;
      out_clr = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 2;
      XX1     = varargin{1};
      OFnames = varargin{2};
      nobj    = size(XX1,2);
      if isempty(OFnames);
        [OFnames, ~] = create_fignames(nobj,'obj');
      end
      if nobj > length(OFnames);
        error(' OFnames are incomplete');
      end      
      [~, OFscales] = create_fignames(nobj,'obj');
      bycolor = unidrnd(nobj);
      bysize  = bycolor;
      while (bycolor == bysize)
        bysize = unidrnd(nobj);
      end
      ncolor = 10;
      nsizes =  6;
      out_clr = [1, 0, 0; 1 1 0; 0 0.65 1];
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 3;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      nobj     = size(XX1,2);
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
      bycolor = unidrnd(nobj);
      bysize  = bycolor;
      while (bysize == bycolor)
        bysize = unidrnd(nobj);
      end
      ncolor = 10;
      nsizes =  6;
      out_clr = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 4;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      nobj     = size(XX1,2);
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
      if (nobj > bycolor) || (bycolor < 1)
        error(' bycol is not a feasible value');
      end
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      bysize = bycolor;
      while (bysize == bycolor)
        bysize = unidrnd(nobj);
      end
      ncolor   = 10;
      nsizes   =  6;
      out_clr  = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 5;
      XX1      = varargin{1};
      nobj     = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
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
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      if isempty(bysize);
        bysize = unidrnd(nobj);
      end
      while (bycolor == bysize);
        bysize = unidrnd(nobj);
      end      
      ncolor   = 10;
      nsizes   =  6;
      out_clr  = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 6;
      XX1      = varargin{1};
      nobj     = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
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
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      if isempty(bysize);
        bysize = unidrnd(nobj);
      end
      while (bycolor == bysize);
        bysize = unidrnd(nobj);
      end
      if isempty(ncolor);
        ncolor = 10;
        disp(' setting >> ncolor = 10');
      end
      nsizes   = 6;
      out_clr  = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 7
      XX1      = varargin{1};
      nobj     = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
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
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      if isempty(bysize);
        bysize = unidrnd(nobj);
      end
      while (bycolor == bysize);
        bysize = unidrnd(nobj);
      end
      if isempty(ncolor);
        ncolor = 10;
        disp(' setting >> ncolor = 10');
      end
      if isempty(nsizes);
        nsizes = 6;
        disp(' setting >> nsizes = 6');
      end      
      out_clr  = create_colors(ncolor);
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 8
      XX1      = varargin{1};
      [~,nobj] = size(XX1);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
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
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      if isempty(bysize);
        bysize = unidrnd(nobj);
      end
      while (bycolor == bysize);
        bysize = unidrnd(nobj);
      end
      if isempty(ncolor);
        ncolor = 10;
        disp(' setting >> ncolor = 10');
      end
      if isempty(nsizes);
        nsizes = 6;
        disp(' setting >> nsizes = 6');
      end
      if isempty(out_clr);
        out_clr = [1.00, 0.00, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      end
      if (size(out_clr,2)~=3) || (size(out_clr,1)==1);
        error(' Colorpalette must have at least 2 rows (colors) and always 3 columns [r, g, b]');
      end
      xbrush   = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      [hh,idx] = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush);
      
    case 9
      XX1      = varargin{1};
      [ndata,nobj] = size(XX1);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
      xbrush   = varargin{9};
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
      if isempty(bycolor);
        bycolor = unidrnd(nobj);
      end
      if isempty(bysize);
        bysize = unidrnd(nobj);
      end
%       while (bycolor == bysize);
%         bysize = unidrnd(nobj);
%       end
      if isempty(ncolor);
        ncolor = 10;
        disp(' setting >> ncolor = 10');
      end
      if isempty(nsizes);
        nsizes = 6;
        disp(' setting >> nsizes = 6');
      end
      if isempty(out_clr);
        out_clr = [1.00, 0.00, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      end
      if (size(out_clr,2)~=3) || (size(out_clr,1)==1);
        error(' Colorpalette must have at least 2 rows (colors) and always 3 columns [r, g, b]');
      end
      
      if isempty(xbrush);
        xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
      end
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
      
      % This is the switch for plotting
      switch nobj
        case {0,1,2};
          error('Too few objectives to display');
        case {3, 4, 5, 6, 7};
          % order of combinations of 3 variables for plotting
          ifigs = nchoosek(1:nobj,3);
          ifigs = sortrows(ifigs,3);
%           ifigs = ifigs(1,:); % single plot
          nsub = size(ifigs,1);
          ncols = ceil(sqrt(nsub));
          nrows = ceil(nsub/ncols);
          disp(' ');
          disp(['   Figure with nrows : ',num2str(nrows),' and ncolumns : ',num2str(ncols)]);
          disp(['   by Color : OF',num2str(bycolor),' ',OFnames{bycolor}]);
          disp(['   by Size  : OF',num2str(bysize),' ',OFnames{bysize}]);
          make_it_tight = true; 
          % make_it_tight = false;
          %                  subtightplot(m, n, p,           gap,       marg_h,       marg_w,  varargin)
          %
          subplot = @(m,n,p) subtightplot(m, n, p, [0.075 0.075], [0.075 0.05], [0.065 0.125]);
          if ~make_it_tight;
            clear subplot;
          end
          
          hh = figure;%(ifig);
          set(hh,'Color',[1.0 1.0 1.0]);
          xmin = min(XX1);
          xmax = max(XX1);

          % Selection of a valid range between two different values
          sel_brush = find((XX1(:,bycolor) > xbrush(1)) & (XX1(:,bycolor) < xbrush(2)) == 1);
          sel_shadow = find((XX1(:,bycolor) < xbrush(1)) | (XX1(:,bycolor) > xbrush(2)) == 1);
          
          nbrush = size(sel_brush,1);
          XXshadow = XX1(sel_shadow,:);
          
          nshadow = size(sel_shadow,1);
          XXbrush = XX1(sel_brush,:);
          
          % Estimate maximum and minimum of selection
          xbrushmin = min(XXbrush);
          xbrushmax = max(XXbrush);
          
          % Check if it finds any solution
          if nbrush == 0;
            error(' The range of the selected variable does not find any feasible solution');
          end
          
          % find best Trade-off point of MOO
          euc = zeros(ndata,1);
          for ii =1:ndata; 
            euc(ii,1)=norm(XX1(ii,:)); 
          end
          [~,BTOpos]=min(euc);
          clear euc;
          
          [out_clr] = create_colors(ncolor,out_clr);
          %[ idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
          [ idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XXbrush,bysize,nsizes,bycolor,ncolor);
          
          for isub = 1:nsub;
            ser = 0;
            x1 = ifigs(isub,1);
            y1 = ifigs(isub,2);
            z1 = ifigs(isub,3);
            subplot(nrows,ncols,isub);
            
            % Draw the series once to get the proper legend with the
            % corresponding sizes. Plot with no color.
            % only do it on the last subplot
            if (isub == nsub);
              for isize = 1:nsizes;
                xclus = find(out_sizes == isize,1); % only find the first value with that marker size
                if ~isempty(xclus);
                  %plot3(XX1(xclus,x1),XX1(xclus,y1),XX1(xclus,z1),...
                  plot3(XXbrush(xclus,x1),XXbrush(xclus,y1),XXbrush(xclus,z1),...
                    'Marker',out_mrk{1},...
                    'lineStyle','No',...
                    'Markersize',(nsizes - isize + 1),...
                    'MarkerEdgecolor',[0.2 0.2 0.2],...
                    'Markerfacecolor','No');
                  hold on;
                  xleg_a = sprintf(strformat,tick_size(nsizes-isize+1));
                  xleg_b = sprintf(strformat,tick_size(nsizes-isize+2));
                  legend_sizes = cat(2,[xleg_a,' - ',xleg_b],legend_sizes);
                end % ~isempty(xclus);
              end % iclus;
              % Plot the  desired point
              plot3(round(xmin(x1),ndecimals),round(xmin(y1),ndecimals),round(xmin(z1),ndecimals),...
                'p','Markersize',12,...
                'MarkerEdgecolor','k',...
                'Markerfacecolor','y');
              legend_sizes = cat(2,legend_sizes,'BOP');   
              
              % Plot the Best Trade-Off point 
%               plot3(XX1(BTOpos,x1),XX1(BTOpos,y1),XX1(BTOpos,z1),...
%                 's','Markersize',6,...
%                 'MarkerEdgecolor','k',...
%                 'linewidth',2,...
%                 'Markerfacecolor','w');
%               hold on;
%               legend_sizes = cat(2,legend_sizes,'BTO');

              % plot the shadow data LAST of all
              if nshadow > 0;
                plot3(XXshadow(:,x1),XXshadow(:,y1),XXshadow(:,z1),...
                  'x','Markersize',2,...
                  'Markeredgecolor',[0.85 0.85 0.85]);
                hold on;
                legend_sizes = cat(2,legend_sizes,'Shadow');
              end
            end % (isub == n);
            
            % plot the shadow data for all other sublplots
            if nshadow > 0;
              plot3(XXshadow(:,x1),XXshadow(:,y1),XXshadow(:,z1),...
                'x','Markersize',2,...
                'Markeredgecolor',[0.85 0.85 0.85]);
            end
            hold on;
            
            % The second time we do the classification is when we really
            % draw each of the data points. Plot with color.
            ntot = 0;
            for isize = nsizes:-1:1;
              for icolor = ncolor:-1:1;
                ser = ser + 1;
                xrange = idx{isize,icolor};
                ndataclus = length(xrange);
                ntot = ntot + ndataclus;
                if ndataclus > 0;
                  %plot3(XX1(xrange,x1),XX1(xrange,y1),XX1(xrange,z1),...
                  plot3(XXbrush(xrange,x1),XXbrush(xrange,y1),XXbrush(xrange,z1),...
                    'Marker',out_mrk{1},...
                    'lineStyle','No',...
                    'Markersize',(nsizes - isize + 1),...
                    'Markerfacecolor',out_clr(icolor,:),...
                    'MarkerEdgecolor',out_clr(icolor,:));                    
%                     'MarkerEdgecolor',[0.2 0.2 0.2]);
                end % ndataclus > 0;
                hold on;
              end % jclus;
            end % iclus;
            
            % Plot the  desired point
            plot3(round(xmin(x1),ndecimals),round(xmin(y1),ndecimals),round(xmin(z1),ndecimals),...
              'p','Markersize',12,...
              'MarkerEdgecolor','k',...
              'Markerfacecolor','y');

            % Plot the Best Trade-Off point
%             plot3(XX1(BTOpos,x1),XX1(BTOpos,y1),XX1(BTOpos,z1),...
%               's','Markersize',6,...
%               'MarkerEdgecolor','k',...
%               'linewidth',2,...
%               'Markerfacecolor','w');
            

            % Set current axis properties
            set(gca,'Fontname',sel_fontname);
            set(gca,'Fontsize',sel_fontsize);
            set(gca,'Fontweight',sel_fontweight);
            
            % Set x limits, tick, labelrotation, scale
            set(gca,'xlim', round([xmin(x1) xmax(x1)],ndecimals));
            set(gca,'xtick',Rmin*round(Rmax*[xmin(x1) 0.5*(xmax(x1)+xmin(x1)) xmax(x1)]));
            set(gca,'xticklabelrotation',-42);
            set(gca,'XScale',OFscales{x1});
            xlab = xlabel(OFnames{x1},'Rotation',5.25);
            set(xlab,'HorizontalAlignment','left');
            set(xlab,'VerticalAlignment','cap'); % 'baseline' 'top' 'cap' 'middle' 'bottom'
            
            % Set y limits, tick, scale
            set(gca,'ylim', round([xmin(y1) xmax(y1)],ndecimals));
            set(gca,'ytick',Rmin*round(Rmax*[xmin(y1) 0.5*(xmax(y1)+xmin(y1)) xmax(y1)]));
            set(gca,'YScale',OFscales{y1});
            ylab = ylabel(OFnames{y1},'Rotation',-40);
            set(ylab,'HorizontalAlignment','center');
            
            % Set z limits, tick, scale
            set(gca,'zlim', round([xmin(z1) xmax(z1)],ndecimals));
            set(gca,'ztick',Rmin*round(Rmax*[xmin(z1) 0.5*(xmax(z1)+xmin(z1)) xmax(z1)]));
            zlab = zlabel(OFnames{z1});
            set(zlab,'HorizontalAlignment','left');
            set(gca,'ZScale',OFscales{z1});
            axis square;
            grid on;
            view(-20,15);
%             hold on;
          end % isub;

          % Set colormap ticks equal number as divisions
          % less than 10 colors is fine
          xclrtick = (0:(1.0/ncolor):1.0);
          xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;
          
          % Create colorbar ticks labels
          xclrlabel = {};
          for ii = 1:length(xclrtick);
%             xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
            xclrlabel = cat(2, xclrlabel, sprintf(strformat,tick_color(ii)));
          end % ii

          % Set colorbar
          hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'Location','manual','Position',[0.930 0.55 0.025 0.20],'colormap',out_clr);

          % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          % something like...
          %
          %      get(hc,'xlabel');
          %      set(get(hc,'xlabel'),'String',['by Color : ',OFnames{bycol}]);
          %
          text('Parent', hc.DecorationContainer, ...
            'String', ['by Color : ',OFnames{bycolor}], ...
            'FontSize', 8.0,...
            'Fontweight','Bold',...
            'Position', [-0.15, 1.05, 0], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom');
          
          % Set legend. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          % something like...
          %
          %    text('String',['by Size : ',OFnames{bysize}],'Position',<...>);
          %
          leg = legend(legend_sizes,'Box','Off',...
            'Position',[0.90 0.275 0.09 0.13]);
          text('Parent', leg.DecorationContainer, ...
            'String', ['by Size : ',OFnames{bysize}], ...
            'FontSize', 8.0,...
            'Fontweight','Bold',...
            'Position', [0.25, 1.05, 0],...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'VerticalAlignment', 'bottom');
          drawnow;
        otherwise
          error(' Too many objectives to display, Max{nobj} = 7');
      end % switch nobj
      hold off;
      grid on;
      drawnow;
      
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
    otherwise
      error(' Too many input arguments');
  end % switch nargin
end % function