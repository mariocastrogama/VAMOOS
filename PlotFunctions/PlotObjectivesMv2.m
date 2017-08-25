function [hh,idx] = PlotObjectivesMv2(varargin)
% function [hh,idx] = PlotObjectivesM(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,OFtypes,ndecimals)
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
% [1] XX1      : matrix of size [ndata, nobjs]
% [2] OFnames  : list of names of the OF's
% [3] OFscales : Scales of each of the time series
% [4] bycolor  : number of the column to order by color
% [5] bysize   : number of the column to order by size
% [6] ncolor   : number of colors
% [7] nsizes   : number of sizes
% [8] out_clr  : a palette of colors to use for interpolation [ni x 3]
% [9] xbrush   : [xmin, xmax] two values which define the limits of the
%             variable bycolor, if the limits are violated no brushing is
%             performed
% [10] OFtypes   : define if function is for maximization or minimizqtion
% [11] ndecimals : number of decimals of each Objective (required for better display)
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
%   2015-11-16, modified some coloring issues
%   2016-05-31, correcting the number of colors that are plot to a limited number. 
%               It will increase speed of computation
%               by limiting the number of legends in the colorbar
%   2016-06-06, fixed a bug with the number of decimals for
%               display in legend and colorbar.
%   2016-09-01, allowed brushing of solutions in specified
%               range with variable xbrush
%
%  % Test
%   nsub     = 500;
%   XX1      = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1)];
%   out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
%   xbrush   = [0.50 1.75];
%  [h3,id3]  = PlotObjectivesMv2(XX1,[],[],3,4,[],[],out_clr,xbrush);
%
  tplot = tic;
%   [out_mrk] = create_markers(1);
  sel_fontname   = 'Arial';
  sel_fontsize   = 13;
  sel_fontweight ='bold';
  sel_fontsize_CLR = 14;
  legend_sizes   = {};
  sel_marker     = 'v'; %'o'; %;'v'; %'o'; '^'; % 'square' | 'diamond' | 'v' | '^' | '>' | '<' 
  % it indicates the numerical precision of results
%   ndecimals = [0 1 2 0 2];
  
  % ndecimals = [3 3 3 3 3];


  % select which monitor to use 1 or 2
  select_screen(1);

  switch nargin
    case 0;
      disp(' ');
      disp(' % Running test [500x5]');
      disp('   XX1 = [1.01+3.55*rand(500,1), 3.14+2.01*rand(500,1), rand(500,1), -1.85+2.03*rand(500,1),eps+rand(nsub,1)];');
      disp('   [hh,idx] = PlotObjectivesMv2(XX1);');
      disp(' ');
      nsub     = 100;
      XX1      = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1),eps+rand(nsub,1)];
      PlotObjectivesMv2(XX1);
    case 1;
      XX1      = varargin{1}; nobj = size(XX1,2);
      
      [OFnames, OFscales] = create_fignames(nobj,'obj');
      bycolor  = unidrnd(nobj);
      bysize   = unidrnd(nobj);
      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 2;
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};

      [~, OFscales] = create_fignames(nobj,'obj');
      bycolor  = unidrnd(nobj);
      bysize   = unidrnd(nobj);
      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 3;
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};

      bycolor  = unidrnd(nobj);
      bysize   = unidrnd(nobj);
      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 4;
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};

      bysize   = bycolor;
      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 5;
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};

      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 6;
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};

      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 7
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};

      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 8
      XX1      = varargin{1}; nobj = size(XX1,2);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
      
      xbrush   = [min(XX1(:,bycolor))-0.001, max(XX1(:,bycolor))+0.001];
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 9
      XX1      = varargin{1}; nobj = size(XX1);
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
      xbrush   = varargin{9};
      
      ndecimals = create_ndecimals(XX1);
      [OFtypes, ~] = create_fignames(nobj,'type');
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
      
    case 10
      XX1       = varargin{1}; nobj = size(XX1);
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = varargin{4};
      bysize    = varargin{5};
      ncolor    = varargin{6};
      nsizes    = varargin{7};
      out_clr   = varargin{8};
      xbrush    = varargin{9};
      ndecimals = varargin{10};
      [OFtypes, ~] = create_fignames(nobj,'type');
      
      [hh,idx] = PlotObjectivesMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals,OFtypes);
    case 11
      XX1       = varargin{1}; [ndata, nobj] = size(XX1);
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = varargin{4};
      bysize    = varargin{5};
      ncolor    = varargin{6};
      nsizes    = varargin{7};
      out_clr   = varargin{8};
      xbrush    = varargin{9};
      ndecimals = varargin{10}; 
      OFtypes   = varargin{11};
      
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
      if isempty(ncolor);
        ncolor = 6;
        disp(' setting >> ncolor = 10');
      end
      if isempty(nsizes);
        nsizes = 6;
        disp(' setting >> nsizes = 6');
      end
      if isempty(out_clr);
        out_clr = create_colors(ncolor);
      else
        if (size(out_clr,2)~=3) || (size(out_clr,1)==1);
          error(' Colorpalette must have at least 2 rows (colors) and always 3 columns [r, g, b]');
        end
      end
      if isempty(xbrush);
        xbrush = [min(XX1(:,bycolor))-1, max(XX1(:,bycolor))+1];
      end
      if isempty(OFtypes);
        [OFtypes, ~] = create_fignames(nobj,'type');
      end
      if isempty(ndecimals);
        ndecimals = create_ndecimals(XX1);
      end
      if length(ndecimals) < nobj
        error('vector of number of decimals is less than required');
      end
      
      xdir = {};
      for itype =1:length(OFtypes)
        if strcmp(OFtypes{itype},'max') ==1;
          xdir = cat(2,xdir,'reverse');
        else
          xdir = cat(2,xdir,'normal');
        end
      end
      
      % this can be set to plot different variables
      strformat = {};
      for iobj = 1:nobj;
        strformat = cat(2,strformat,['%4.',num2str(ndecimals(iobj)),'f']);
      end
      strformat
      
      
      
      % This is the switch for plotting
      switch nobj
        case {0,1,2,3};
          error('Too few objectives to display using color and size classification');
        case {4, 5, 6, 7, 8, 9, 10};
          % order of combinations of 2 variables for plotting
          ifigs = nchoosek(1:nobj,2);
          ifigs = sortrows(ifigs,2);
          
          % single plot
          % ifigs = ifigs(1,:); 
          
          % Only plot the figures which do not include the bycolor and
          % bysize columns
          for ii=size(ifigs,1):-1:1
            r1 = length(find(ifigs(ii,:)==bycolor));
            r2 = length(find(ifigs(ii,:)==bysize ));
            irej = r1+r2;
            if (irej > 0) % get rid of that figure
              ifigs(ii,:) = [];
            end
          end
          clear r1 r2 irej;

          nsub  = size(ifigs,1);
          ncols = 2;%ceil(sqrt(nsub));
          nrows = 2;%ceil(nsub/ncols);
          disp(' ');
          disp('   Plot of 2D trade-offs among objectives');
          disp(['   Figure with nrows : ',num2str(nrows),' and ncolumns : ',num2str(ncols)]);
          disp(['   by Color : OF',num2str(bycolor),' ',OFnames{bycolor}]);
          disp(['   by Size  : OF',num2str(bysize),' ',OFnames{bysize}]);
          if (xbrush(2) < min(XX1(:,bycolor)));
            disp('   Max value of xbrush < min(X(:,bycolor))');
            disp('   Reverting to min and max of bycolor');
            xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
          end
          if (xbrush(1) > max(XX1(:,bycolor)));
            disp('   min value of xbrush > max(X(:,bycolor))');
            disp('   Reverting to min and max of bycolor');
            xbrush = [min(XX1(:,bycolor)), max(XX1(:,bycolor))];
          end
          make_it_tight = true;
          % make_it_tight = false;
          %                    subtightplot(m, n, p,          gap,         [marg_bottom, marg_top],  [marg_left, marg right],varargin)
          %           subplot = @(m,n,p) subtightplot(m, n, p, [0.03 0.025], [0.0825 0.025], [0.035 0.095]);
          subplot = @(m,n,p) subtightplot(m, n, p, [0.0875 0.105], [0.075 0.05], [0.115 0.25]);
          if ~make_it_tight,
            clear subplot;
          end
          
          hh = figure;%(ifig);
          set(hh,'Color',[1.0 1.0 1.0]);
          set(hh,'Position',[63  1  1538  833]);%[70 10 900 810]);
          xmin = min(XX1);
          xmax = max(XX1);
          
          % Select the Best Optimal Point depending if OF is 
          % maximization or minimization
%           xBOP = zeros(1,nobj);
          xBOP = xmin;
%           for iobj =1:nobj
%             if strcmp(xtype{iobj},'max') == 1;
%               xBOP(iobj) = xmax(iobj);
%             else
%               xBOP(iobj) = xmin(iobj);
%             end
%           end

          % Selection of a valid range between two different values
          sel_brush = find((XX1(:,bycolor) > xbrush(1)) & (XX1(:,bycolor) < xbrush(2)) == 1);
          sel_shadow = find((XX1(:,bycolor) <= xbrush(1)) | (XX1(:,bycolor) >= xbrush(2)) == 1);
          
          nbrush   = size(sel_brush,1);
          XXbrush  = XX1(sel_brush,:);
          
          nshadow  = size(sel_shadow,1);
          XXshadow = XX1(sel_shadow,:);
                    
          % Estimate maximum and minimum of selection
          xbrushmin = min(XXbrush);
          xbrushmax = max(XXbrush);
          
          % Check if it found any feasible solution inside xbrush, otherwise exit
          if nbrush == 0;
            error(' The range of the selected variable does not find any feasible solution');
          end
          
          % find best Trade-off point of MOO
          norm_type  = 2;
          [euc]      = norm_forall(XXbrush,norm_type);
          [~,BTOpos] = min(euc);
          clear euc;
          
          [out_clr] = create_colors(ncolor,out_clr);
          %[ idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XX1,bysize,nsizes,bycolor,ncolor);
          [ idx, tick_size, tick_color, out_sizes, out_colors ] = create_sizes_v2(XXbrush,bysize,nsizes,bycolor,ncolor);
          
          % determine the right value for a MAX or MIN objective
          switch xdir{bysize}
            case 'reverse';
              ksize = -1;
              arrow_bysize = '\uparrow';
            case 'normal';
              ksize = +1;
              arrow_bysize = '\downarrow';
          end
            
          for isub = 1:nsub;
            x1 = ifigs(isub,1);
            y1 = ifigs(isub,2);
            % there is no Z value in 2D
            subplot(nrows,ncols,isub);
            
            % Draw the series once to get the proper legend with the
            % corresponding sizes. Plot with no color.
            % only do it on the last subplot
            if (isub == nsub);
              for isize = 1:nsizes;
                xclus = find(out_sizes == isize,1); % only find the first value with that marker size
                if ~isempty(xclus);
                  %plot(XX1(xclus,x1),XX1(xclus,y1),...
                  plot(XXbrush(xclus,x1),XXbrush(xclus,y1),...
                    'Marker',sel_marker,...  %'Marker',out_mrk{1},...
                    'lineStyle','No',...
                    'Markersize',(nsizes - isize + 1),...
                    'MarkerEdgecolor',[0.2 0.2 0.2],...
                    'Markerfacecolor','No'); 
                  hold on;
                  xleg_a = sprintf(strformat{bysize},ksize*tick_size(nsizes-isize+1));
                  xleg_b = sprintf(strformat{bysize},ksize*tick_size(nsizes-isize+2));
                  legend_sizes = cat(2,[xleg_a,' to ',xleg_b],legend_sizes);
                end % ~isempty(xclus);
              end % iclus;
              % Plot the  desired point BOP
              plot(round(xBOP(x1),ndecimals(x1)),round(xBOP(y1),ndecimals(y1)),...
                'p','Markersize',13,...
                'linewidth',1,...
                'MarkerEdgecolor','k',...
                'Markerfacecolor','y');
              hold on;
              legend_sizes = cat(2,legend_sizes,'Ideal Point');
              
%               % Plot the Best Trade-Off point 
%               plot(XXbrush(BTOpos,x1),XXbrush(BTOpos,y1),...
%                 'h','Markersize',12,...
%                 'linewidth',2,...
%                 'MarkerEdgecolor','k',...
%                 'Markerfacecolor','w');
%               hold on;
%               legend_sizes = cat(2,legend_sizes,'BTO');
              
              % plot the shadow data LAST of all
              if nshadow > 0;
                plot(XXshadow(:,x1),XXshadow(:,y1),...
                  'x','Markersize',2,...
                  'Markeredgecolor',[0.85 0.85 0.85]);
                hold on;
                legend_sizes = cat(2,legend_sizes,'Brushed data');
              end
            end % (isub == n);
            
            % plot the shadow data for all other sublplots
            if nshadow > 0;
              plot(XXshadow(:,x1),XXshadow(:,y1),...
                'x','Markersize',2,...
                'Markeredgecolor',[0.85 0.85 0.85]);
            end
            hold on;
            
            % The second time we do the classification is when we really
            % draw each of the data points. Plot with color.
            ntot = 0;
            ser  = 0;
            for isize = nsizes:-1:1;
              for icolor = ncolor:-1:1;
                ser = ser + 1;
                xrange = idx{isize,icolor};
                ndataclus = length(xrange);
                ntot = ntot + ndataclus;
                if ndataclus > 0;
                  %plot(XX1(xrange,x1),XX1(xrange,y1),...
                  plot(XXbrush(xrange,x1),XXbrush(xrange,y1),...
                    'Marker',sel_marker,...  %'Marker',out_mrk{1},...
                    'lineStyle','No',...
                    'Markersize',(nsizes - isize + 1),...
                    'Markerfacecolor',out_clr(icolor,:),...
                    'MarkerEdgecolor',out_clr(icolor,:));
%                     'MarkerEdgecolor',[0.2 0.2 0.2]);
                end % ndataclus > 0;
                hold on;
              end % icolor;
            end % isize;
            
            % Plot the  desired point            
            plot(round(xBOP(x1),ndecimals(x1)),round(xBOP(y1),ndecimals(y1)),...
              'p','Markersize',13,...
              'linewidth',1,...
              'MarkerEdgecolor','k',...
              'Markerfacecolor','y');
            
%             % Plot the Best Trade-Off point 
%             plot(XXbrush(BTOpos,x1),XXbrush(BTOpos,y1),...
%               'h','Markersize',12,...
%               'linewidth',2,...
%               'MarkerEdgecolor','k',...
%               'Markerfacecolor','w');
            
            % Set current axis properties
            set(gca,'Fontname',sel_fontname);
            set(gca,'Fontsize',sel_fontsize);
            set(gca,'Fontweight',sel_fontweight);

            % Set x limits, tick, labelrotation, scale
            ntx = 5;
            xvals = round(linspace(xmin(x1),xmax(x1),ntx),ndecimals(x1));
            set(gca,'xlim',  round([xmin(x1) xmax(x1)],ndecimals(x1)));
            set(gca,'xtick', xvals);
            xsa = {}; 
            
            switch xdir{x1}
              case 'reverse';
                kxlabel = -1;
                prefix_xlabel = '';
                sufix_xlabel = '\rightarrow';
              case 'normal';
                kxlabel = +1;
                prefix_xlabel = '\leftarrow';
                sufix_xlabel = '';
            end
            for is = 1:ntx; 
              xs = sprintf(strformat{x1},kxlabel*xvals(is)); 
              xsa = cat(2,xsa,xs);
            end;
            set(gca,'xticklabel',xsa);
            
            set(gca,'XScale',OFscales{x1});
            set(gca,'XDir',xdir{x1});
            xlab = [prefix_xlabel,' ',OFnames{x1},' ',sufix_xlabel];
            xlabel(xlab,'HorizontalAlignment','left','VerticalAlignment','cap') % 'baseline' 'top' 'cap' 'middle' 'bottom'
            
            % Set y limits, tick, scale
            nty = 5;
            yvals = round(linspace(xmin(y1),xmax(y1),nty),ndecimals(y1));
            set(gca,'ylim',  round([xmin(y1) xmax(y1)],ndecimals(y1)));
            set(gca,'ytick', yvals);
            ysa = {}; 
            switch xdir{y1}
              case 'reverse';
                kylabel = -1;
                prefix_ylabel = '';
                sufix_ylabel = '\rightarrow';
              case 'normal';
                kylabel = +1;
                prefix_ylabel = '\leftarrow';
                sufix_ylabel = '';
            end
            for is = 1:nty; 
              ys = sprintf(strformat{y1},kylabel*yvals(is)); 
              ysa = cat(2,ysa,ys); 
            end;
            set(gca,'yticklabel',ysa);

            set(gca,'YScale',OFscales{y1});
            set(gca,'YDir',xdir{y1});
            ylab = [prefix_ylabel,' ',OFnames{y1},' ',sufix_ylabel];
            ylabel(ylab,'HorizontalAlignment','left');
            axis square;
            grid on;
          end % isub

          % Set colormap ticks equal number as divisions
          % less than 10 colors is fine
          xclrtick = linspace(0,1.0,ncolor+1);
%           xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;
          
          
          % Create colorbar ticks labels
          xclrlabel = {};
          switch xdir{bycolor}
            case 'reverse'
              kclr = -1;
              arrow_bycolor = '\uparrow';
            case 'normal'
              kclr = +1;
              arrow_bycolor = '\downarrow';
          end
          for ii = 1:ncolor+1;
%             xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
            xclrlabel = cat(2, xclrlabel, sprintf(strformat{bycolor},kclr*tick_color(ii)));
          end % ii

          % Set colorbar
          hc = colorbar('Ticks',xclrtick,...
            'TickLabels',xclrlabel,...
            'Location','manual',...
            'Direction','reverse',...
            'Position',[0.85 0.55 0.025 0.30],...%[0.560 0.115 0.050 0.300],...
            'FontSize', sel_fontsize_CLR,...
            'Fontweight','Bold',...
            'colormap',out_clr);

          % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          % something like...
          %
          %      get(hc,'xlabel');
          %      set(get(hc,'xlabel'),'String',['by Color : ',OFnames{bycol}]);
          %
          text('Parent', hc.DecorationContainer, ...
            'String', ['by Color : ',arrow_bycolor,' ',OFnames{bycolor}],...
            'FontSize', sel_fontsize,...
            'Fontweight','Bold',...
            'Position', [-0.15, -0.05, 0], ... 
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
            'Position',[0.83 0.060 0.088 0.36]);
          text('Parent', leg.DecorationContainer,...
            'String', ['by Size : ',arrow_bysize,' ',OFnames{bysize}],...
            'FontSize', sel_fontsize,...
            'Fontweight','Bold',...
            'Position', [0.05, 1.025, 0],... 
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
        disp(['   Brushing selection : ',num2str(nbrush),' out of ',num2str(ndata),' samples (',sprintf('%2.1f',nbrush/ndata*100),'%)']);
        for ii = 1:length(idx);
          ki        = length(idx{ii});
          if (ki > 0);
            idx{ii} = sel_brush(1:ki);
            sel_brush(1:ki) = [];
          end
        end
      end
      disp(['   plotting time ',sprintf('%2.1f',toc(tplot)),' [s]']);
%       pause;
    otherwise
      error(' Too many input arguments');
  end % switch nargin
end % function
