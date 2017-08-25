function [hh,idx] = PlotObjectives3DMv2(varargin)
% version 2
% function [hh] = PlotObjectives3DM(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals)
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
% [1] XX1      : matrix of size [ndata, nobjs]
% [2] OFnames  : list of names of the OF's
% [3] OFscales : Scales of each of the time series
% [4] bycolor  : number of the column to order by color
% [5] bysize   : number of the column to order by size 
% [6] ncolor   : number of colors
% [7] nsizes   : number of sizes
% [8] out_clr  : a palette of colors to use for interpolation [ni x 3]
% [9] xbrush   : [xmin, xmax] two values which define the limits of the variable bycolor, 
%                if the limits are violated NO brushing is performed
% [10] ndecimals : number of decimals of each Objective (required for better display)
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
%   2016-10-01, fixed a bug with the decimals for display for each objective and linked to legend and colorbar
%   2016-11-01, allow to create a GIF, very useful for presentations
%   2016-12-01, fixed the squares behind data to improve contrast of
%               brushed data. However, squares so not move when generating a GIF 
%
%  % Test
%   nsub     = 500;
%   XX1      = [1.01+3.55*rand(nsub,1), 3.14+2.01*rand(nsub,1), rand(nsub,1), -1.85+2.03*rand(nsub,1)];
%   out_clr  = [1, 0, 0; 1 1 0; 0 0.65 1];
%   xbrush   = [3.50 3.75];
%  [h3,idx] = PlotObjectives3DMv2(XX1,[],[],2,3,[],[],out_clr,xbrush);
%
  tplot = tic;
  [out_mrk] = create_markers(1);
  sel_fontname = 'Arial';
  sel_fontsize     = 13;
  sel_fontsize_CLR = 12;
  sel_fontsize_LEG = 12;
  sel_fontweight   = 'bold';
  sel_FaceAlpha    = 0.5;
  sel_brushcolor   = [0.5 0.5 0.5];
  legend_sizes = {};
  sel_marker     = 'v'; %'o'; %'o'; '^'; % 'square' | 'diamond' | 'v' | '^' | '>' | '<' 
  
  % number of decimals of each Objective Function
%   ndecimals = [4, 0, 2, 2, 2];  
  
  % Animation of visualization
  % GIF filename
  filename = 'WDN_ns_09_case_10_6.gif';  
  % filename = 'WATER_3D_RGB_v5.gif';
  % filename = 'GAA_3D_RGB_v5.gif';
  % rotation angles of GIF
  rot_angles = 0;%[0:20:359]; % 0:1; %359; %359;
  
  % select which monitor to use 1 or 2
  select_screen(1);
  
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
      bycolor   = unidrnd(nobj);
      bysize    = bycolor;
      ncolor    = 6;
      nsizes    = 6;
      out_clr   = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 2;
      XX1       = varargin{1};
      nobj      = size(XX1,2);
      OFnames   = varargin{2};
      [~, OFscales] = create_fignames(nobj,'obj');
      bycolor   = unidrnd(nobj);
      bysize    = bycolor;
      ncolor    = 6;
      nsizes    = 6;
      out_clr   = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 3;
      XX1       = varargin{1};
      nobj      = size(XX1,2);
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = unidrnd(nobj);
      bysize    = bycolor;
      ncolor    = 6;
      nsizes    = 6;
      out_clr   = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 4;
      XX1       = varargin{1};
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = varargin{4};
      bysize    = bycolor;
      ncolor    = 6;
      nsizes    = 6;
      out_clr   = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 5;
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = 6;
      nsizes   = 6;
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 6;
      XX1       = varargin{1};
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = varargin{4};
      bysize    = varargin{5};
      ncolor    = varargin{6};
      nsizes    = 6;
      out_clr   = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 7
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = [1.00, 0.20, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 8
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
      xbrush    = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
      
    case 9
      XX1      = varargin{1};
      OFnames  = varargin{2};
      OFscales = varargin{3};
      bycolor  = varargin{4};
      bysize   = varargin{5};
      ncolor   = varargin{6};
      nsizes   = varargin{7};
      out_clr  = varargin{8};
      xbrush   = varargin{9};
      ndecimals = create_ndecimals(XX1);
      [hh,idx]  = PlotObjectives3DMv2(XX1,OFnames,OFscales,bycolor,bysize,ncolor,nsizes,out_clr,xbrush,ndecimals);
    case 10
      XX1      = varargin{1};
      [ndata,nobj] = size(XX1);
      OFnames   = varargin{2};
      OFscales  = varargin{3};
      bycolor   = varargin{4};
      bysize    = varargin{5};
      ncolor    = varargin{6};
      nsizes    = varargin{7};
      out_clr   = varargin{8};
      xbrush    = varargin{9};
      ndecimals = varargin{10};
            
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
        out_clr = [1.00, 0.00, 0.00; 1.00 1.00 0.00; 0.00 0.65 1.00];
      end
      if (size(out_clr,2)~=3) || (size(out_clr,1)==1);
        error(' Colorpalette must have at least 2 rows (colors) and always 3 columns [r, g, b]');
      end
      
      if isempty(xbrush);
        xbrush = [min(XX1(:,bycolor))-0.0001, max(XX1(:,bycolor))+0.0001];
      end
      
      if length(ndecimals) < nobj
        error('vector of number of decimals is less than required');
      end
      if isempty(ndecimals);
        ndecimals = create_ndecimals(XX1);
      end
      % this can be set to plot different variables
      % it indicates the numerical precision of results
      strformat = {};
      for iobj = 1:nobj;
        strformat = cat(2,strformat,['%4.',num2str(ndecimals(iobj)),'f']);
      end
      strformat
      
      % This is the switch for plotting
      switch nobj
        case {0,1,2};
          error('Too few objectives to display');
        case {3, 4, 5, 6, 7, 10};
          % order of combinations of 3 variables for plotting
          ifigs = nchoosek(1:nobj,3);
          ifigs = sortrows(ifigs,3);
          
          % single plot
          % ifigs = ifigs(1,:); 

          % Only plot the figures which do not include the bycolor and
          % bysize columns
          for ii = size(ifigs,1):-1:1
            r1 = length(find(ifigs(ii,:)==bycolor));
            r2 = length(find(ifigs(ii,:)==bysize ));
            irej = r1 + r2;
            if (irej > 0) % get rid of that figure
              ifigs(ii,:) = [];
            end
          end
          clear r1 r2 irej;
          
          nsub  = size(ifigs,1);
          ncols = ceil(sqrt(nsub));
          nrows = ceil(nsub/ncols);
          disp(' ');
          disp('   Plot of 3D trade-offs among objectives');
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
          %                    subtightplot(m, n, p,  [gap_vert, gap_horz],  [marg_bottom, marg_top], [marg_left, marg right],varargin)
          %
          subplot = @(m,n,p) subtightplot(m, n, p, [0.05, 0.075], [0.050, 0.025], [0.075, 0.25]);
          if ~make_it_tight;
            clear subplot;
          end
          
          hh = figure;%(ifig);
          set(hh,'Color',[1.0 1.0 1.0]);
          set(hh,'Position',[63 1 1538 833]);%[70 10 1400 810]);
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
          norm_type = 2;
          [euc]      = norm_forall(XX1,norm_type);
          [~,BTOpos] = min(euc);
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
            % corresponding sizes. Plot with no color
            % only do it on the last subplot
            if (isub == nsub);
              for isize = 1:nsizes;
                xclus = find(out_sizes == isize,1); % only find the first value with that marker size
                if ~isempty(xclus);
                  %plot3(XX1(xclus,x1),XX1(xclus,y1),XX1(xclus,z1),...
                  plot3(XXbrush(xclus,x1),XXbrush(xclus,y1),XXbrush(xclus,z1),...
                    'Marker',sel_marker,... 'Marker',out_mrk{1},...
                    'lineStyle','No',...
                    'Markersize',(nsizes - isize + 1),...
                    'MarkerEdgecolor',[0.2 0.2 0.2],...
                    'Markerfacecolor','No');
                  hold on;
                  xleg_a = sprintf(strformat{bysize},tick_size(nsizes-isize+1));
                  xleg_b = sprintf(strformat{bysize},tick_size(nsizes-isize+2));
                  legend_sizes = cat(2,[xleg_a,' to ',xleg_b],legend_sizes);
                end % ~isempty(xclus);
              end % iclus;
              % Plot the  desired point
              plot3(xmin(x1),xmin(y1),xmin(z1),...
                'p','Markersize',20,...
                'MarkerEdgecolor','k',...
                'Markerfacecolor','y');
              legend_sizes = cat(2,legend_sizes,'Ideal Point');   
              
              % Plot the Best Trade-Off point 
%               plot3(XX1(BTOpos,x1),XX1(BTOpos,y1),XX1(BTOpos,z1),...
%                 's','Markersize',6,...
%                 'MarkerEdgecolor','k',...
%                 'linewidth',2,...
%                 'Markerfacecolor','w');
%               hold on;
%               legend_sizes = cat(2,legend_sizes,'BTO');

              % plot the brushed data LAST of all
              if nshadow > 0;
                plot3(XXshadow(:,x1),XXshadow(:,y1),XXshadow(:,z1),...
                  '.','Markersize',2,...
                  'Markeredgecolor',sel_brushcolor);
                hold on;
                legend_sizes = cat(2,legend_sizes,'Brushed data');
              end
            end % (isub == n);
            
            % plot the shadow data for all other sublplots
            if nshadow > 0;
              plot3(XXshadow(:,x1),XXshadow(:,y1),XXshadow(:,z1),...
                '.','Markersize',2,...
                'Markeredgecolor',sel_brushcolor);
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
                if (ndataclus > 0);
                  %plot3(XX1(xrange,x1),XX1(xrange,y1),XX1(xrange,z1),...
                  plot3(XXbrush(xrange,x1),XXbrush(xrange,y1),XXbrush(xrange,z1),...
                    'Marker',sel_marker,...'Marker',out_mrk{1},...
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
            plot3(xmin(x1),xmin(y1),xmin(z1),...
              'p','Markersize',20,...
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
            nlab = 3;
            xta = [];
            xsa = {};
            for is = 1:nlab; 
              xt = xmin(x1)+(is-1)/(nlab-1)*(xmax(x1)-xmin(x1));
              xs = sprintf(strformat{x1},round(xt,ndecimals(x1)));
%               if (is == 3); xsa{nlab} = cell2mat([xs,' ',OFnames{x1}]); else
              xsa = cat(2,xsa,xs);
              xta = cat(2,xta,xt);
%               end
            end
            set(gca,'xlim', [xmin(x1) xmax(x1)]);
            set(gca,'xtick',xta);
            set(gca,'xticklabel',xsa);      % set(gca,'xticklabelrotation',-42);
            set(gca,'XScale',OFscales{x1}); % set(gca,'xticklabel',xs);

            xlab = xlabel(OFnames{x1});%,'Rotation',5.25);
%             set(xlab,'HorizontalAlignment','left');
%             set(xlab,'VerticalAlignment','cap'); % 'baseline' 'top' 'cap' 'middle' 'bottom'
            
            % Set y limits, tick, scale
            ysa = {}; 
            yta = [];
            for is = 1:nlab; 
              yt = xmin(y1)+(is-1)/(nlab - 1)*(xmax(y1)-xmin(y1));
              ys = sprintf(strformat{y1},round(yt,ndecimals(y1)));
%               if (is == 3); ysa{nlab} = cell2mat([ OFnames{y1},' ',ys]); else
              ysa = cat(2,ysa,ys);
              yta = cat(2,yta,yt);
%               end
            end;
            set(gca,'ylim', [xmin(y1) xmax(y1)]);
            set(gca,'ytick',yta);
            set(gca,'yticklabel',ysa);
            set(gca,'YScale',OFscales{y1}); % set(gca,'yticklabel',ys);
            ylab = ylabel(OFnames{y1}); %,'Rotation',-40);
%             set(ylab,'HorizontalAlignment','center');
            
            % Set z limits, tick, scale
            zsa = {}; 
            zta = [];
            for is = 1:nlab; 
              zt = xmin(z1)+(is-1)/(nlab - 1)*(xmax(z1)-xmin(z1));
              zs = sprintf(strformat{z1},round(zt,ndecimals(z1))); 
              zsa = cat(2,zsa,zs);
              zta = cat(2,zta,zt);
            end
            set(gca,'zlim', [xmin(z1) xmax(z1)]);
            set(gca,'ztick',zta);
            set(gca,'zticklabel',zsa);
            zlab = zlabel(OFnames{z1});
            set(zlab,'HorizontalAlignment','left');
            set(gca,'ZScale',OFscales{z1});
            axis square;
%             set(gca,'Position',[0.01 0.10 0.80 0.85]);
            
            grid on;
            view(330,15);
            
            % plot the squared faced behind data for visualization of brushed data
            face_color = 0.975*ones(1,3);
            edge_color = (2/3)*face_color;
            frdlt = 80;
            xdlt = min(diff(xta))/frdlt; % xdlt = xdlt(1)/frdlt;
            ydlt = min(diff(yta))/frdlt; % ydlt = ydlt(1)/frdlt;
            zdlt = min(diff(zta))/frdlt; % zdlt = zdlt(1)/frdlt;
            
            for isq = 1:nlab-1;
              for jsq = 1:nlab-1;
                xface1 = [xta(isq)+xdlt, xta(isq+1)-xdlt, xta(isq+1)-xdlt, xta(isq)+xdlt,   xta(isq)+xdlt];
                yface1 = [yta(jsq)+ydlt, yta(jsq)+ydlt,   yta(jsq+1)-ydlt, yta(jsq+1)-ydlt, yta(jsq)+ydlt];
                patch(xface1,yface1,zta(1)*ones(1,5),face_color,'EdgeColor',edge_color,'FaceAlpha',sel_FaceAlpha);
                hold on;
                xface2 = [xta(isq)+xdlt, xta(isq+1)-xdlt, xta(isq+1)-xdlt, xta(isq)+xdlt,   xta(isq)+xdlt];
                zface2 = [zta(jsq)+zdlt, zta(jsq)+zdlt,   zta(jsq+1)-zdlt, zta(jsq+1)-zdlt, zta(jsq)+zdlt];
                patch(xface2,yta(end)*ones(1,5),zface2,face_color,'EdgeColor',edge_color,'FaceAlpha',sel_FaceAlpha);
                hold on;
                yface3 = [yta(isq)+ydlt, yta(isq+1)-ydlt, yta(isq+1)-ydlt, yta(isq)+ydlt,   yta(isq)+ydlt];
                zface3 = [zta(jsq)+zdlt, zta(jsq)+zdlt,   zta(jsq+1)-zdlt, zta(jsq+1)-zdlt, zta(jsq)+zdlt];
                patch(xta(end)*ones(1,5),yface3,zface3,face_color,'EdgeColor',edge_color,'FaceAlpha',sel_FaceAlpha);
                hold on;
              end
            end
            clear isq jsq
          end % isub;

          % Set colormap ticks equal number as divisions
          % less than 10 colors is fine
          xclrtick = (0:(1.0/ncolor):1.0);
          xclrbar = xmin(bycolor) + (xmax(bycolor)-xmin(bycolor))*xclrtick;
          
          % Create colorbar ticks labels
          xclrlabel = {};
          for ii = 1:length(xclrtick);
%             xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
            xclrlabel = cat(2, xclrlabel, sprintf(strformat{bycolor},tick_color(ii)));
          end % ii

          % Set colorbar
          hc = colorbar('Ticks',xclrtick,...
            'TickLabels',xclrlabel,...
            'Location','manual',...
            'Direction','reverse',...
            'Position',[0.85 0.55 0.025 0.30],...
            'FontSize', sel_fontsize_CLR,...
            'Fontweight','Bold',...
            'colormap',out_clr);
%           [0.712 0.150 0.019 0.825] [0.711904761904762 0.15 0.0190476190476191 0.825]
          % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
          % this changed recently between versions you may just need to use
          % something like...
          %
          %      get(hc,'xlabel');
          %      set(get(hc,'xlabel'),'String',['by Color : ',OFnames{bycol}]);
          %
          text('Parent', hc.DecorationContainer, ...
            'String', ['by Color : ',OFnames{bycolor}], ...
            'FontSize', sel_fontsize_CLR,...
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
            'Position',[0.85 0.025 0.09 0.33],...%[0.705 0.50 0.09 0.33],...
            'FontSize', sel_fontsize_LEG,...
            'Fontweight','Bold');
          text('Parent', leg.DecorationContainer, ...
            'String', ['by Size : ',OFnames{bysize}], ...
            'FontSize', sel_fontsize_LEG,...
            'Fontweight','Bold',...
            'Position', [0.25, 1.05, 0],...
            'Units', 'normalized',...
            'HorizontalAlignment', 'left',...
            'VerticalAlignment', 'bottom');
          
         
          time_delay = 0.01;
          FlagChangeLabels = 0;
          FlagStore = 0;
          FlagStoreNewGIF  = 0;
          for irot = rot_angles;% 0:1;%359;%359;
            theta = mod(315+irot,360);
            for isub = 1:nsub;
              subplot(nrows,ncols,isub);
              view(theta,30);
%               title(['\theta = ',num2str(theta)],'FontSize', sel_fontsize,'Fontweight','Bold');
              x = get(gca,'xticklabel');
              y = get(gca,'yticklabel');
              theta_cuad = mod(theta,180);
              % Change last element on X axis
%               if and(((0 <= theta_cuad) && (theta_cuad < 90)),FlagChangeLabels)
%                 x{3} = [x{3}(length(x{1})+2:end),' ',x{3}(1:length(x{1}))];
%                 set(gca,'xticklabel',x);
%                 y{3} = [y{3}(end-length(y{1})+1:end),' ', y{3}(1:end-length(y{1}))];
%                 set(gca,'yticklabel',y); 
%                 FlagChangeLabels = 0;
%               end
%               if and(((90 <= theta_cuad) && ( theta_cuad < 180)),~FlagChangeLabels)
%                 x{3} = [x{3}(end-length(x{1})+1:end),' ', x{3}(1:end-length(x{1}))];
%                 set(gca,'xticklabel',x);
%                 y{3} = [y{3}(length(y{1})+1:end),' ',y{3}(1:length(y{1}))];
%                 set(gca,'yticklabel',y);
%                 FlagChangeLabels = 1;
%               end
            end
            drawnow;
            % Store each frame and build a gif
            if (FlagStore == 1)
              frame = getframe(1);
              im    = frame2im(frame);
              [imind, cm] = rgb2ind(im,256);
              if (FlagStoreNewGIF == 1);              
                imwrite(imind,cm,filename,'gif','DelayTime',time_delay,'Loopcount',inf);
                FlagStoreNewGIF = 0;
              else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_delay);
              end
            end
          end % irot
        otherwise
          error(' Too many objectives to display, Max{nobj} = 7');
      end % switch nobj
      hold off;
      grid on;
%       drawnow;
      
      % extract the id's of the samples which were brushed IN
      if nbrush < size(XX1,1);
        disp(['   Brushing selection : ',num2str(nbrush),' out of ',num2str(ndata),' samples (',sprintf('%2.1f',nbrush/ndata*100),'%)']);
        for ii = 1:length(idx);
          ki        = length(idx{ii});
          if (ki > 0)
            idx{ii} = sel_brush(1:ki);
            sel_brush(1:ki) = [];
          end
        end
      end
      disp(['   plotting time ',sprintf('%2.1f',toc(tplot)),' [s]']);
    otherwise
      error(' Too many input arguments');
  end % switch nargin
end % function
