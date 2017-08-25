function [hh] = PlotObjectives3D(islands,OFscales,OFnames,bycol,bysize)
% plots the scatters of trade-off between objective functions obtained 
% for every island as 3D surfaces
% It helps to visualize any hiden trade-off not visible in low dimensional
% Pareto fronts
%
  global nobj;
  global nislands;
  global npop;
  global ngen;
  global igen;
  
%   [out_dot] = create_attributes(nislands);
  [out_col] = create_colors(nislands);
  [out_mrk] = create_markers(nislands);
  sel_fontname = 'Arial';
  sel_fontsize = 9;
  sel_fontweight = 'bold';
  ifig = 4;
  IslandsNames = {};
  legend_sizes = {};
  if isstruct(islands)
    switch nobj
      case {0,1,2}
        disp('Too few objectives to display as 3D');
      case 3
        figure(ifig);
        clf;
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.ObjectiveFunctions],nobj, npop(r)));
          plot3(XX1(1,:), XX1(2,:), XX1(3,:),'Marker',out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
          set(gca,'Fontname',sel_fontname);
          set(gca,'Fontsize',sel_fontsize);
          set(gca,'Fontweight',sel_fontweight);
          set(gca,'XScale',OFscales{1});
          set(gca,'YScale',OFscales{2});
          set(gca,'ZScale',OFscales{3});
          xlabel(OFnames(1));
          ylabel(OFnames(2));
          zlabel(OFnames(3));
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end
        title(['Objective Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames);
      case {4,5,6,7}
        ifigs = nchoosek(1:nobj,3);
        n = length(ifigs);
        ncols = ceil(sqrt(n));
        nrows = ceil(n/ncols);
        figure(ifig);
        set(gcf,'Color',[1 1 1])
        clf;
        make_it_tight = true; % make_it_tight = false;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.035 0.025], [0.065 0.025], [0.025 0.025]);
        if ~make_it_tight; 
          clear subplot; 
        end
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.ObjectiveFunctions],nobj, npop(r)))';
%           size(XX1)
          xmin = min(XX1);
          xmax = max(XX1);
          isub = 0;
          for irow = 1:nrows
            for icol = 1:ncols
              isub = isub + 1;
              if isub <= n
                subplot(nrows,ncols,isub);
                x1 = ifigs(isub,1);
                y1 = ifigs(isub,2);
                z1 = ifigs(isub,3);
                hh = plot3(XX1(:,x1),XX1(:,y1),XX1(:,z1),'Marker',out_mrk{r},'Markersize',4,...
                'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
                set(gca,'Fontname',sel_fontname);
                set(gca,'Fontsize',sel_fontsize);
                set(gca,'Fontweight',sel_fontweight);
                set(gca,'xlim', 0.01*round(100*[xmin(x1) xmax(x1)]));
                set(gca,'xtick',0.01*round(100*[xmin(x1) 0.5*(xmax(x1)+xmin(x1)) xmax(x1)]));
                xlabel(OFnames(x1));
                set(gca,'ylim', 0.01*round(100*[xmin(y1) xmax(y1)]));
                set(gca,'ytick',0.01*round(100*[xmin(y1) 0.5*(xmax(y1)+xmin(y1)) xmax(y1)]));
                ylabel(OFnames(y1));
                set(gca,'zlim', 0.01*round(100*[xmin(z1) xmax(z1)]));
                set(gca,'ztick',0.01*round(100*[xmin(z1) 0.5*(xmax(z1)+xmin(z1)) xmax(z1)]));
                zlabel(OFnames(z1));
                axis square;
                grid on;
              end
            end % icol
          end % irow
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end % r
        set(gcf,'tag',['Objective Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames,'Position',[0.940 0.815 0.030 0.100]);
%         title(['Objective Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
%         hold off;
      otherwise
        disp('Too many objectives to display');
    end
    hold off;
    grid on;
    drawnow;
  else % Not an structure from my optimization algorithm it is only a matrix
     switch nobj
      case {0,1,2}
        disp('Too few decision variables to display');
      case 3
        figure(ifig);
        clf;
        XX1 = islands;
        xmin = min(XX1);
        xmax = max(XX1);
        hh = plot3(XX1(:,1), XX1(:,2), XX1(:,3),'Marker',out_mrk{1},'Markersize',4,'MarkerEdgecolor','k','Markerfacecolor',out_col(1,:),'lineStyle','No'); hold on;
        set(gca,'Fontname',sel_fontname);
        set(gca,'Fontsize',sel_fontsize);
        set(gca,'Fontweight',sel_fontweight);
        set(gca,'xlim', 0.01*round(100*[xmin(1) xmax(1)]));
        set(gca,'xtick',0.01*round(100*[xmin(1) 0.5*(xmax(1)+xmin(1)) xmax(1)]));
        set(gca,'ylim', 0.01*round(100*[xmin(2) xmax(2)]));
        set(gca,'ytick',0.01*round(100*[xmin(2) 0.5*(xmax(2)+xmin(2)) xmax(2)]));
        set(gca,'zlim', 0.01*round(100*[xmin(3) xmax(3)]));
        set(gca,'ztick',0.01*round(100*[xmin(3) 0.5*(xmax(3)+xmin(3)) xmax(3)]));
        set(gca,'XScale',OFscales{1});
        set(gca,'YScale',OFscales{2});
        set(gca,'ZScale',OFscales{3});
        axis square;
        view(30,30);
        xlabel(OFnames(1));
        ylabel(OFnames(2));
        zlabel(OFnames(3));
        IslandsNames = cat(2,IslandsNames,num2str(1));
        title('Objective Functions');%, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames);
      case {5,6,7}
        ifigs = nchoosek(1:nobj,3);
        n = length(ifigs);
        ncols = 5;%ceil(sqrt(n));
        nrows = 2;%ceil(n/ncols);
        make_it_tight = true; % make_it_tight = false;
        %                  subtightplot(m, n, p,           gap,        marg_h,       marg_w,varargin)
        subplot = @(m,n,p) subtightplot(m, n, p, [0.055 0.065], [0.09 0.05], [0.065 0.105]);
        if ~make_it_tight,  
          clear subplot;  
        end
        hh = figure(ifig);
        set(gcf,'Color',[1.0 1.0 1.0]);
        clf;
        XX1 = islands;
        ndata = size(XX1,1);
        [out_clr] = create_colors(6);
        [out_clr] = create_colors(ndata,out_clr);
        XX1 = sortrows(XX1,bycol);
        xmin = min(XX1);
        xmax = max(XX1);
        nsizes = 6;
        [out_sizes,out_sizestick] = create_sizes(XX1,bysize,nsizes);
        isub = 0;
        for irow = 1:nrows
          for icol = 1:ncols
            isub = isub + 1;
            if isub <= n;
              subplot(nrows,ncols,isub);
              x1 = ifigs(isub,1);
              y1 = ifigs(isub,2);
              z1 = ifigs(isub,3);
%               plot3(XX1(:,x1),XX1(:,y1),XX1(:,z1),'Marker',out_mrk{1},...
%                 'Markersize',3,'MarkerEdgecolor','k','Markerfacecolor',out_col(1,:),'lineStyle','No');
%               for idata = 1:ndata;
%                 plot3(XX1(idata,x1),XX1(idata,y1),XX1(idata,z1),'Marker',out_mrk{1},...
%                 'Markersize',out_sizes(idata),'MarkerEdgecolor','k','Markerfacecolor',out_clr(idata,:),'lineStyle','No'); hold on;
%               end
              
              % Draw the series once to get the proper legend with the
              % corresponding sizes
              for iclus = 1:nsizes;
                xclus = find(out_sizes==iclus,1);
                if ~isempty(xclus) && (isub == n);
                  plot3(XX1(xclus,x1),XX1(xclus,y1),XX1(xclus,z1),'Marker',out_mrk{1},'lineStyle','No',...
                    'Markersize',iclus,'MarkerEdgecolor','k','Markerfacecolor','No'); hold on;
                  xleg_a = sprintf('%3.2f',out_sizestick(iclus));
                  xleg_b = sprintf('%3.2f',out_sizestick(iclus+1));
                  legend_sizes = cat(2,[xleg_a,'-',xleg_b],legend_sizes);
                end
              end % iclus
              
              % The second time we do the classification is when we really
              % draw each of the data points
%               for iclus = 1:nsizes
%                 xclus = find(out_sizes==iclus);
%                 ndataclus = length(xclus);
%                 for idata = 2:ndataclus
%                   if ndataclus >0;
%                     plot3(XX1(xclus(idata),x1),XX1(xclus(idata),y1),XX1(xclus(idata),z1),'Marker',out_mrk{1},'lineStyle','No',...
%                     'Markersize',iclus,'MarkerEdgecolor','k','Markerfacecolor',out_clr(xclus(idata),:)); hold on;
%                   end
%                 end
%               end % iclus


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
%             axis square;
            grid on;
          end % icol
        end % irow
        
        % Set colormap
        colormap(out_clr);
        xclrtick = (0:0.20:1.0);
        xclrbar = xmin(bycol) + (xmax(bycol)-xmin(bycol))*xclrtick;
        
        % Create colorbar ticks
        xclrlabel = {};
        for ii=1:length(xclrtick);
          xclrlabel = cat(2,xclrlabel, sprintf('%4.3f',xclrbar(ii)));
        end
        
        % Set colorbar
        hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'Location','manual','Position',[0.930 0.675 0.025 0.20]);

        % Set colorbar title
        get(hc,'xlabel');
        set(get(hc,'xlabel'),'String',['OF_',num2str(bycol)]);
        
        % Set legend
        legend(legend_sizes,'Box','Off','Position',[0.90 0.35 0.09 0.13]);
        leg = findobj(gcf,'tag','legend');
        drawnow;
      otherwise
        disp('Too many objectives to display');
     end % switch nobj
    hold off;
    grid on;
    drawnow;
  end % isstruct
end