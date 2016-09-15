function [hh]=PlotDecision(islands,DVnames,DVscales)
% plots the scatters of trade-off between decision variables obtained 
% for every island
  global nislands
  global VarMin
  global VarMax
  global npop
  global igen
  global ngen
  
  if ~exist('nislands','var');
    nislands = 1;
  end
  
  [out_dot] = create_attributes(nislands);
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  ifig = 1;
  IslandsNames = {};
  hh = [];
  if isstruct(islands)
    nvar = size(islands(1).pop(1).DecisionVariables,2);
    [out_col] = create_colors(nislands);
    [out_mrk] = create_markers(nislands);
    switch nvar
      case {0,1}
        disp('Too few decision variables to display');
      case 2
        hh=figure(ifig);
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.DecisionVariables],nvar,npop(r)))';
          % 'Marker',
          plot(XX1(:,1), XX1(:,2),out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:)); hold on;%out_dot{r}); hold on;
          set(gca,'XScale',DVscales{1});
          set(gca,'YScale',DVscales{2});
          xlim([VarMin(1) VarMax(1)]);
          ylim([VarMin(2) VarMax(2)]);
          xlabel(DVnames(1));
          ylabel(DVnames(2));
          axis square;
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end
        title(['Decision Variables, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames);
        hold off;
        grid on;
        drawnow;
      case 3
        hh=figure(ifig);
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.DecisionVariables],nvar,npop(r)))';
          plot3(XX1(:,1), XX1(:,2),XX1(:,3),out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:)); hold on;
          set(gca,'XScale',DVscales{1});
          set(gca,'YScale',DVscales{2});
          set(gca,'ZScale',DVscales{3});
          xlim([VarMin(1) VarMax(1)]);
          ylim([VarMin(2) VarMax(2)]);
          zlim([VarMin(3) VarMax(3)]);
          xlabel(DVnames(1));
          ylabel(DVnames(2));
          zlabel(DVnames(3));
          axis square;
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end
        title(['Decision Variables, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames);
        hold off;
        grid on;
        drawnow;
      case {4, 5, 6, 7}
        hh=figure(ifig);
        clf;
        make_it_tight = true; % make_it_tight = false;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.025], [0.0825 0.075], [0.0625 0.075]);
        if ~make_it_tight; 
          clear subplot; 
        end
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.DecisionVariables],nvar, npop(r)))';
          xmin = min(XX1);
          xmax = max(XX1);
          isub = 0;
          for jvar = 1:nvar
            for ivar = 1:nvar
              isub = isub + 1;
              if (ivar ~= jvar)
                subplot(nvar,nvar,isub);
                hh = plot(XX1(:,ivar),XX1(:,jvar),'Marker',out_mrk{r},'Markersize',5,...
                'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
                set(gca,'Fontname',sel_fontname);
                set(gca,'Fontsize',sel_fontsize);
                set(gca,'Fontweight',sel_fontweight);
                set(gca,'xlim', 0.01*round(100*[xmin(ivar) xmax(ivar)]));
                set(gca,'xtick',0.01*round(100*[xmin(ivar) 0.5*(xmax(ivar)+xmin(ivar)) xmax(ivar)]));
                if isub<=(nvar*(nvar-1))
                  set(gca,'xticklabel',{});
                else
                  xlabel(DVnames(ivar));
                end
                set(gca,'ylim', 0.01*round(100*[xmin(jvar) xmax(jvar)]));
                set(gca,'ytick',0.01*round(100*[xmin(jvar) 0.5*(xmax(jvar)+xmin(jvar)) xmax(jvar)]));
                if mod(isub,nvar)~=1;  
                  set(gca,'yticklabel',{});
                else
                  ylabel(DVnames(jvar));
                end
                grid on;
              end
            end % ivar
          end % jvar
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end % r
        set(gcf,'tag',['Decision Variables, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames,'Position',[0.940 0.815 0.030 0.100]);
      otherwise
          disp('NO PLOT: Too many variables to display');
    end % switch
  else % Non structure just a MATRIX
    nvar = size(islands,2);
    [out_col] = create_colors(1);
    [out_mrk] = create_markers(1);
    switch nvar
      case 2
        hh = figure(ifig);
        clf;
        XX1 = islands;
        xmin = min(XX1);
        xmax = max(XX1);
        plot(XX1(:,1), XX1(:,2),out_mrk{1},'MarkerEdgecolor','k','Markerfacecolor',out_col(1,:),'Markersize',5); hold on;
        set(gca,'Fontname','Arial');
        set(gca,'Fontweight','Bold');
        set(gca,'Fontsize',14);
        set(gca,'XScale',DVscales{1});
        set(gca,'YScale',DVscales{2});
        if abs(xmin(1) - xmax(1)) > 0.01;
          set(gca,'xlim', 0.01*round(100*[xmin(1) xmax(1)]));
          set(gca,'xtick',0.01*round(100*[xmin(1) 0.5*(xmax(1)+xmin(1)) xmax(1)]));
        end
        xlabel(DVnames(1));
        if abs(xmin(2) - xmax(2)) > 0.01;
          set(gca,'ylim', 0.01*round(100*[xmin(2) xmax(2)]));
          set(gca,'ytick',0.01*round(100*[xmin(2) 0.5*(xmax(2)+xmin(2)) xmax(2)]));
        end
        ylabel(DVnames(2));
        axis square;
        grid on;
      case 3
        hh=figure(ifig);
        clf;
        XX1 = islands;
        xmin = min(XX1);
        xmax = max(XX1);
        plot3(XX1(:,1), XX1(:,2),XX1(:,3),'o','MarkerEdgecolor','k','Markerfacecolor',out_col(1,:),'Markersize',5); hold on;
        set(gca,'Fontname','Arial');
        set(gca,'Fontweight','Bold');
        set(gca,'Fontsize',14);
        set(gca,'XScale',DVscales{1});
        set(gca,'YScale',DVscales{2});
        set(gca,'ZScale',DVscales{3});
        xlim([xmin(1) xmax(1)]);
        ylim([xmin(2) xmax(2)]);
        zlim([xmin(3) xmax(3)]);
        if abs(xmin(1) - xmax(1)) > 0.01;
          set(gca,'xlim', 0.01*round(100*[xmin(1) xmax(1)]));
          set(gca,'xtick',0.01*round(100*[xmin(1) 0.5*(xmax(1)+xmin(1)) xmax(1)]));
        end
        xlabel(DVnames(1));
        if abs(xmin(2) - xmax(2)) > 0.01;
          set(gca,'ylim', 0.01*round(100*[xmin(2) xmax(2)]));
          set(gca,'ytick',0.01*round(100*[xmin(2) 0.5*(xmax(2)+xmin(2)) xmax(2)]));
        end
        ylabel(DVnames(2));
        if abs(xmin(3) - xmax(3)) > 0.011;
          set(gca,'zlim', 0.01*round(100*[xmin(3) xmax(3)]));
          set(gca,'ztick',0.01*round(100*[xmin(3) 0.5*(xmax(3)+xmin(3)) xmax(3)]));
        end
        zlabel(DVnames(3));
        axis square;
        grid on;
%         IslandsNames = cat(2,IslandsNames,num2str(r));
    end
  end
end