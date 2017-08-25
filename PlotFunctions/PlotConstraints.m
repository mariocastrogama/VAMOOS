function [hh,AX,BigAX]=PlotConstraints(islands,scales,xnames)
% plots the scatters of trade-off between constraints obtained 
% for every island
  
  global ncon
  global nislands
  global npop
  global igen
  global ngen
  [out_dot] = create_attributes(nislands);
  [out_col] = create_colors(nislands);
  [out_mrk] = create_markers(nislands);
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  IslandsNames = {};
  ifig = 2;  
  switch ncon
    case {0,1}
      disp('Too few constraints to display');
    case 2
      figure(ifig);
      clf;
      for r = 1:nislands
        XX1 = (reshape([islands(r).pop.Constraints],ncon, npop(r)));
        plot(XX1(1,:), XX1(2,:),out_mrk{r},'Marker',out_mrk{r},...
          'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'line','No'); hold on;
        set(gca,'XScale',scales{1});
        set(gca,'YScale',scales{2});
        xlabel(xnames(1));
        ylabel(xnames(2));
        IslandsNames = cat(2,IslandsNames,num2str(r));
      end
      title(['Constraint Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
      legend(IslandsNames);
    case 3
      figure(ifig);
      clf;
      for r = 1:nislands
        XX1 = (reshape([islands(r).pop.Constraints],ncon, npop(r)));
        plot3(XX1(1,:), XX1(2,:), XX1(3,:),'Marker',out_mrk{r},...
          'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'line','No'); hold on;
        set(gca,'XScale',scales{1});
        set(gca,'YScale',scales{2});
        set(gca,'ZScale',scales{3});
        xlabel(xnames(1));
        ylabel(xnames(2));
        zlabel(xnames(3));
        IslandsNames = cat(2,IslandsNames,num2str(r));
      end
      title(['Constraint Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
      legend(IslandsNames);
    case {4,5,6,7}
      figure(ifig);
      clf;
      make_it_tight = true; % make_it_tight = false;
      subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.025], [0.0825 0.075], [0.0625 0.075]);
      if ~make_it_tight; 
        clear subplot; 
      end
        
      for r = 1:nislands
        XX1 = (reshape([islands(r).pop.Constraints],ncon, npop(r)))';
%           size(XX1)
        xmin = min(XX1);
        xmax = max(XX1);
        isub = 0;
        for jcon = 1:ncon
          for icon = 1:ncon
            isub = isub + 1;
            if (icon ~= jcon)
              subplot(ncon,ncon,isub);
              hh = plot(XX1(:,icon),XX1(:,jcon),'Marker',out_mrk{r},'Markersize',5,...
              'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
              set(gca,'Fontname',sel_fontname);
              set(gca,'Fontsize',sel_fontsize);
              set(gca,'Fontweight',sel_fontweight);
              if (xmin(icon) ~= xmax(icon))
                set(gca,'xlim', 0.01*round(100*[xmin(icon) xmax(icon)]));
                set(gca,'xtick',0.01*round(100*[xmin(icon) 0.5*(xmax(icon)+xmin(icon)) xmax(icon)]));
              end
              if isub<=(ncon*(ncon-1))
                set(gca,'xticklabel',{});
              else
                xlabel(xnames(icon));
              end
              if (xmin(jcon) ~= xmax(jcon))
                set(gca,'ylim', 0.01*round(100*[xmin(jcon) xmax(jcon)]));
                set(gca,'ytick',0.01*round(100*[xmin(jcon) 0.5*(xmax(jcon)+xmin(jcon)) xmax(jcon)]));
              end
              if mod(isub,ncon)~=1;  
                set(gca,'yticklabel',{});
              else
                ylabel(xnames(jcon));
              end
              grid on;
            end % (icon ~= jcon)
          end % iobj
        end % jobj
        IslandsNames = cat(2,IslandsNames,num2str(r));
      end % r
      set(gcf,'tag',['Constraints, Gen: ',num2str(igen),' of ',num2str(ngen)]);
      legend(IslandsNames,'Position',[0.940 0.815 0.030 0.100]);
    otherwise
      disp('Too many constraints to display');
  end
  if (ncon > 1)
    drawnow;
    hold off;
    grid on;
%     title(['Constraints, Gen: ',num2str(igen),' of ',num2str(ngen)]);
  end
end