function [hh] = PlotObjectives(islands,xnames,scales)
% plots the scatters of trade-off between objective functions obtained 
% for every island
  global nobj;
  global nislands;
  global npop;
  global ngen;
  global igen;
  
%   [out_dot] = create_attributes(nislands);
  if ~exist('nislands','var');
    nislands = 1;
  end
  disp(nislands);
  
  sel_fontname = 'Arial';
  sel_fontsize = 14;
  sel_fontweight = 'Bold';
  ifig = 3;
  IslandsNames = {};
  if isstruct(islands)
    [out_col] = create_colors(nislands);
    [out_mrk] = create_markers(nislands);
    switch nobj
      case {0,1}
        disp('Too few decision variables to display');
      case 2
        figure(ifig);
        clf;
        for r = 1:nislands;
          XX1 = reshape([islands(r).pop.ObjectiveFunctions],nobj, numel(islands(r).pop));%npop(r)));
          plot(XX1(1,:), XX1(2,:),'Marker',out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;%out_dot{r}); hold on;
          set(gca,'XScale',scales{1});
          set(gca,'YScale',scales{2});
          xlabel(xnames(1));
          ylabel(xnames(2));
          IslandsNames = cat(2,IslandsNames,num2str(r));
          title(['Objective Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
          legend(IslandsNames);
        end
      case 3
        figure(ifig);
        clf;
        for r = 1:nislands
          XX1 = (reshape([islands(r).pop.ObjectiveFunctions],nobj, npop(r)));
          plot3(XX1(1,:), XX1(2,:), XX1(3,:),'Marker',out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
          set(gca,'XScale',scales{1});
          set(gca,'YScale',scales{2});
          set(gca,'ZScale',scales{3});
          xlabel(xnames(1));
          ylabel(xnames(2));
          zlabel(xnames(3));
          IslandsNames = cat(2,IslandsNames,num2str(r));
        end
        title(['Objective Functions, Gen: ',num2str(igen),' of ',num2str(ngen)]);
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
          XX1 = (reshape([islands(r).pop.ObjectiveFunctions],nobj, npop(r)))';
          xmin = min(XX1);
          xmax = max(XX1);
          isub = 0;
          for jobj = 1:nobj
            for iobj = 1:nobj
              isub = isub + 1;
              if (iobj ~= jobj)
                subplot(nobj,nobj,isub);
                hh = plot(XX1(:,iobj),XX1(:,jobj),'Marker',out_mrk{r},'Markersize',4,...
                'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;
                set(gca,'Fontname',sel_fontname);
                set(gca,'Fontsize',sel_fontsize);
                set(gca,'Fontweight',sel_fontweight);
                set(gca,'xlim', 0.01*round(100*[xmin(iobj) xmax(iobj)]));
                set(gca,'xtick',0.01*round(100*[xmin(iobj) 0.5*(xmax(iobj)+xmin(iobj)) xmax(iobj)]));
                if isub<=(nobj*(nobj-1))
                  set(gca,'xticklabel',{});
                else
                  xlabel(xnames(iobj));
                end
                set(gca,'ylim', 0.01*round(100*[xmin(jobj) xmax(jobj)]));
                set(gca,'ytick',0.01*round(100*[xmin(jobj) 0.5*(xmax(jobj)+xmin(jobj)) xmax(jobj)]));
                if mod(isub,nobj)~=1;  
                  set(gca,'yticklabel',{});
                else
                  ylabel(xnames(jobj));
                end
                grid on;
              end
            end % iobj
          end % jobj
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
    [out_col] = create_colors(1);
    [out_mrk] = create_markers(1);
    switch nobj
      case {0,1}
        disp('Too few decision variables to display');
      case 2
        figure(ifig);
        r = 1;
        XX1 = islands;
        hh = plot(XX1(:,1), XX1(:,2),'Marker',out_mrk{r},'MarkerEdgecolor','k','Markerfacecolor',out_col(r,:),'lineStyle','No'); hold on;%out_dot{r}); hold on;
        set(gca,'Fontname','Arial');
        set(gca,'Fontweight','Bold');
        set(gca,'Fontsize',14);
        set(gca,'XScale',scales{1});
        set(gca,'YScale',scales{2});
        xlabel(xnames(1));
        ylabel(xnames(2));
        IslandsNames = cat(2,IslandsNames,num2str(r));
        title('Objective Functions');%, Gen: ',num2str(igen),' of ',num2str(ngen)]);
        legend(IslandsNames);
      case {3, 4, 5, 6, 7}
        make_it_tight = true;         
        % make_it_tight = false;
        %                  subtightplot(m, n, p,          gap,         marg_h,      marg_w,varargin)
        subplot = @(m,n,p) subtightplot(m, n, p, [0.03 0.025], [0.0825 0.025], [0.035 0.0725]);
        if ~make_it_tight,  
          clear subplot;  
        end
        hh = figure(ifig);
        bycol  = 3;
        bysize = 2;
        clf;
        XX1 = islands;
        ndata = size(XX1,1);
        [out_clr] = create_colors(6);
        [out_clr] = create_colors(ndata,out_clr);
        XX1 = sortrows(XX1,bycol);
        xmin = min(XX1);
        xmax = max(XX1);
        out_size = create_sizes(XX1,bysize,6);
        size(XX1);
        isub = 0;
        for jobj = 1:nobj
          for iobj = 1:nobj
            isub = isub + 1;
            subplot(nobj,nobj,isub);
            if iobj ~= jobj
              for idata = 1:ndata 
                plot(XX1(idata,iobj),XX1(idata,jobj),'Marker',out_mrk{1},'Markersize',out_size(idata),'MarkerEdgecolor','k','Markerfacecolor',out_clr(idata,:),'lineStyle','No'); hold on;
              end
            else
              set(gca,'xlim', 0.01*round(100*[xmin(iobj) xmax(iobj)]));
              set(gca,'xtick',0.01*round(100*[xmin(iobj) xmax(iobj)]));
              set(gca,'ylim',0.01*round(100*[xmin(jobj) xmax(jobj)]));
              set(gca,'ytick',0.01*round(100*[xmin(jobj) xmax(jobj)]));
              text(0.5*(xmax(iobj)+xmin(iobj)), 0.5*(xmax(jobj)+xmin(jobj)), xnames(jobj),...
                'HorizontalAlignment','center',...
                'VerticalAlignment','middle');
            end % iobj ~= jobj
            set(gca,'Fontname',sel_fontname);
            set(gca,'Fontsize',sel_fontsize);
            set(gca,'Fontweight',sel_fontweight);
            set(gca,'xlim', 0.01*round(100*[xmin(iobj) xmax(iobj)]));
            set(gca,'xtick',0.01*round(100*[xmin(iobj) 0.5*(xmax(iobj)+xmin(iobj)) xmax(iobj)]));
            set(gca,'ylim', 0.01*round(100*[xmin(jobj) xmax(jobj)]));
            set(gca,'ytick',0.01*round(100*[xmin(jobj) 0.5*(xmax(jobj)+xmin(jobj)) xmax(jobj)]));
            set(gca,'XScale',scales{iobj});
            set(gca,'YScale',scales{jobj});
            if isub <=(nobj*(nobj-1))
              set(gca,'xticklabel',{});
            end
            if mod(isub,nobj)~=1;
              set(gca,'yticklabel',{});
            end
            grid on;
          end % iobj
        end % jobj
      otherwise
        disp('Too many objectives to display');
    end
    hold off;
    grid on;
    drawnow;
  end
end