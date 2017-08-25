function [hh]=PlotObjectivesTrend(islands,xnames,bycol)
% plots the scatters of trade-off between objective functions obtained 
% for every island
  global nobj;
  global nislands;
  global npop;
  global ngen;
  global igen;
  
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'Bold';
%   ifig = 3;
%   IslandsNames = {};
  
  if (isstruct(islands))
    for r = 1:nislands;
      [out_col] = create_colors(npop(r));
      XX = (reshape([islands(r).pop.ObjectiveFunctions],nobj, npop(r)))';
      XX = sortrows(XX,bycol);
      XX = rescale(XX);
      for ipop = 1:npop(r);
        plot((1:nobj)',XX(ipop,:)','color',out_col(ipop,:)); hold on;
      end
    end
  else
    XX = islands;
    [nsamples, ncols] = size(XX);
    [out_col] = create_colors(nsamples);
    XX = sortrows(XX,bycol);
    XX = rescale(XX);
    for ipop = 1:nsamples;
      plot((1:ncols)',XX(ipop,:)','color',out_col(ipop,:)); hold on;
    end
  end
  options.Fontname = sel_fontname;
  set(gca,options); 
%   set(gca,'Fontname',);
  set(gca,'Fontsize',sel_fontsize);
  set(gca,'Fontweight',sel_fontweight);
  set(gca,'xtick',(1:nobj));
  set(gca,'xticklabel',xnames);
  set(gca,'ytick',[0.0, 1.0]);
  set(gca,'yticklabel',{'min','max'});
  set(gca,'Position',[0.065 0.085 0.90 0.85]);
  xlabel('Objective Function');
  title('Parallel Objectives');
  hold off;
end