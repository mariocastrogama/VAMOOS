function [ hh ] = tradeoff_index(XX)
%TRADEOFF_INDEX Summary of this function goes here
%   Detailed explanation goes here
  
  
  [nsol,nobj] = size(XX);
  
  % either receive names of Objectives or create new ones
  [OFnames, ~]=create_fignames(nobj,'obj');
  
  % Receive the number of colors required
  ncolor = 10;
  
  % Receive the color palette to display the values
  out_clr = [1.0 0.2 0; 1 1 1;0 0.2 1];
  
  % Plotting options
  sel_fontname = 'Arial';
  sel_fontsize = 10;
  sel_fontweight = 'bold';
  ndecimals = 3;
  strformat = ['%4.',num2str(ndecimals),'f'];

  % select which monitor to use 1 or 2
  select_screen(1);
  
  % here starts the analysis of trade-offs
  xmin = min(XX);
  xmax = max(XX);
  XX2  = rescale(XX);
  ntradeoffs = nchoosek(nobj,2);
  nrows      = ceil(sqrt(nsol));
  % Main cycle of calculation of tradeoff index
  nlambdas = nsol*(nsol-1)/2;
  
  % create empty matrices to store results
  beta   = zeros(nsol-1,nsol-1,ntradeoffs);
  lambda = zeros(nrows,nrows,ntradeoffs);
  
  to_names = {};
  ito = 0;
  for of1 = 1:nobj-1
    for of2 = (of1+1):nobj;
      ito = ito + 1;
      L = zeros(1,nrows*nrows);
      for si = 1:nsol;
        for sj = 1:nsol;
          if (si ~= sj)
            xp = (XX2(si,of1) - XX2(sj,of1))*( XX2(si,of2) - XX2(sj,of2));
            if xp < 0; % there is a crossing
              beta(si,sj,ito) = 1.0;
            end
          end % (ii ~= jj)
        end % jj
        L(1,si) = sum(beta(si,:,ito))/(nsol-1);
      end % ii
      lambda(:,:,ito) = reshape(L,nrows,nrows)';
      ncrossings = length(find(beta(:,:,ito)==1))/2;
      percentage = (ncrossings/nlambdas)*100;
      str = [OFnames{of1},' - ',OFnames{of2},' (',sprintf('%3.1f',percentage),'%)'];
      to_names = cat(2,to_names,str);
    end % of2
  end % of1
  clear beta
  
  %% Plot the Mosaic of lambdas representing the tradeoff indexes
  hh = figure;
  set(hh,'Color',[1.0 1.0 1.0]);
  ito = 0;
  xfalse = (0.5:nrows-0.5)/nrows;
  for of1 = 1:nobj-1;
    for of2 = (of1+1):nobj;
      ito = ito + 1;
      yp = (nobj-(of1+1)) + xfalse;
      xp = (of2-2) + xfalse;
      surf(xp,yp,lambda(:,:,ito)); hold on;
      t1 = text(mean(xp), mean(yp),2,to_names{ito},...
        'BackgroundColor',[1 1 1],'edgecolor',[0 0 0],...
        'FontSize',sel_fontsize,'fontweight',sel_fontweight,...
        'horizontalalignment','center','verticalalignment','middle');
    end % of2
  end % of1
  plot([0 (nobj-1) (nobj-1) 0 0],[0 0 (nobj-1) (nobj-1) 0],'k-','linewidth',2);
  set(gca,'xlim',[0 (nobj-1)]);
  set(gca,'xtick',[0 (nobj-1)]);
  set(gca,'xticklabel',{'',''});
  set(gca,'ylim',[0 (nobj-1)]);
  set(gca,'ytick',[0 (nobj-1)]);
  set(gca,'yticklabel',{'',''});
  shading flat;
  [out_clr]=create_colors(ncolor,out_clr);
  colormap(out_clr); 
  view(0,90); 
  axis square; 

  title('Trade-off Index \lambda', ...
    'FontSize', sel_fontsize+6,...
    'Fontweight','Bold');
  % Set colormap ticks equal number as divisions
  % less than 10 colors is fine
  xclrtick = (0:(1.0/ncolor):1.0);
  xclrbar = xclrtick;

  % Create colorbar ticks labels
  xclrlabel = {};
  for ii = 1:length(xclrtick);
    xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
  end % ii

  % Set colorbar
  hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,...
    'Location','manual','Position',[0.750 0.25 0.05 0.50],'colormap',out_clr,...
    'Fontsize',sel_fontsize+2,'Fontweight','Bold');

  % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
  % this changed recently between versions you may just need to use
  % something like...
  %
  %      get(hc,'xlabel');
  %      set(get(hc,'xlabel'),'String',['by Color : ',OFnames{bycol}]);
  %
  text('Parent', hc.DecorationContainer, ...
    'Interpreter','tex','String', '\lambda', ...
    'FontSize', sel_fontsize+4,...
    'Fontweight','Bold',...
    'Position', [1.35, 1.08, 0], ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

end % function

