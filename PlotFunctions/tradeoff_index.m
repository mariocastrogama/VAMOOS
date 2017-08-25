function [hh, lambdas, to_names] = tradeoff_index(XX,ncolor,out_clr)
% function [ hh ] = tradeoff_index(XX,ncolor,out_clr)
% This function estimates the Beta and Lambda of the tradeoff between Objective Functions
% 
% Input Argument 
% XX       = matrix containing Objectives 
% ncolor   = number of colors to display the mosaic of tradeoff indexes
% out_clr  = selected color palette
%
% Output Arguments
% hh       = handle of the figure
% lambdas  = This matrix contains the lambdas which are estimated for the
%            comparison among objecteives 2 at a time
% to_names = total trade-off of the Objectives
%
% References
% [1] Unal M, Warn G, Simpson TW (2016) Quantifying tradeoffs to 
%     reduce the dimensionality of complex design optimization problems and 
%     expedite trade space exploration. Struct. Multidisc. Optim.  54: 233. 
%     DOI 10.1007/s00158-015-1389-7
% [2] Unal M, Warn G, Simpson TW (2015) Introduction of a tradeoff index for
%     efficient trade space exploration. ASME International Design
%     Engineering Technical Conferences, Boston, Massachusetts,
%     IDETC/CIE 2015
%
% Created by
%  MSc Mario Castro-Gama
%  PhD Researcher UNESCO-IHE / TU Delft
%  2016-09-05
%
% Still to do 
% - Input verification and include more variables such as 
% - return also the lambdas for further evaluation (if required)
%
  
  [nsols,nobj] = size(XX);
  
  % either receive names of Objectives or create new ones
  [OFnames, ~] = create_fignames(nobj,'obj');
  
  % Plotting options
  sel_fontname   = 'Arial';
  sel_fontsize   = 14;
  sel_fontweight = 'bold';
  ndecimals      = 3;
  strformat      = ['%2.',num2str(ndecimals),'f'];

  % select which monitor to use 1 or 2
  select_screen(1);
  
  % here starts the analysis of trade-offs
  XX2            = rescale(XX);
  ntradeoffs     = nchoosek(nobj,2);
  nrows          = ceil(sqrt(nsols));
  clear XX
  % create empty matrices to store results
  lambdas = zeros(nrows,nrows,ntradeoffs);
  
  to_names = {};
  ito = 0;
  for of1 = 1:nobj-1
    for of2 = (of1+1):nobj;
      ito = ito + 1;
      L = zeros(1,nrows*nrows);
      for si = 1:nsols;
        beta = zeros(1,nsols);
        for sj = 1:nsols;
          if (si ~= sj) % only important in different solutions
            xp = (XX2(si,of1) - XX2(sj,of1))*( XX2(si,of2) - XX2(sj,of2));
            if (xp < 0); % there is a crossing
              beta(1,sj) = 1.0; % values are already zero no need to put else
            end
          end % (si ~= sj)
        end % sj
        L(1,si) = sum(beta(1,:))/(nsols-1);
      end % si
      lambdas(:,:,ito) = reshape(L,nrows,nrows)';
      percentage = (sum(L)/nsols)*100;
      str = {[OFnames{of1},' - ',OFnames{of2}],['(',sprintf('%3.1f',percentage),'%)']};
      to_names = cat(2,to_names,{str});
    end % of2
  end % of1
  clear beta L of1 of2 si sj str xp percentage

  %% Plot the Mosaic of lambdas representing the tradeoff indexes
  hh = figure;
  set(hh,'Color',[1.0 1.0 1.0]);
  set(hh,'Position',[ 80 20 930 800]);
  ito = 0;
  xfalse = (0.5:nrows-0.5)/nrows;
  text_size = 18 - nobj; % this seems to solve a main configuration problem
  for of1 = 1:nobj-1;
    for of2 = (of1+1):nobj;
      ito = ito + 1;
      yp = (nobj-(of1+1)) + xfalse;
      xp = (of2-2) + xfalse;
      surf(xp,yp,lambdas(:,:,ito)); hold on;
      text(mean(xp), mean(yp),2,to_names{ito},...
        'BackgroundColor',[1 1 1],'edgecolor',[0 0 0],...
        'Fontname',sel_fontname,'FontSize',text_size,'fontweight',sel_fontweight,...
        'horizontalalignment','center','verticalalignment','middle');
    end % of2
  end % of1
  clear of1 of2 ito xp yp xfalse
  
  % Plot external square
  plot3([0 (nobj-1) (nobj-1) 0 0],[0 0 (nobj-1) (nobj-1) 0],3+zeros(1,5),'k-','linewidth',2);
  
  % set axis properties
  set(gca,'Position',[ 0.005  0.025  0.900 0.900]);
  set(gca,'xlim',[0 (nobj-1)]);
  set(gca,'xtick',[0 (nobj-1)]);
  set(gca,'xticklabel',{'',''});
  set(gca,'ylim',[0 (nobj-1)]);
  set(gca,'ytick',[0 (nobj-1)]);
  set(gca,'yticklabel',{'',''});
  shading flat;
  view(0,90); 
  axis square; 
  
  % Set title of Figure
  title('Trade-off Index \lambda', ...
    'FontSize', sel_fontsize+6,...
    'Fontweight','Bold');
  
  % Set colormap ticks equal number as divisions
  % less than 10 colors is fine, more can be an issue for visualization
  xclrtick = (0.0:(1.0/ncolor):1.0);
  xclrbar = xclrtick;

  % Create colorbar ticks labels
  xclrlabel = {};
  for ii = 1:length(xclrtick);
    xclrlabel = cat(2, xclrlabel, sprintf(strformat,xclrbar(ii)));
  end % ii
  
  % create the color palette
  [out_clr] = create_colors(ncolor,out_clr);
  colormap(out_clr); 
  caxis([0, 1]);
  
  % Set colorbar
  hc = colorbar('Ticks',xclrtick,'TickLabels',xclrlabel,'colormap',out_clr,...
    'Location','manual','Position',[0.86 0.250 0.05 0.50],...
    'Fontname','Arial','Fontsize',sel_fontsize+2,'Fontweight','Bold');

  % Set colorbar title. Thanks to the guys of http://undocumentedmatlab.com/
  % this changed recently between versions you may just need to use
  % something like...
  %
  %      get(hc,'xlabel');
  %      set(get(hc,'xlabel'),'String',['by Color : ',OFnames{bycol}]);
  %
  text('Parent', hc.DecorationContainer, ...
    'Interpreter','tex','String', '\lambda', ...
    'FontSize', sel_fontsize + 4,...
    'Fontweight','Bold',...
    'Position', [1.45, 1.08, 0], ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

end % function

