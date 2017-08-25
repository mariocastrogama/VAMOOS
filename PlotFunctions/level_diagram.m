function [hh] = level_diagram(XX)
% function [ hh ] = level_diagram(XX)
% This function estimates the norm of the Objective Functions w.r.t. the
% best normalized solution
% 
% Input Argument 
% XX = matrix containing Objectives 
%
% Output Arguments
% hh = handle of the figure
%
% Referecences
% [1] Blasco X, Herrero JM, Sanchis J, Martínez M (2008) A new graphical
%     visualization of n-dimensional Pareto front for decision-making in
%     multiobjective optimization. Inf Sci 178(20):3908–3924
%
% Created by
% MSc Mario Castro-Gama
% PhD Researcher UNESCO-IHE / TU Delft
% 2016-09-05
%
% Still to do 
% - Input verification and include more variables such as 
%  + return also the levels for further evaluation (if required)
%
  [ndata, nobj] = size(XX);
  
  % Plotting options
  [OFnames,OFscales] = create_fignames(nobj,'obj');
  sel_fontname       = 'Arial';
  sel_fontsize       = 12;
  sel_fontweight     = 'bold';
  ndecimals          = 3;
  strformat          = ['%3.',num2str(ndecimals),'f'];

  % select which monitor to use 1 or 2
  select_screen(1);

  xmin = min(XX);
  xmax = max(XX);
  XX2  = rescale(XX);

  % find the norm of the rescaled MOO solutions
  norm_type        = 2;
  [euc]            = norm_forall(XX2,norm_type);
  BTOmax           = max(euc);
  [BTOmin, BTOpos] = min(euc);

 %% Plot the level diagram of each objective w.r.t. the norm
  make_it_tight = true;
  % make_it_tight = false;
  %                    subtightplot(m, n, p,          gap,         marg_h,      marg_w,varargin)
  %           subplot = @(m,n,p) subtightplot(m, n, p, [0.03 0.025], [0.0825 0.025], [0.035 0.095]);
  subplot = @(m,n,p) subtightplot(m, n, p, [0.060 0.030], [0.0725 0.05], [0.02 0.01]);
  if ~make_it_tight,
    clear subplot;
  end

  hh = figure;
  set(hh,'Color',[1 1 1]);
  ncols = ceil(nobj/2);
  for isub = 1:nobj
    subplot(2,ncols,isub);
    plot(XX(:,isub),euc,'b.'); hold on;
    plot(round(XX(BTOpos,isub),ndecimals),round(BTOmin,ndecimals),'rs','Markerfacecolor','y','Markersize',10);
    % Set current axis properties
    set(gca,'Fontname',sel_fontname);
    set(gca,'Fontsize',sel_fontsize);
    set(gca,'Fontweight',sel_fontweight);

    % Set x limits, tick, labelrotation, scale
    set(gca,'xlim',  round([xmin(isub) xmax(isub)],ndecimals));
    set(gca,'xtick', round(xmin(isub) + (0.00:0.50:1.00).*(xmax(isub)-xmin(isub)),ndecimals));

    set(gca,'XScale',OFscales{isub});
    xlab = xlabel(OFnames{isub});
    set(xlab,'HorizontalAlignment','left');
    set(xlab,'VerticalAlignment','cap'); % 'baseline' 'top' 'cap' 'middle' 'bottom'

    % Set y limits, tick, scale

    if mod(isub,ncols) == 1
      set(gca,'ylim',  round([BTOmin, BTOmax],ndecimals));
      set(gca,'ytick', round(BTOmin + (0.00:0.25:1.00).*(BTOmax - BTOmin),ndecimals));
      set(gca,'YScale','log');
      ylab = ylabel('2-norm');
      set(ylab,'HorizontalAlignment','left');
    else
      set(gca,'ylim',  round([BTOmin, BTOmax],ndecimals));
      set(gca,'ytick', round(BTOmin + (0.00:0.25:1.00).*(BTOmax - BTOmin),ndecimals));
      set(gca,'yticklabel',{'','','','',''});
    end
    axis square;
    grid on;
  end % isub
end % function

