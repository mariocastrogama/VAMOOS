function [hh] = HRV_method(XX)
% function [hh] = HRV_method(XX)
% This function estiamtes the Hyper-Radial normalization of 
% a set of objectives
% 
% 0) Normalizes the matrix of objectives
% 1) Divides de objectives into 2 groups sel1 and sel2. 
%    There are ncombos = 2^(nobj-1)-1 ways to split the objectives into two
%    vectors. It creates all the possible such cases by random generation
%    and use of 'unique.m'.
% 2) estimates the norm of each subgroup of objecties with 'norm_forall.m' 
%    and normalizes the norm by the number of objectives. (1/nobj^0.5). 
%    Such norms 'euc1' and 'euc2' are always between [0-1]
% 3) It then estimates the norm between 'euc1' and 'euc2' and finds the best
%    possible solution close to the origin {0,0}
% 4) plots all the different Hyper-Radial norms
% 
% Created by: 
%   Mario Castro Gama
%   m.castrogama@unesco-ihe.org
%   PhD Researcher IWSG, UNESCO-IHE
%   Last Update: 2016.09.22
%
  [nobj] = size(XX,2);
  
  % Plotting options
%   [OFnames,OFscales]=create_fignames(nobj,'obj');
  sel_fontname = 'Arial';
  sel_fontweight = 'bold';
  ndecimals = 2;
  strformat = ['%3.',num2str(ndecimals),'f'];

  % select which monitor to use 1 or 2
  select_screen(1);

  xmin = min(XX);
  xmax = max(XX);
  XX2 = rescale(XX);
  
  norm_type = 2;
  nsub = nchoosek(nobj,2);
  
  make_it_tight = true;
  % make_it_tight = false;
  %                    subtightplot(m, n, p,          gap,         marg_h,      marg_w,varargin)
  %           subplot = @(m,n,p) subtightplot(m, n, p, [0.03 0.025], [0.0825 0.025], [0.035 0.095]);
  subplot = @(m,n,p) subtightplot(m, n, p, [0.05 0.05], [0.050 0.025], [0.025 0.025]);
  if ~make_it_tight,
    clear subplot;
  end
  
  hh = figure;
  set(hh,'Color',[1 1 1]);
  ntrials = (2^nobj)^2;
  combos_all = unique(sortrows(randi(2,ntrials,nobj)-1),'rows');
  xp = find(combos_all(:,1)==0);
  combos_all = combos_all(xp,:);
  if unique(combos_all(1,:)) == 0;
    combos_all = combos_all(2:end,:);
  end
  [~,y]=sort(sum(combos_all,2));
  combos_all = combos_all(y,:);
  nsub = size(combos_all,1);
  ncols = round(nsub^0.5)+1;
  nrows = ceil(nsub/(round(nsub^0.5)+1));
  sel_fontsize = 18-ncols;
  text_size    = 16-ncols;
  for isub = 1:nsub;
    
    subplot(nrows,ncols,isub);
    sel1 = find(combos_all(isub,:) == 0);
    sel2 = find(combos_all(isub,:) == 1);
    s1   = length(sel1);
    s2   = length(sel2);
    euc1 = norm_forall(XX2(:,sel1),norm_type);
    euc1 = euc1/(s1^.5);
    euc2 = norm_forall(XX2(:,sel2),norm_type);
    euc2 = euc2/(s2^.5);
    max_norm = max([euc1(:);euc2(:)]);
    euc3 = norm_forall([euc1, euc2],norm_type);
    [eucbest,xbest] = min(euc3);
  
    str = ['HR_{',sprintf('%d',sel2),'} vs HR_{',sprintf('%d',sel1),'} : S_{',num2str(xbest),'}'];
    pos = [-eucbest -eucbest 2*eucbest 2*eucbest];
    r1 = rectangle('Position',pos,'Curvature',[1 1],'LineWidth',1,'EdgeColor',[1 0 0]); hold on;
    plot(euc1,euc2,'.'); hold on;
    plot(euc1(xbest),euc2(xbest),'ro','markerfacecolor','y','Markersize',10,'linewidth',1);
    disp(str);
    text(0.60,0.9,str,...
      'BackgroundColor',[1 1 1],'edgecolor',[0 0 0],...
        'Fontname',sel_fontname,'FontSize',text_size,'fontweight',sel_fontweight,...
        'horizontalalignment','center','verticalalignment','middle');
    % Set current axis properties
    set(gca,'Fontname',sel_fontname);
    set(gca,'Fontsize',sel_fontsize);
    set(gca,'Fontweight',sel_fontweight);
    set(gca,'xlim',[ 0.0 1.01]);
    set(gca,'xtick',(0.0:0.50:1.0));
    set(gca,'xticklabel',{'0.0','0.5','1.0'});
    set(gca,'ylim',[ 0.0 1.01]);
    set(gca,'ytick',(0.0:0.50:1.0));
    set(gca,'yticklabel',{'0.0','0.5','1.0'});
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
    axis square;% tight;
    grid on;
  end
  
  
end