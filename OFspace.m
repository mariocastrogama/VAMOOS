function [h] = OFspace(OF,DV,OFsel,OFnames,DVnames,objtype)
% [h] = OFspace(OF,DV,OFsel,OFnames,OF,objtype)
%  This function generates the plot OF the objective space as scatter plot.
%
% Input Arguments:
%  OF         : Matrix containing the Objective functions           [nsamples, nobj]
%  DV         : Matrix containing the DV/Decision Variables [nsamples, npar]
%  OFsel      : The number OF the column OF the Objective Functions you want
%              to display
%  OFnames    : list containing the names OF the Objective Functions
%  DVnames    : list containing the names OF the DV
%  objtype    : list containing the type OF Objective Function 
%              if it is maximization ('max') then keep the axes as they are, while 
%              if it is minimization ('min') reverse the zaxis.
%
% Output Arguments
%  h          : handle OF the figure
%
% Other requirements 
%  nchoosek.m  <specfun>
% 
%
%
% Developed by
% MSc. Mario Castro Gama
% PhD Researcher UNESCO-IHE, IWSG.
% 2015-09-15
%
%
  nargin
  switch nargin 
    case 2;
      [~, nobj] = size(OF);
      OFsel = randi(nobj);
      [OFnames, OFscales] = create_fignames(nobj,'obj');
      [DVnames, DVscales] = create_fignames(nobj,'dec');
      objtype = {'min','min'};
    case 3;

    case 4

    case 5
      
    case 6

    otherwise
      error(' Can not display the Objectives vs Paramters');
  end
  npar = size(DV,2);
  nobj = size(OF,2);
  nsubs = nchoosek(npar,2);
  isub = 0;
  nsubx = ceil(sqrt(nsubs));
  if mod(nsubs,nsubx)==0;
    nsuby = nsubs / nsubx;
  else
    nsuby = nsubx + 1;
  end
  if (OFsel < 1)  || (OFsel > nobj)
    error(' OFsel is not a valid value');
  end
  h = figure(2);
  set(h,'Color',[1.0 1.0 1.0]);
  if strcmp(objtype(OFsel),'min');
    [BestOF,BestOFid] = min(OF(:,OFsel));
    OFlim = min(OF(:,OFsel)); %0.0;
  else
    [BestOF,BestOFid] = max(OF(:,OFsel));
    OFlim = 1.0;
  end
  
     
  for ii = 1:npar;
    for jj = ii:npar;
      if (ii ~= jj)
        isub = isub + 1;
        subplot(nsuby,nsubx,isub);
        alldata = [DV(:,ii), DV(:,jj), OF(:,OFsel)];
        x = alldata(BestOFid,1);
        y = alldata(BestOFid,2);
        z = alldata(BestOFid,3);
        x1 = [x, x, x, x, min(DV(:,ii))];
        y1 = [min(DV(:,jj)), y, y, y, y];
        z1 = [OFlim, OFlim, BestOF, OFlim, OFlim];
        
        plot3(alldata(:,1),alldata(:,2),alldata(:,3),'.');      hold on;
        plot3(x1,y1,z1,'r--');                                  hold on; % the line
        plot3(x,y,z,'rs','MarkerFaceColor','y','Markersize',8); hold off;
        set(gca,'Fontname','Arial');
        set(gca,'Fontsize',12);
        set(gca,'Fontweight','Bold');
        xlabel(DVnames(ii));
        ylabel(DVnames(jj));
        zlabel(OFnames(OFsel))
        if strcmp(objtype(OFsel),'min');
%           set(gca,'zaxis','reverse');
          set(gca,'zdir','reverse');
        end
        grid on;
%         xlabel(['par',num2str(ii)]);
%         ylabel(['par',num2str(jj)]);
      end
    end
  end
  
end