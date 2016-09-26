function [h] = OFspace(XX,YY,OFsel,OFnames,DVnames,objtype)
% [h] = OFspace(XX,YY,OFsel,OFnames,OF,objtype)
%  This function generates the plot OF the objective space as scatter plot.
%
% Input Arguments:
%  XX         : Matrix containing the Objective functions           [nsamples, nobj]
%  YY         : Matrix containing the DV/Decision Variables [nsamples, npar]
%  OFsel      : The number OF the column OF the Objective Functions you want
%               to display
%  OFnames    : list containing the names OF the Objective Functions
%  DVnames    : list containing the names OF the DV
%  objtype    : list containing the type OF Objective Function 
%               if it is maximization ('max') then keep the axes as they are, while 
%               if it is minimization ('min') reverse the zaxis.
%
% Output Arguments
%  h          : handle OF the figure
%
% Other requirements 
%  nchoosek.m  <specfun>
% 
% Developed by
% MSc. Mario Castro Gama
% PhD Researcher UNESCO-IHE, IWSG.
% 2015-09-15
%
%
  switch nargin 
    case 2;
      [~, nobj] = size(XX);
      OFsel     = randi(nobj);
      [OFnames, ~] = create_fignames(nobj,'obj');
      [DVnames, ~] = create_fignames(nobj,'dec');
      objtype = {'min'};
      [h] = OFspace(XX,YY,OFsel,OFnames,DVnames,objtype);
    case 3;
      [~, nobj] = size(XX);
      [OFnames, ~] = create_fignames(nobj,'obj');
      [DVnames, ~] = create_fignames(nobj,'dec');
      objtype = {'min'};
      [h] = OFspace(XX,YY,OFsel,OFnames,DVnames,objtype);
    case 4
      [~, nobj] = size(XX);
      [DVnames, ~] = create_fignames(nobj,'dec');
      objtype = {'min'};
      [h] = OFspace(XX,YY,OFsel,OFnames,DVnames,objtype);
    case 5
      objtype = {'min'};
      [h] = OFspace(XX,YY,OFsel,OFnames,DVnames,objtype);
    case 6
      npar = size(YY,2);
      nobj = size(XX,2);
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
      if isempty(OFnames)
        [OFnames, ~] = create_fignames(nobj,'obj');
      end
      if isempty(DVnames)
        [DVnames, ~] = create_fignames(nobj,'dec');
      end
      
      h = figure(2);
      set(h,'Color',[1.0 1.0 1.0]);
      if strcmp(objtype(OFsel),'min');
        [BestOF,BestOFid] = min(XX(:,OFsel));
        OFlim = min(XX(:,OFsel)); %0.0;
      else
        [BestOF,BestOFid] = max(XX(:,OFsel));
        OFlim = 1.0;
      end

      for ii = 1:npar;
        for jj = ii:npar;
          if (ii ~= jj)
            isub = isub + 1;
            subplot(nsuby,nsubx,isub);
            alldata = [YY(:,ii), YY(:,jj), XX(:,OFsel)];
            x = alldata(BestOFid,1);
            y = alldata(BestOFid,2);
            z = alldata(BestOFid,3);
            x1 = [x, x, x, x, min(YY(:,ii))];
            y1 = [min(YY(:,jj)), y, y, y, y];
            z1 = [OFlim, OFlim, BestOF, OFlim, OFlim];

            plot3(alldata(:,1),alldata(:,2),alldata(:,3),'.');      hold on;
            plot3(x1,y1,z1,'r--');                                  hold on; % the line
            plot3(x,y,z,'rs','MarkerFaceColor','y','Markersize',8); hold off;
            set(gca,'Fontname','Arial');
            set(gca,'Fontsize',12);
            set(gca,'Fontweight','Bold');
            set(gca,'xlim',[])
            xlabel(DVnames(ii));
            ylabel(DVnames(jj));
            zlabel(OFnames(OFsel))
            if strcmp(objtype(OFsel),'min');
              set(gca,'zdir','reverse');
            end
            grid on;
          end
        end
      end
    otherwise
      error(' Can not display the Objectives vs Paramters');
  end % switch nargin
end % function
