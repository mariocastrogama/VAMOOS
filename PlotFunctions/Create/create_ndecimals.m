function [ndec] = create_ndecimals(XX)
  [~,nobj] = size(XX);
  logdelta = log10(abs(max(XX)-min(XX)));
  ndec = zeros(1,nobj);
  for ii = 1:nobj
    if logdelta(ii) > 2
      ndec(ii) = 0;
    else
      if logdelta(ii) > 1
        ndec(ii) = 1;
      else
        ndec(ii) = abs(floor(logdelta(ii)))+1;
      end
    end
  end
end