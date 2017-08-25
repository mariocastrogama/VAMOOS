function [euc] = norm_forall(XX,norm_type)
% function [euc] = norm_forall(XX,norm_type)
% Estimate the norm of several vectors
%
% Input arguments:
% XX        = [ns,xi] matrix containing the vectors (xi) for which the norm must be estimated 
% norm_type = valid cases are '1', '2' (default), 'Inf', 'Fro'  
%
% Output arguments:
% euc       = [ns,1] matrix containing the norms of each vector
%
% Created by: 
%   Mario Castro Gama
%   m.castrogama@unesco-ihe.org
%   PhD Researcher IWSG, UNESCO-IHE
%   Last Update: 2016.09.22
%
  switch nargin
    case 0
      XX = rand(100,3);
      norm_type = 2;
      euc = norm_forall(XX,norm_type);
    case 1
      norm_type = 2;
      euc = norm_forall(XX,norm_type);
    case 2
      [ndata] = size(XX,1);
      % find the norm of the rescaled MOO solutions
      euc = zeros(ndata,1);
      for ii = 1:ndata; 
        euc(ii,1) = norm(XX(ii,:),norm_type); 
      end
  end
end