function [DVscales,COscales,OFscales,DVnames,OFnames,COnames]=create_namesall
  global nvar
  global ncon
  global nobj
  DVscales = {};
  DVnames =  {};
  for ivar = 1:nvar
    DVscales = cat(2,DVscales,'linear');
    DVnames  = cat(2,DVnames,['X_',num2str(ivar)]);
  end
  clear ivar
  COscales = {};
  COnames =  {};
  for icon =1:ncon
    COscales = cat(2,COscales,'linear');
    COnames  = cat(2,COnames,['C_',num2str(icon)]);
  end
  clear icon
  OFscales = {};
  OFnames  = {};
  for iobj=1:nobj
    OFscales = cat(2,OFscales,'linear');
    OFnames  = cat(2,OFnames,['OF_',num2str(iobj)]);
  end
  clear iobj
end