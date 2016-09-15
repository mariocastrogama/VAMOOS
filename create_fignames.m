function [OUTnames,OUTscales]=create_fignames(n,nametype)
  OUTnames  = {};
  OUTscales = {};
  switch nametype
    case 'obj'
      ntext = 'OF_';
    case 'con'
      ntext = 'C_';
    case 'dec'
      ntext = 'X_';
  end
  % create names for figures
  for i = 1:n
    OUTscales = cat(2,OUTscales,'linear');
    OUTnames  = cat(2,OUTnames,[ntext,num2str(i)]);
  end
  clear i
end