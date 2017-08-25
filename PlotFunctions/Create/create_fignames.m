function [OUTnames,OUTscales]=create_fignames(n,nametype)
  OUTnames  = {};
  OUTscales = {};
  switch nametype
    case 'obj'
      ntext = 'OF_{';
    case 'con'
      ntext = 'C_{';
    case 'dec'
      ntext = 'X_{';
    case 'type'
      ntext = 'min';
  end
  % create names for figures
  for i = 1:n
    OUTscales = cat(2,OUTscales,'linear');
    if strcmp(nametype,'type')
      OUTnames  = cat(2,OUTnames,ntext);
    else
      OUTnames  = cat(2,OUTnames,[ntext,num2str(i),'}']);
    end
  end
  clear i
end % function
