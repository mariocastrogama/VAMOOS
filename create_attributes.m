function [out]=create_attributes(n)
  icolors = {'r';... % red
             'b';... % blue                  
             'g';... % green
             'c';... % cyan         
             'm';... % magenta
             'y';... % yellow        
             'k'};... % black
          %'w'}; % white     
  imarkers = {'o';... % circle 
              'x';... % x-mark     
              '+';... % plus
              '*';... % star        
              's';... % square
              'd';... % diamond
              '.';... % point
              'v';... % triangle (down)
              '^';... % triangle (up)
              '<';... % triangle (left)
              '>';... % triangle (right)
              'p';... % pentagram
              'h'}; % hexagram
  if (n > 13*6); % maximum number of series reached
    error(['n > ',num2str(13*6), ', maximum number of series reached']);
    return;
  end
  out ={};
  kk=0;
  for jj=1:13;
    for ii=1:7;
      kk=kk+1;
      out{kk} = [icolors{ii},imarkers{jj}]; % many things can be done even rand;
    end
  end
  out=out(1:n);
end