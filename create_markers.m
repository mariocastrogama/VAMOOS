function [out]=create_markers(n)
  out = {'o';... % circle
    's';... % square
    'v';... % triangle (down)
    '^';... % triangle (up)
    '<';... % triangle (left)
    '>';... % triangle (right)
    'd';... % diamond
    'p';... % pentagram
    'h'; % hexagram
    'x';... % x-mark
    '+';... % plus
    '*';... % star
    '.'};... % point
              
  out=out(1:n);
end