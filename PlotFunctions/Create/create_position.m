function [Pix_SS] = create_position()
  set(0,'units','pixels') 
  Pix_SS = get(0,'screensize');
  Pix_SS = Pix_SS + [0 29 0 -90];
end