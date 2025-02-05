function lx = interval_x(I)
%L_X 此处显示有关此函数的摘要
% 检测图形X方向的像素数目
  Nw=sum(I(:));
  Nwall=0;
  i=1;
  while Nwall~=2*Nw
      i=i+1;
      I1=imtranslate(I,[i, 0]);
      II=I | I1;      
      Nwall=sum(II(:));      
  end
  lx=i;
end
