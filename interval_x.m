function lx = interval_x(I)
%L_X �˴���ʾ�йش˺�����ժҪ
% ���ͼ��X�����������Ŀ
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
