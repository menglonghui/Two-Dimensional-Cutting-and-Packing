function Imid = trans2mid(I,jd)
%��ͼ��������м�λ��
I = trans2lefttop(I);
nr = numel(I(1,:)) - l_x(I);
nb = numel(I(:,1)) - l_y(I);
% nl=0;%��߿ճ���
% nr=0;%�ұ߿ճ���
% nt=0;%�ϱ߿ճ���
% nb=0;%�±߿ճ���

% for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
%     if sum(I(:,i))==0
%         nl=nl+1;
%     else
%         break
%     end
% end

% for i=numel(I(1,:)):-1:1
%     if sum(I(:,i))==0
%         nr=nr+1;
%     else
%         break
%     end
% end

% for i=1:numel(I(:,1))
%     if sum(I(i,:))==0
%         nt=nt+1;
%     else
%         break
%     end
% end

% for i=numel(I(:,1)):-1:1
%     if sum(I(i,:))==0
%         nb=nb+1;
%     else
%         break
%     end
%
if jd == 0
    
    Imid=imtranslate(I,[round(nr/2), round(nb/2)]);
    
elseif jd == 1
    
    if mod(nr,2) == 0
        x_translate = nr / 2;
        x_add = 0;
    else
        x_translate = (nr+1) / 2;
        x_add = 1;
    end
    
    if mod(nb,2) == 0
        y_translate = nb / 2;
        y_add = 0;
    else
        y_translate = (nb+1) / 2;
        y_add = 1;
    end
    
    I_base = zeros(numel(I(:,1)) + y_add, numel(I(1,:)) + x_add);
    Imid = paste_img(I_base,I,1,1);
    Imid = trans2lefttop(Imid);
    Imid = imtranslate(Imid,[x_translate, y_translate]);
    
elseif jd == 2
    
    Imid=imtranslate(I,[nr/2, nb/2]);
    
end

end

