function [nl,nr,nt,nb] = margin(I)
%MARGIN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
nl=0;
nr=0;
nt=0;
nb=0;


for i=1:numel(I(1,:)) %newRec�ĳ�ʼֵΪI1
    if sum(I(:,i))==0
        nl=nl+1;
    else
        break
    end
end  %���ڷ���ɹ����ͼ�񣬼���nl��nr��nt��nb�ĸ�ֵ

for i=numel(I(1,:)):-1:1
    if sum(I(:,i))==0
        nr=nr+1;
    else
        break
    end
end



for i=1:numel(I(:,1)) %newRec�ĳ�ʼֵΪI1
    if sum(I(i,:))==0
        nt=nt+1;
    else
        break
    end
end  %���ڷ���ɹ����ͼ�񣬼���nl��nr��nt��nb�ĸ�ֵ

for i=numel(I(:,1)):-1:1
    if sum(I(i,:))==0
        nb=nb+1;
    else
        break
    end
end

end

