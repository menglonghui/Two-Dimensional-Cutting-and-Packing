function [nl,nr,nt,nb] = margin(I)
%MARGIN 此处显示有关此函数的摘要
%   此处显示详细说明
nl=0;
nr=0;
nt=0;
nb=0;


for i=1:numel(I(1,:)) %newRec的初始值为I1
    if sum(I(:,i))==0
        nl=nl+1;
    else
        break
    end
end  %基于分离成功后的图像，计算nl，nr，nt，nb四个值

for i=numel(I(1,:)):-1:1
    if sum(I(:,i))==0
        nr=nr+1;
    else
        break
    end
end



for i=1:numel(I(:,1)) %newRec的初始值为I1
    if sum(I(i,:))==0
        nt=nt+1;
    else
        break
    end
end  %基于分离成功后的图像，计算nl，nr，nt，nb四个值

for i=numel(I(:,1)):-1:1
    if sum(I(i,:))==0
        nb=nb+1;
    else
        break
    end
end

end

