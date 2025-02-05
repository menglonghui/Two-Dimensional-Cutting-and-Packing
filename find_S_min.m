function S_min_sita = find_S_min(I)
%FIND_S_MIN 此处显示有关此函数的摘要
%   此处显示详细说明
S=zeros(180,1);
for sita=1:180
    nl=0;%左边空出的
    nr=0;%右边空出的
    nt=0;%上边空出的
    nb=0;%下边空出的
    
    I_t=imrotate(I,sita,'crop','bicubic');
    %         shape=zeros(180,1);
    
    for i=1:numel(I_t(1,:)) %newRec的初始值为I1
        if sum(I_t(:,i))==0
            nl=nl+1;
        else
            break
        end
    end
    
    for i=numel(I_t(1,:)):-1:1
        if sum(I_t(:,i))==0
            nr=nr+1;
        else
            break
        end
    end
    
    for i=1:numel(I_t(:,1))
        if sum(I_t(i,:))==0
            nt=nt+1;
        else
            break
        end
    end
    
    for i=numel(I_t(:,1)):-1:1
        if sum(I_t(i,:))==0
            nb=nb+1;
        else
            break
        end
    end
    S(sita)=(numel(I(1,:))-nl-nr)*(numel(I(:,1))-nt-nb);
end

[~, so] = sort(S, 'ascend'); %so为排序后的序号
% S_min=S(1);
S_min_sita=so(1);

end

