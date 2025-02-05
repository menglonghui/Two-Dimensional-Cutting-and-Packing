function S_min_sita = find_S_min(I)
%FIND_S_MIN �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
S=zeros(180,1);
for sita=1:180
    nl=0;%��߿ճ���
    nr=0;%�ұ߿ճ���
    nt=0;%�ϱ߿ճ���
    nb=0;%�±߿ճ���
    
    I_t=imrotate(I,sita,'crop','bicubic');
    %         shape=zeros(180,1);
    
    for i=1:numel(I_t(1,:)) %newRec�ĳ�ʼֵΪI1
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

[~, so] = sort(S, 'ascend'); %soΪ���������
% S_min=S(1);
S_min_sita=so(1);

end

