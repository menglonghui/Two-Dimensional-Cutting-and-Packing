function [sita_best,width_min] = find_I2_t_xmin(I)
%FIND_I2_T_XMIN 此处显示有关此函数的摘要
%   此处显示详细说明
    width=zeros(1,1);%初始化
    I=trans2mid(I);
    for sita=1:2
        I_temp=imrotate(I,(sita-1)*90,'crop','bicubic');
        nl=0;
        nr=0;
        for i=1:numel(I_temp(1,:)) %newRec的初始值为I1
            if sum(I_temp(:,i))==0
                nl=nl+1;
            else
                break
            end
        end

        for i=numel(I_temp(1,:)):-1:1
            if sum(I_temp(:,i))==0
                nr=nr+1;
            else
                break
            end
        end

        width(sita,1)=sita;
        width(sita,2)=numel(I_temp(1,:))-nl-nr;

    end
    
    [~, so] = sort(width(:,2), 'ascend'); %so为降序后的序号排序
    sita_best=(width(so(1),1)-1)*90;
    width_min=width(so(1),2);

end

