function [ I2_t_temp0 ] = getmid(I0,I2_t0)
%GETMID �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    I2_t_temp0=zeros(size(I0));
    I2_t0=trans2lefttop(I2_t0);
    
    for i=1:numel(I0(1,:))
        I2_t_temp0(:,i)=I2_t0(:,i);    
    end
end

