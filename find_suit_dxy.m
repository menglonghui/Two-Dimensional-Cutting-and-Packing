function [dx1,dy1,dx2,dy2,dx3,dy3,sita1] = find_suit_dxy(ZP,dx1,dy1,dx2,dy2,dx3,dy3,sita1)
%FIND_SUIT_DXY �˴���ʾ�йش˺�����ժҪ
%   ���еĲ��������ղ���ȷ��dx3�ķ��ţ����ǿ���ȷ��dy3�ķ���Ϊ�Ǹ�
if ZP == 1 %��ƫ
    if dx1 < 0
        dx1 = -dx1;
        dy1 = -dy1;
        sita1 = sita1 + 180;
    end
    
    if dy3 > 0%˵������ƫ�����ǲ�֪��ƫ���ǲ��Ǻ���
        
        while dy3 > 0 && dy2 ~= 0
            dy3 = dy3 + dy2;
            dx3 = dx3 + dx2;
        end
        dy3 = dy3 - dy2;%��֤�ұ�����ƫ��С�ľ���
        dx3 = dx3 - dx2;
        
    elseif dy3 < 0
        
        while dy3 < 0 && dy2 ~= 0%˵������ƫ������������ƫ����ô����Ϊ���������̵�ʱ���ұߵ����±��ڷֽ���һ�£��������˷ѿռ�
            dy3 = dy3 - dy2;
            dx3 = dx3 - dx2;
        end
        
    end
    
else %��ƫ
    
    if dx1 < 0
        dx1 = -dx1;
        dy1 = -dy1;
        sita1 = sita1 + 180;
    end
    %dx3��dy3�ֱ�Ϊ�ұߵĿ��������ߵĿ�����λ�ƣ������Ǵ��������̣���Ҫ֪��������ߵĿ�������ұߵĿ��λ��
    %�����Ҫ��������ȡ�෴���Ĳ���
    %     dx3 = -dx3;
    %     dy3 = -dy3;
    
    if dy3 > 0%˵������ƫ�����ǲ�֪��ƫ���ǲ��Ǻ���
        
        while dy3 > 0 && dy2 ~= 0
            dy3 = dy3 + dy2;
            dx3 = dx3 + dx2;
        end
        
        dy3 = dy3 - dy2;%��֤�ұ�����ƫ��С�ľ���
        dx3 = dx3 - dx2;
        
    elseif dy3 < 0
        
        while dy3 < 0 && dy2 ~= 0%˵������ƫ������ʹ������ƫ
            dy3 = dy3 - dy2;
            dx3 = dx3 - dx2;
            
        end        
    end    
end

if sita1 > 180
    sita1 = sita1 - 360;
elseif sita1 < -180
    sita1 = sita1 + 360;
end

end

