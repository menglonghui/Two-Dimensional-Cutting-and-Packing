function [dx1,dy1,dx2,dy2,dx3,dy3,sita1] = find_suit_dxy(ZP,dx1,dy1,dx2,dy2,dx3,dy3,sita1)
%FIND_SUIT_DXY 此处显示有关此函数的摘要
%   所有的操作，最终不能确定dx3的符号，但是可以确定dy3的符号为非负
if ZP == 1 %左偏
    if dx1 < 0
        dx1 = -dx1;
        dy1 = -dy1;
        sita1 = sita1 + 180;
    end
    
    if dy3 > 0%说明其下偏，但是不知道偏的是不是合适
        
        while dy3 > 0 && dy2 ~= 0
            dy3 = dy3 + dy2;
            dx3 = dx3 + dx2;
        end
        dy3 = dy3 - dy2;%保证右边是下偏最小的距离
        dx3 = dx3 - dx2;
        
    elseif dy3 < 0
        
        while dy3 < 0 && dy2 ~= 0%说明其上偏，必须让其下偏，这么做是为了在向右铺的时候，右边的最下边在分界线一下，不至于浪费空间
            dy3 = dy3 - dy2;
            dx3 = dx3 - dx2;
        end
        
    end
    
else %右偏
    
    if dx1 < 0
        dx1 = -dx1;
        dy1 = -dy1;
        sita1 = sita1 + 180;
    end
    %dx3和dy3分别为右边的块相对于左边的块的相对位移，现在是从右往左铺，需要知道的是左边的快相对于右边的块的位移
    %因此需要进行如下取相反数的操作
    %     dx3 = -dx3;
    %     dy3 = -dy3;
    
    if dy3 > 0%说明其下偏，但是不知道偏的是不是合适
        
        while dy3 > 0 && dy2 ~= 0
            dy3 = dy3 + dy2;
            dx3 = dx3 + dx2;
        end
        
        dy3 = dy3 - dy2;%保证右边是下偏最小的距离
        dx3 = dx3 - dx2;
        
    elseif dy3 < 0
        
        while dy3 < 0 && dy2 ~= 0%说明其上偏，必须使得其下偏
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

