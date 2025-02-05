function [I_final] = expand_mix(I,N)
%EXPAND_MIX 此处显示有关此函数的摘要
%   此处显示详细说明
I_final = zeros(size(I));
leftpoint = 1;
rightpoint = numel(I(1,:));
toppoint = 1;
bottompoint = numel(I(:,1));

for i = toppoint : bottompoint
    for j = leftpoint : rightpoint
        
        if (i-1 >= toppoint) && (j-1 >= leftpoint) && (i-1 <= bottompoint) && (j-1 <= rightpoint)
            N1 = I(i-1,j-1);
        else
            N1 = 0;
        end
        
        if (i-1 >= toppoint) && (j >= leftpoint) && (i-1 <= bottompoint) && (j <= rightpoint)
            N2 = I(i-1,j);
        else
            N2 = 0;
        end
        
        if (i-1 >= toppoint) && (j+1 >= leftpoint) && (i-1 <= bottompoint) && (j+1 <= rightpoint)
            N3 = I(i-1,j+1);
        else
            N3 = 0;
        end
        
        if (i >= toppoint) && (j-1 >= leftpoint) && (i <= bottompoint) && (j-1 <= rightpoint)
            N4 = I(i,j-1);
        else
            N4 = 0;
        end
        
        if (i >= toppoint) && (j+1 >= leftpoint) && (i <= bottompoint) && (j+1 <= rightpoint)
            N5 = I(i,j+1);
        else
            N5 = 0;
        end
        
        if (i+1 >= toppoint) && (j-1 >= leftpoint) && (i+1 <= bottompoint) && (j-1 <= rightpoint)
            N6 = I(i+1,j-1);
        else
            N6 = 0;
        end
        
        if (i+1 >= toppoint) && (j >= leftpoint) && (i+1 <= bottompoint) && (j <= rightpoint)
            N7 = I(i+1,j);
        else
            N7 = 0;
        end
        
        if (i+1 >= toppoint) && (j+1 >= leftpoint) && (i+1 <= bottompoint) && (j+1 <= rightpoint)
            N8 = I(i+1,j+1);
        else
            N8 = 0;
        end
        
        N_around = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;
        
        top_lim = max(toppoint,i - N);
        bot_lim = min(bottompoint,i + N);
        lef_lim = max(leftpoint,j - N);
        rig_lim = min(rightpoint,j + N);
        
        if I(i,j) == 1 && (N_around < 8)
            
            for p = top_lim : bot_lim
                for q = lef_lim : rig_lim
                    if ((p - i)^2 + (q - j)^2) <= N^2 / 4
                        I_final(p,q) = 1;
                    end
                end
            end
        end
        
    end
end

I_final = I_final | I;

end