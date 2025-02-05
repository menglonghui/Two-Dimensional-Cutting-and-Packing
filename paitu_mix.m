function [coordinate,I_mix,N] = paitu_mix(Nmax,key,sita1,I2_t,I_mix,coordinate_single_pre,tx,find_best_dxy5_used,sita_xz) %#ok<INUSL>
%PAITU 此处显示有关此函数的摘要
%   此处显示详细说明
%I_mix用来铺零件，该零件不与自己做比较判断，只与之前排好的零件做比较,之前排好零件的板子就是I_ref
coordinate_single_pre(:,3) = coordinate_single_pre(:,3) - sita_xz;
N = 0;
changdu = numel(I_mix(:,1));
kuandu = numel(I_mix(1,:));

if mod(changdu,2) == 0 && mod(kuandu,2) == 0
    BS = 2;
else
    BS = 1;
end

x_m = numel(I2_t(1,:)) / 2;
y_m = numel(I2_t(:,1)) / 2;
I2_t = trans2lefttop(I2_t);
I2_t = imcrop(I2_t,[0, 0, l_x(I2_t), l_y(I2_t)]);
I_base = zeros(size(I_mix));
I2_t = paste_img(I_base,I2_t,1,1);

% I_mix = paste_img(I_base,I_mix,4,4);
I2_t = double(logical(I2_t));
I_mix = double(logical(I_mix));

I2_t = trans2mid(I2_t,2);
I10 = imrotate(I2_t,sita1,'crop','bicubic');
if find_best_dxy5_used == 0
    
    if tx == 0
        I20 = imrotate(I2_t,sita1 + 180,'crop','bicubic');
    else
        I20 = I10;
    end
    
else
    
    I20 = imrotate(I2_t,90,'crop','bicubic');
    
end

I10 = imresize(I10,1/BS);
I20 = imresize(I20,1/BS);
I_mix = imresize(I_mix,1/BS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%先进行缩小
I10(I10 < 0.5) = 0;
I20(I20 < 0.5) = 0;
I_mix(I_mix < 0.5) = 0;

I10 = logical(I10);
I20 = logical(I20);
I_mix = logical(I_mix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_ref = I_mix;
I_ref = logical(I_ref);

N_I_ref = sum(I_ref(:));
% N1 = sum(I10(:));
% N2 = sum(I20(:));
mm = 0;
% coordinate1 = [];
% coordinate2 = [];

for k = 1 : length(coordinate_single_pre)
    
    if coordinate_single_pre(k,3) == sita1
        dx=(coordinate_single_pre(k,1) - x_m) / BS;
        dy=-(coordinate_single_pre(k,2) - y_m) / BS;
        I1 = imtranslate(I10,[round(dx),round(dy)]);
        N1 = sum(I1(:));
        I_combine = I1 | I_ref;
        
        if abs(sum(I_combine(:)) - N1 - N_I_ref) <= 2 %说明并未重叠
            I_mix = I_mix | I1;
            N = N + 1;
            mm=mm+1;
            coordinate(mm,:) = coordinate_single_pre(k,:);  %#ok<AGROW>
            if eval(['N == Nmax.' key])
                break
            end
        end
        
    elseif (coordinate_single_pre(k,3) == sita1 + 180 && find_best_dxy5_used == 0) || (coordinate_single_pre(k,3) == 90 && find_best_dxy5_used == 1)
        dx=(coordinate_single_pre(k,1) - x_m) / BS;
        dy=-(coordinate_single_pre(k,2) - y_m) / BS;
        I2 = imtranslate(I20,[round(dx),round(dy)]);
        N2 = sum(I2(:));
        I_combine = I2 | I_ref;
        
        if abs(sum(I_combine(:)) - N2 - N_I_ref) <= 2%说明并未重叠
            I_mix = I_mix | I2;
            N = N + 1;
            mm=mm+1;
            coordinate(mm,:) = coordinate_single_pre(k,:); %#ok<AGROW>
            if eval(['N == Nmax.' key])
                break
            end
        end
        
    else %上面追加的零件
        
        I30 = imrotate(I2_t,coordinate_single_pre(k,3),'crop','bicubic');
        I30(I30 < 0.5) = 0;
        I30 = logical(I30);
        I30 = imresize(I30, 1 / BS);
        
        dx=(coordinate_single_pre(k,1) - x_m) / BS;
        dy=-(coordinate_single_pre(k,2) - y_m) / BS;
        I3 = imtranslate(I30,[round(dx),round(dy)]);
        N3 = sum(I3(:));
        I_combine = I3 | I_ref;
        
        if abs(sum(I_combine(:)) - N3 - N_I_ref) <= 2 %说明并未重叠
            I_mix = I_mix | I3;
            N = N + 1;
            mm=mm+1;
            coordinate(mm,:)=coordinate_single_pre(k,:);  %#ok<AGROW>
            if eval(['N == Nmax.' key])
                break
            end
        end
    end
end

if exist('coordinate','var')==0%如果该值不存在
    coordinate = [];
else
    coordinate(:,3) = coordinate(:,3) + sita_xz;
end
I_mix =double(I_mix);
I_mix = imresize(I_mix,BS);
I_mix(I_mix < 0.5) = 0;
I_mix = logical(I_mix);
end