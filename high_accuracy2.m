function [dx222,dy222,I2_hr] = high_accuracy2(I1_hr,alg2,key,dx2,dy2,imrf,ijmax,RF) %#ok<*INUSL>
%HIGH_ACCURACY2 此处显示有关此函数的摘要
%   此处显示详细说明
dx2_1_last = 0;
dy2_1_last = 0;

if ijmax < 1/RF
    ijmax = 1/RF;
end
% found1 = 1;
% found3 = 1;
I1_hr = logical(I1_hr);
% N_s_c = sum(I1_hr(:));
% Alg2 = eval(['alg2.' key ';']);
if isstruct(alg2)
    Alg2 = eval(['alg2.' key]);
else
    Alg2 = alg2;
end
%第1次寻找
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dx2 >= 0
    %     IC1 = I1_hr;
    %     IC1_r = imrotate(IC1,180);
    
    IC1 = trans2leftbottom(I1_hr);
    IC1_r = IC1;
    %     IC1_r = trans2leftbottom(I1_hr);
    
    IC1 = imtranslate(IC1,[10*2*imrf+1,-10*2*imrf-1]);
    IC1_r = imtranslate(IC1_r,[10*2*imrf+1,-10*2*imrf-1]);
    
    N1 = sum(IC1(:));
    N2 = sum(IC1_r(:));
else
    %     IC1 = I1_hr;
    %     IC1_r = imrotate(IC1,180);
    
    IC1 = trans2rightbottom(I1_hr);
    IC1_r = IC1;
    %     IC1_r = trans2rightbottom(IC1_r);
    
    IC1=imtranslate(IC1,[-10*2*imrf-1,-10*2*imrf-1]);
    IC1_r=imtranslate(IC1_r,[-10*2*imrf-1,-10*2*imrf-1]);
    
    N1=sum(IC1(:));
    N2=sum(IC1_r(:));
end

if dy2 > 0
    IC1 = imtranslate(IC1,[0,-dy2*imrf]);
    IC1_r = imtranslate(IC1_r,[0,-dy2*imrf]);
end

IC1_r=imtranslate(IC1_r,[round(dx2*imrf),round(dy2*imrf)]);

while sum(IC1_r(:)) < sum(IC1(:))
    
    I_base = zeros(2*size(IC1));
    IC1 = paste_img(I_base,IC1,1,1);
    
    if sum(IC1(:)) == sum(I1_hr(:))
        IC1 = paste_img(I_base,IC1,1,1);
    else
        IC1 = paste_img(I_base,I1_hr,1,1);
    end    
    
    if dx2 >= 0     
        
        IC1 = trans2leftbottom(IC1);
        IC1_r = IC1;
        
        IC1 = imtranslate(IC1,[10*2*imrf+1,-10*2*imrf-1]);
        IC1_r = imtranslate(IC1_r,[10*2*imrf+1,-10*2*imrf-1]);
        
        N1 = sum(IC1(:));
        N2 = sum(IC1_r(:));
        
    else        
        
        IC1 = trans2rightbottom(IC1);
        IC1_r = IC1;
        
        IC1=imtranslate(IC1,[-10*2*imrf-1,-10*2*imrf-1]);
        IC1_r=imtranslate(IC1_r,[-10*2*imrf-1,-10*2*imrf-1]);
        
        N1=sum(IC1(:));
        N2=sum(IC1_r(:));
    end
    
    if dy2 > 0
        IC1 = imtranslate(IC1,[0,-dy2*imrf]);
        IC1_r = imtranslate(IC1_r,[0,-dy2*imrf]);
    end
    
    IC1_r=imtranslate(IC1_r,[round(dx2*imrf),round(dy2*imrf)]);
end

S_c = zeros(1,3);
k=0;

% IC1 = logical(IC1);
% IC1_r = logical(IC1_r);

I=IC1 | IC1_r;
% imshow(I)
i_max = ceil(min(round(l_x(I)/(3*imrf)),ijmax)/2);
j_max = ceil(min(round(l_y(I)/(3*imrf)),ijmax)/2);

%分段第一次
for i=-i_max:i_max
    for j=-j_max:j_max
        k=k+1;
        IC1_r_t= imtranslate(IC1_r,[i*imrf*2,j*imrf*2]);
        IC_t=IC1 | IC1_r_t;
        %imshow(IC_t);
        S_c(k,1)=i;
        S_c(k,2)=j;
        
        if sum(IC_t(:)) == N1 + N2 %说明图像完全分离
            if Alg2==1
                S_c(k,3)=(-1)*fun_XY(IC_t);
            elseif Alg2==2
                S_c(k,3)=(-1)*fun_YX(IC_t);
            elseif Alg2==3
                S_c(k,3)=(-1)*fun_bs2(IC_t,IC1)+(-1)*fun_Variance(IC_t)/100000000;
            elseif Alg2==4
                S_c(k,3)=(-1)*fun_Variance(IC_t)+(-1)*fun_bs2(IC_t,IC1)/100000000;
            end
        else
            S_c(k,3)=NaN;
        end
    end
end

m = 1;
while all(isnan(S_c(:,3)))%判断第三列是不是全是NaN
    %首先将画布放大
    I_base = zeros(numel(IC1(:,1)) + 8 * imrf,numel(IC1(1,:)) + 8 * imrf);
    IC1 = paste_img(I_base,IC1,4 * imrf + 1,4 * imrf + 1); IC1 = logical(IC1);
    IC1_r = paste_img(I_base,IC1_r,4 * imrf + 1,4 * imrf + 1); IC1_r = logical(IC1_r);
    
    for i = (-i_max - 4 * m) : (i_max + 4 * m)
        for j = (-j_max - 4 * m) : (j_max + 4 * m)
            found = 0;%在原先的S_c有没有找到雷同的i和j的标志
            for n = 1: numel(S_c(:,1))
                if (S_c(n,1)) == i && (S_c(n,2)) == j
                    found = 1;
                    break
                end
            end
            
            if found == 0
                k = k + 1;
                IC1_r_t= imtranslate(IC1_r,[i * imrf,j * imrf]);
                IC_t=IC1 | IC1_r_t;
                %imshow(IC_t);
                S_c(k,1)=i;
                S_c(k,2)=j;
                if sum(IC_t(:)) == N1 + N2 %说明图像完全分离
                    if Alg2==1
                        S_c(k,3)=(-1)*fun_XY(IC_t);
                    elseif Alg2==2
                        S_c(k,3)=(-1)*fun_YX(IC_t);
                    elseif Alg2==3
                        S_c(k,3)=(-1)*fun_bs2(IC_t,IC1)+(-1)*fun_Variance(IC_t)/100000000;
                    elseif Alg2==4
                        S_c(k,3)=(-1)*fun_Variance(IC_t)+(-1)*fun_bs2(IC_t,IC1)/100000000;
                    end
                else
                    S_c(k,3)=NaN;
                end
            end
        end
    end
    m = m + 1;
end

row_has_nan=any(isnan(S_c),2);
S_c(row_has_nan,:)=[];
[~, so] = sort(S_c(:,3), 'descend'); %so为降序后的序号排序

dx2_1=S_c(so(1),1);
dy2_1=S_c(so(1),2);

IC1_r=imtranslate(IC1_r,[dx2_1*imrf*1,dy2_1*imrf*1]);

% I=IC1 | IC1_r;
% imshow(I)

%分段第二次
while (dx2_1_last ~= dx2_1) || (dy2_1_last ~= dy2_1)
    
    dx2_1_last = dx2_1;
    dy2_1_last = dy2_1;
    
    S_c = zeros(1,3);
    k=0;
    for i=-i_max:i_max
        for j=-j_max:j_max
            k=k+1;
            IC1_r_t= imtranslate(IC1_r,[i*imrf*1,j*imrf*1]);
            IC_t=logical(IC1 | IC1_r_t);
            %             imshow(IC_t);
            S_c(k,1)=i;
            S_c(k,2)=j;
            
            if sum(IC_t(:)) == N1 + N2 %说明图像完全分离
                if Alg2==1
                    S_c(k,3)=(-1)*fun_XY(IC_t);
                elseif Alg2==2
                    S_c(k,3)=(-1)*fun_YX(IC_t);
                elseif Alg2==3
                    S_c(k,3)=(-1)*fun_bs2(IC_t,IC1)+(-1)*fun_Variance(IC_t)/100000000;
                elseif Alg2==4
                    S_c(k,3)=(-1)*fun_Variance(IC_t)+(-1)*fun_bs2(IC_t,IC1)/100000000;
                end
            else
                S_c(k,3)=NaN;
            end
        end
    end
    row_has_nan=any(isnan(S_c),2);
    S_c(row_has_nan,:)=[];
    [~, so] = sort(S_c(:,3), 'descend'); %so为降序后的序号排序
    
    dx2_1_t=S_c(so(1),1);
    dy2_1_t=S_c(so(1),2);
    
    IC1_r=imtranslate(IC1_r,[dx2_1_t * imrf * 1,dy2_1_t * imrf * 1]);
%     I=IC1 | IC1_r;
%     imshow(I)
    
    dx2_1 = dx2_1 + dx2_1_t;
    dy2_1 = dy2_1 + dy2_1_t;
end
clear S_c row_has_nan
%第2次寻找
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if found1 ==1
%     S_c = zeros(1,3);
%     k=0;
%     for i=-2:2
%         for j=-2:2
%             k=k+1;
%             IC1_r_t= imtranslate(IC1_r,[i*imrf,j*imrf]);
%             IC_t=IC1 | IC1_r_t;
%             %             imshow(IC_t);
%             S_c(k,1)=i;
%             S_c(k,2)=j;
%             if sum(IC_t(:))==2*N_s_c %说明图像完全分离
%                 if Alg2==0
%                     S_c(k,3)=(-1)*fun_bs(IC_t);
%                 elseif Alg2==1
%                     S_c(k,3)=(-1)*fun_bao2(IC_t);
%                 elseif Alg2==2
%                     S_c(k,3)=(-1)*fun_YX(IC_t);
%                 elseif Alg2==3
%                     S_c(k,3)=(-1)*fun_XY(IC_t);
%                 elseif Alg2==4
%                     S_c(k,3)=(-1)*fun_Variance(IC_t);
%                 elseif Alg2==5
%                     S_c(k,3)=(-1)*fun_S(IC_t);
%                 else
%                     S_c(k,3)=(-1)*fun_Lxy(IC_t);
%                 end
%             else
%                 S_c(k,3)=NaN;
%             end
%         end
%     end
%     row_has_nan=any(isnan(S_c),2);
%     S_c(row_has_nan,:)=[];
%
%     [~, so] = sort(S_c(:,3), 'descend'); %so为降序后的序号排序
%     dx2_2=S_c(so(1),1);
%     dy2_2=S_c(so(1),2);
%
%     IC1_r=imtranslate(IC1_r,[dx2_2*imrf,dy2_2*imrf]);
%
%     I=IC1 | IC1_r;
%     imshow(I)
%     clear S_c row_has_nan
% end

%第3次寻找
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_c = zeros(1,3);
k=0;
% i_max = min(round(l_x(I)/(8*imrf)),10);
% j_max = min(round(l_y(I)/(8*imrf)),10);
i_max = imrf;
j_max = imrf;
for i=-i_max:i_max
    for j=-j_max:j_max
        k=k+1;
        IC1_r_t= imtranslate(IC1_r,[i*1,j*1]);
        IC_t=IC1 | IC1_r_t;
        %imshow(IC_t);
        S_c(k,1)=i;
        S_c(k,2)=j;
        if sum(IC_t(:)) == N1 + N2 %说明图像完全分离
            if Alg2==1
                S_c(k,3)=(-1)*fun_XY(IC_t);
            elseif Alg2==2
                S_c(k,3)=(-1)*fun_YX(IC_t);
            elseif Alg2==3
                S_c(k,3)=(-1)*fun_bs2(IC_t,IC1)+(-1)*fun_Variance(IC_t)/100000000;
            elseif Alg2==4
                S_c(k,3)=(-1)*fun_Variance(IC_t)+(-1)*fun_bs2(IC_t,IC1)/100000000;
            end
            %             if Alg2==0
            %                 S_c(k,3)=(-1)*fun_bs(IC_t);
            %             elseif Alg2==1
            %                 S_c(k,3)=(-1)*fun_bao2(IC_t,IC1);
            %             elseif Alg2==2
            %                 S_c(k,3)=(-1)*fun_YX(IC_t);
            %             elseif Alg2==3
            %                 S_c(k,3)=(-1)*fun_XY(IC_t);
            %             elseif Alg2==4
            %                 S_c(k,3)=(-1)*fun_Variance(IC_t);
            %             elseif Alg2==5
            %                 S_c(k,3)=(-1)*fun_S(IC_t);
            %             else
            %                 S_c(k,3)=(-1)*fun_Lxy(IC_t);
            %             end
        else
            S_c(k,3)=NaN;
        end
    end
end

row_has_nan=any(isnan(S_c),2);
S_c(row_has_nan,:)=[];

% if isempty(S_c)%如果是空的，即没有找到
%     found3 = 0;
%     dx2_3 = 0;
%     dy2_3 = 0;
%
%     %     dx2_4 = 0;
%     %     dy2_4 = 0;
%
% end
row_has_nan=any(isnan(S_c),2);
S_c(row_has_nan,:)=[];
[~, so] = sort(S_c(:,3), 'descend'); %so为降序后的序号排序

dx2_3=S_c(so(1),1);
dy2_3=S_c(so(1),2);

IC1_r=imtranslate(IC1_r,[dx2_3*1,dy2_3*1]);
I=IC1 | IC1_r;
% imshow(I)
clear S_c row_has_nan

%第4次寻找
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if found3 ==1
%     S_c = zeros(1,3);
%     k=0;
%     for i=-2:2
%         for j=-2:2
%             k=k+1;
%             IC1_r_t = imtranslate(IC1_r,[i,j]);
%             IC_t=IC1 | IC1_r_t;
%             %             imshow(IC_t);
%             S_c(k,1)=i;
%             S_c(k,2)=j;
%             if sum(IC_t(:))==2*N_s_c %说明图像完全分离
%                 if Alg2==0
%                     S_c(k,3)=(-1)*fun_bs(IC_t);
%                 elseif Alg2==1
%                     S_c(k,3)=(-1)*fun_bao2(IC_t);
%                 elseif Alg2==2
%                     S_c(k,3)=(-1)*fun_YX(IC_t);
%                 elseif Alg2==3
%                     S_c(k,3)=(-1)*fun_XY(IC_t);
%                 elseif Alg2==4
%                     S_c(k,3)=(-1)*fun_Variance(IC_t);
%                 elseif Alg2==5
%                     S_c(k,3)=(-1)*fun_S(IC_t);
%                 else
%                     S_c(k,3)=(-1)*fun_Lxy(IC_t);
%                 end
%             else
%                 S_c(k,3)=NaN;
%             end
%         end
%     end
%     row_has_nan=any(isnan(S_c),2);
%     S_c(row_has_nan,:)=[];
%
%     [~, so] = sort(S_c(:,3), 'descend'); %so为降序后的序号排序
%     dx2_4=S_c(so(1),1);
%     dy2_4=S_c(so(1),2);

%     IC1_r=imtranslate(IC1_r,[dx2_4,dy2_4]);

% end
I = trans2mid(I,0);
I2_hr = double(I);

dx2_2 = 0;
dy2_2 = 0;
dx2_4 = 0;
dy2_4 = 0;
% dx222=round(dx2*imrf)/imrf+(dx2_1*1+dx2_2)+(dx2_3*1+dx2_4)/imrf;
% dy222=round(dy2*imrf)/imrf+(dy2_1*1+dy2_2)+(dy2_3*1+dy2_4)/imrf;
dx222 = dx2+(dx2_1*1+dx2_2)+(dx2_3*1+dx2_4)/imrf;
dy222 = dy2+(dy2_1*1+dy2_2)+(dy2_3*1+dy2_4)/imrf;
end

