function [dx333,dy333,I3_hr] = high_accuracy3(I3_hr,dx3,dy3,sita0,imrf,ijmax,RF)
%HIGH_ACCURACY3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
dx3_1_last = 0;
dy3_1_last = 0;

if ijmax < 1/RF
    ijmax = 1/RF;
end
I3_hr = logical(I3_hr);
% N_s_c = sum(I3_hr(:));
% found1 = 1;
% found3 = 1;
I3_hr = trans2mid(I3_hr,1);
I3_hr = logical(I3_hr);
%��1��Ѱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dx3>=0
    
    IC1 = trans2leftbottom(I3_hr);
    IC1_r = imrotate(I3_hr,180);
    IC1_r = trans2leftbottom(IC1_r);
    
    IC1=imtranslate(IC1,[12*2*imrf+1,-12*2*imrf-1]);
    IC1_r=imtranslate(IC1_r,[12*2*imrf+1,-12*2*imrf-1]);
    
    N1=sum(IC1(:));
    N2=sum(IC1_r(:));
    
else
    
    IC1 = trans2rightbottom(I3_hr);
    IC1_r = imrotate(I3_hr,180);
    IC1_r = trans2rightbottom(IC1_r);
    %     IC1_r = trans2rightbottom(I3_hr);
    
    IC1=imtranslate(IC1,[-12*2*imrf-1,-12*2*imrf-1]);
    IC1_r=imtranslate(IC1_r,[-12*2*imrf-1,-12*2*imrf-1]);
    
    N1=sum(IC1(:));
    N2=sum(IC1_r(:));
end

if dy3 > 0
    IC1 = imtranslate(IC1,[0,-dy3*imrf]);
    IC1_r = imtranslate(IC1_r,[0,-dy3*imrf]);
end

IC1_r=imtranslate(IC1_r,[round(dx3*imrf),round(dy3*imrf)]);

while sum(IC1_r(:)) < sum(IC1(:))
    
    I_base = zeros(2*size(IC1));
    IC1 = paste_img(I_base,IC1,1,1);
    
    if dx3>=0
        
        IC1 = trans2leftbottom(IC1);
        IC1_r = imrotate(IC1,180);
        IC1_r = trans2leftbottom(IC1_r);
        
        IC1=imtranslate(IC1,[12*2*imrf+1,-12*2*imrf-1]);
        IC1_r=imtranslate(IC1_r,[12*2*imrf+1,-12*2*imrf-1]);
        
        N1=sum(IC1(:));
        N2=sum(IC1_r(:));
        
    else
        
        IC1 = trans2rightbottom(IC1);
        IC1_r = imrotate(IC1,180);
        IC1_r = trans2rightbottom(IC1_r);
        
        IC1=imtranslate(IC1,[-12*2*imrf-1,-12*2*imrf-1]);
        IC1_r=imtranslate(IC1_r,[-12*2*imrf-1,-12*2*imrf-1]);
        
        N1=sum(IC1(:));
        N2=sum(IC1_r(:));
    end
    
    if dy3 > 0
        IC1 = imtranslate(IC1,[0,-dy3*imrf]);
        IC1_r = imtranslate(IC1_r,[0,-dy3*imrf]);
    end
    
    IC1_r=imtranslate(IC1_r,[round(dx3*imrf),round(dy3*imrf)]);
    
end

S_c = zeros(1,3);
k=0;

I=IC1 | IC1_r;
%  imshow(I)
i_max = ceil(min(round(l_x(I)/(3*imrf)),ijmax)/2);
j_max = ceil(min(round(l_y(I)/(3*imrf)),ijmax)/2);

%�ֶε�һ��
for i = -i_max:i_max
    for j = -j_max:j_max
        k = k+1;
        IC1_r_t = imtranslate(IC1_r,[i*imrf*2,j*imrf*2]);
        IC_t = IC1 | IC1_r_t;
        %         imshow(IC_t);
        S_c(k,1) = i;
        S_c(k,2) = j;
        
        %IC_t = trans2mid(IC_t);
        %imshow(IC_t)
        if sum(IC_t(:)) == N1 + N2 %˵��ͼ����ȫ����
            S_c(k,3) = (-1) * fun_XY_ro(IC_t,sita0);
        else
            S_c(k,3) = NaN;
        end
    end
end

m = 1;
while all(isnan(S_c(:,3)))%�жϵ������ǲ���ȫ��NaN
    %���Ƚ������Ŵ�
    I_base = zeros(numel(IC1(:,1)) + 8 * imrf,numel(IC1(1,:)) + 8 * imrf);
    IC1 = paste_img(I_base, IC1, 4 * imrf + 1, 4 * imrf + 1); IC1 = logical(IC1);
    IC1_r = paste_img(I_base, IC1_r, 4 * imrf + 1, 4 * imrf + 1); IC1_r = logical(IC1_r);
    
    for i = (-i_max - 4 * m) : (i_max + 4 * m)
        for j = (-j_max - 4 * m) : (j_max + 4 * m)
            found = 0;%��ԭ�ȵ�S_c��û���ҵ���ͬ��i��j�ı�־
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
                %               imshow(IC_t);
                S_c(k,1)=i;
                S_c(k,2)=j;
                if sum(IC_t(:)) == N1 + N2 %˵��ͼ����ȫ����
                    S_c(k,3)=(-1) * fun_XY_ro(IC_t,sita0);
                else
                    S_c(k,3)=NaN;
                end
            end
        end
    end
    m = m + 1;
end


row_has_nan=any(isnan(S_c),2);
S_c(row_has_nan,:) = [];
[~, so] = sort(S_c(:,3), 'descend'); %soΪ�������������

dx3_1 = S_c(so(1),1);
dy3_1 = S_c(so(1),2);

IC1_r=imtranslate(IC1_r,[dx3_1*imrf*1,dy3_1*imrf*1]);

% I=IC1 | IC1_r;
% imshow(I)

%�ֶεڶ���
while (dx3_1_last ~= dx3_1) || (dy3_1_last ~= dy3_1)
    
    dx3_1_last = dx3_1;
    dy3_1_last = dy3_1;
    
    S_c = zeros(1,3);
    k=0;
    for i=-i_max:i_max
        for j=-j_max:j_max
            k=k+1;
            IC1_r_t= imtranslate(IC1_r,[i*imrf*1,j*imrf*1]);
            IC_t=IC1 | IC1_r_t;
            %         imshow(IC_t);
            S_c(k,1)=i;
            S_c(k,2)=j;
            
            %IC_t = trans2mid(IC_t);
            %imshow(IC_t)
            if sum(IC_t(:)) == N1 + N2 %˵��ͼ����ȫ����
                S_c(k,3)=(-1)*fun_XY_ro(IC_t,sita0);
            else
                S_c(k,3)=NaN;
            end
        end
    end
    
    row_has_nan=any(isnan(S_c),2);
    S_c(row_has_nan,:)=[];
    [~, so] = sort(S_c(:,3), 'descend'); %soΪ�������������
    
    dx3_1_t=S_c(so(1),1);
    dy3_1_t=S_c(so(1),2);
    
    IC1_r=imtranslate(IC1_r,[dx3_1_t * imrf * 1,dy3_1_t * imrf * 1]);
%     I=IC1 | IC1_r;
%     imshow(I)
    
    dx3_1 = dx3_1 + dx3_1_t;
    dy3_1 = dy3_1 + dy3_1_t;
end
clear S_c row_has_nan
%��2��Ѱ��
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
%
%             %IC_t = trans2mid(IC_t);
%             %imshow(IC_t)
%             if sum(IC_t(:))==2*N_s_c %˵��ͼ����ȫ����
%                 S_c(k,3)=(-1)*fun_XY_ro(IC_t,sita0);
%             else
%                 S_c(k,3)=NaN;
%             end
%         end
%     end
%     row_has_nan=any(isnan(S_c),2);
%     S_c(row_has_nan,:)=[];
%
%     [~, so] = sort(S_c(:,3), 'descend'); %soΪ�������������
%     dx3_2=S_c(so(1),1);
%     dy3_2=S_c(so(1),2);
%
%     IC1_r=imtranslate(IC1_r,[dx3_2*imrf,dy3_2*imrf]);
%
%     I=IC1 | IC1_r;
%     imshow(I)
%     clear S_c row_has_nan
% end

%��3��Ѱ��
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
        %         imshow(IC_t);
        S_c(k,1)=i;
        S_c(k,2)=j;
        
        %IC_t = trans2mid(IC_t);
        %imshow(IC_t)
        if sum(IC_t(:)) == N1 + N2 %˵��ͼ����ȫ����
            S_c(k,3)=(-1)*fun_XY_ro(IC_t,sita0);
        else
            S_c(k,3)=NaN;
        end
    end
end

row_has_nan=any(isnan(S_c),2);
S_c(row_has_nan,:)=[];
[~, so] = sort(S_c(:,3), 'descend'); %soΪ�������������

dx3_3=S_c(so(1),1);
dy3_3=S_c(so(1),2);

IC1_r=imtranslate(IC1_r,[dx3_3*1,dy3_3*1]);
I=IC1 | IC1_r;
% imshow(I)
clear S_c row_has_nan

%��4��Ѱ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if found3 ==1
%     S_c = zeros(1,3);
%     k=0;
%     for i=-2:2
%         for j=-2:2
%             k=k+1;
%             IC1_r_t= imtranslate(IC1_r,[i,j]);
%             IC_t=IC1 | IC1_r_t;
%             %             imshow(IC_t);
%             S_c(k,1)=i;
%             S_c(k,2)=j;
%
%             %         IC_t = trans2mid(IC_t);
%             %         imshow(IC_t)
%             if sum(IC_t(:))==2*N_s_c %˵��ͼ����ȫ����
%                 S_c(k,3)=(-1)*fun_XY_ro(IC_t,sita0);
%             else
%                 S_c(k,3)=NaN;
%             end
%         end
%     end
%     row_has_nan=any(isnan(S_c),2);
%     S_c(row_has_nan,:)=[];
%
%     [~, so] = sort(S_c(:,3), 'descend'); %soΪ�������������
%     dx3_4=S_c(so(1),1);
%     dy3_4=S_c(so(1),2);
%
%     IC1_r=imtranslate(IC1_r,[dx3_4,dy3_4]);

% end

I = trans2mid(I,0);
I3_hr = double(I);

dx3_2 = 0;
dy3_2 = 0;
dx3_4 = 0;
dy3_4 = 0;
% dx333=round(dx3*imrf)/imrf+(dx3_1*1+dx3_2)+(dx3_3*1+dx3_4)/imrf;
% dy333=round(dy3*imrf)/imrf+(dy3_1*1+dy3_2)+(dy3_3*1+dy3_4)/imrf;
dx333 = dx3+(dx3_1*1+dx3_2)+(dx3_3*1+dx3_4)/imrf;
dy333 = dy3+(dy3_1*1+dy3_2)+(dy3_3*1+dy3_4)/imrf;
end

