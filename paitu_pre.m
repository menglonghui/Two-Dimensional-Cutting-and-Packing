function [coordinate1,coordinate2,N] = paitu_pre(ZP,dx1,dx2,dx3,dy1,dy2,dy3,sita1,N_dx3,N_dx3_max,N_shu,I0,lx_single,ly_single,x_py,y_py,sita_xz) %#ok<*INUSL>
%PAITU 此处显示有关此函数的摘要
%   此处显示详细说明
%I_mix_tt用来铺零件，该零件不与自己做比较判断，只与之前排好的零件做比较,之前排好零件的板子就是I_ref
% I=I_s_m0;
N_shu = N_shu + 1;
N_dx3 = N_dx3 + 1;
N_dx3_max = N_dx3_max + 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%先进行缩小
N = 0;
% dy3 = dy3 / BS;

% I0 = imresize(I0,1/BS);
% I_s_m0 = imresize(I_s_m0,1/BS);
% I_s_m_r0 = imresize(I_s_m_r0,1/BS);
% I_mix = imresize(I_mix,1/BS);
% I_mix_tt = imresize(I_mix_tt,1/BS);

% I_s_m0(I_s_m0 < 0.55)=0;
% I_s_m_r0(I_s_m_r0 < 0.55)=0;
% I_mix(I_mix < 0.55)=0;
% I_mix_tt(I_mix_tt < 0.55)=0;

% I_s_m0 = logical(I_s_m0);
% I_s_m_r0 = logical(I_s_m_r0);

% lx_single = lx_single / BS;
% ly_single = ly_single / BS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I_s_m0 = logical(I_s_m0);
% I_s_m_r0 = logical(I_s_m_r0);

% imshow(I_ref);
%l_x_initial = l_x(I_ref)-lx_single;
mm = 0;
nn = 0;
coordinate1 = [];
coordinate2 = [];
% if ZP == 1 %左偏
extend_max = 5;

if dx1 <= 0 && dy1 <= 0
    leftdoor = - abs(dx1);
    leftdoor_r = 0;
    rightdoor = numel(I0(1,:)) - lx_single - abs(dx1);
    rightdoor_r = numel(I0(1,:)) - lx_single;
    bottomdoor = 0;
    bottomdoor_r = abs(dy1);
    topdoor = - numel(I0(:,1)) + ly_single;
    topdoor_r = - numel(I0(:,1)) + ly_single + abs(dy1);
elseif dx1 <= 0 && dy1 > 0
    leftdoor = - abs(dx1);
    leftdoor_r = 0;
    rightdoor = numel(I0(1,:)) - lx_single - abs(dx1);
    rightdoor_r = numel(I0(1,:)) - lx_single;
    bottomdoor = abs(dy1);
    bottomdoor_r = 0;
    topdoor = - numel(I0(:,1)) + ly_single + abs(dy1);
    topdoor_r = - numel(I0(:,1)) + ly_single;
elseif dx1 > 0 && dy1 <= 0
    leftdoor = 0;
    leftdoor_r = - abs(dx1);
    rightdoor = numel(I0(1,:)) - lx_single;
    rightdoor_r = numel(I0(1,:)) - lx_single - abs(dx1);
    bottomdoor = 0;
    bottomdoor_r = abs(dy1);
    topdoor = - numel(I0(:,1)) + ly_single;
    topdoor_r = - numel(I0(:,1)) + ly_single + abs(dy1);
else % dx1>0 && dy1>0
    leftdoor = 0;
    leftdoor_r = - abs(dx1);
    rightdoor = numel(I0(1,:)) - lx_single;
    rightdoor_r = numel(I0(1,:)) - lx_single - abs(dx1);
    bottomdoor = abs(dy1);
    bottomdoor_r = 0;
    topdoor = - numel(I0(:,1)) + ly_single + abs(dy1);
    topdoor_r = - numel(I0(:,1)) + ly_single;
end

if ZP == 1
    i = - extend_max;
    i_max = N_dx3_max + extend_max;
else
    i = N_dx3 - N_dx3_max - extend_max;
    i_max = N_dx3 + extend_max;
end

while i <= i_max% && N_flag == 0 %i<=N_dx3_max
    
    i = i + 1;
    
    for j = 0 : N_shu  %i为dx3，dy3的增加次数，而j为dx2，dy2的增加次数
        
        dx = j * dx2 + (i - 1) * dx3 + x_py;
        dy = j * dy2 + (i - 1) * dy3 + y_py;
        
        if dx >= leftdoor && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
            N=N+1;
            mm=mm+1;
            coordinate1(mm,:)=[dx, dy, sita1, dx, dy]; %#ok<*AGROW>
        end
        
        if dx >= leftdoor_r && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
            N=N+1;
            nn=nn+1;
            coordinate2(nn,:)=[dx, dy, sita1 + 180, dx, dy];
        end
    end
end
if isempty(coordinate1) == 0
    coordinate1(:,3) = coordinate1(:,3) + sita_xz;
end
if isempty(coordinate2) == 0
    coordinate2(:,3) = coordinate2(:,3) + sita_xz;
end
end


