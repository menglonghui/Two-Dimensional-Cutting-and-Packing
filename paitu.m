function [coordinate1,coordinate2,I_mix,I_mix_tt,N_flag,N] = paitu(flag_N_h,ZP,n_p,dx1,dx2,dx3,dy1,dy2,dy3,sita1,Nmax,key,...
    N_flag,N_dx3,N_dx3_max,N_shu,I0,bujin,I_mix,~,I_mix_tt,N,I_s_m0,I_s_m_r0,~,~,lx_single,ly_single,single,n_w,p,N_part,after_try,lx_add) %#ok<*INUSL>
%PAITU �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%I_mix_tt���������������������Լ����Ƚ��жϣ�ֻ��֮ǰ�źõ�������Ƚ�,֮ǰ�ź�����İ��Ӿ���I_ref
% I=I_s_m0;
N_shu = round(N_shu * 1.2);
BS = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Ƚ�����С
dx1 = dx1 / BS;
dx2 = dx2 / BS;
dx3 = dx3 / BS;

dy1 = dy1 / BS;
dy2 = dy2 / BS;
% dy3 = dy3 / BS;

% I0 = imresize(I0,1/BS);
% I_s_m0 = imresize(I_s_m0,1/BS);
% I_s_m_r0 = imresize(I_s_m_r0,1/BS);
% I_mix = imresize(I_mix,1/BS);
% I_mix_tt = imresize(I_mix_tt,1/BS);

I_s_m0(I_s_m0 < 0.55)=0;
I_s_m_r0(I_s_m_r0 < 0.55)=0;
I_mix(I_mix < 0.55)=0;
I_mix_tt(I_mix_tt < 0.55)=0;

I_s_m0 = logical(I_s_m0);
I_s_m_r0 = logical(I_s_m_r0);

lx_single = lx_single / BS;
ly_single = ly_single / BS;
bujin = bujin / BS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_ref = I_mix;

I_ref = logical(I_ref);
% I_s_m0 = logical(I_s_m0);
% I_s_m_r0 = logical(I_s_m_r0);

N_I_ref = sum(I_ref(:));
N_single = sum(I_s_m0(:));
N_single_r = sum(I_s_m_r0(:));
% imshow(I_ref);
%l_x_initial = l_x(I_ref)-lx_single;
mm = 0;
nn = 0;
% coordinate1=zeros(1,3);
% coordinate2=zeros(1,3);
coordinate1 = [];
coordinate2 = [];
% if ZP == 1 %��ƫ
extend_max = max(round(N_dx3_max * 0.1),4);

if after_try == 0 % || n_w == 1
    
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
    
    if n_p == 1%��һ����ͼ
        
        if ZP == 1
            i = 0;
            i_max = N_dx3_max + extend_max;
        else
            i = N_dx3 - N_dx3_max - extend_max;
            i_max = N_dx3 + extend_max;
        end
        
        while i <= i_max && N_flag == 0 %i<=N_dx3_max
            i = i + 1;
            for j = 0 : N_shu  %iΪdx3��dy3�����Ӵ�������jΪdx2��dy2�����Ӵ���
                dx = j * dx2 + (i - 1) * dx3;
                dy = j * dy2;% + (i - 1) * dy3;
                if (dx <= numel(I0(1,:))) && (dx >= -abs(dx1)) && dy >= -numel(I0(:,1)) && dy <= abs(dy1)
                    
                    if dy1 <= 0
                        
                        if dx >= leftdoor && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m | I_mix;
                            if single == 1 || (p == N_part && n_w == 1) %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)]; %#ok<*AGROW>
                                %I_mix=I_tt;
                            else
                                I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix_temp=I_tt;
                                else
                                    I_tt=I_s_m | I_ref;%��һ����Ϊ�˺����ж��ǲ����ص���׼����
                                    %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                    if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        mm=mm+1;
                                        coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        if dx >= leftdoor_r && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m_r | I_mix;
                            if single == 1 || (p == N_part && n_w == 1)   %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m_r | I_ref;
                                    %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                    if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        nn=nn+1;
                                        coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m_r;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        
                    else
                        
                        if dx >= leftdoor_r && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m_r | I_mix;
                            if single == 1 || (p == N_part && n_w == 1)   %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m_r | I_ref;
                                    %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                    if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        nn=nn+1;
                                        coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m_r;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        
                        if dx >= leftdoor && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m | I_mix;
                            if single == 1 || (p == N_part && n_w == 1) %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)]; %#ok<*AGROW>
                                %I_mix=I_tt;
                            else
                                I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix_temp=I_tt;
                                else
                                    I_tt=I_s_m | I_ref;%��һ����Ϊ�˺����ж��ǲ����ص���׼����
                                    %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                    if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        mm=mm+1;
                                        coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    else%���ǵ�һ����ͼ
        
        if ZP == 1
            i = 0;
            i_max = N_dx3_max + extend_max;
        else
            i = N_dx3 - N_dx3_max - extend_max;
            i_max = N_dx3 + extend_max;
        end
        
        while i <= i_max && N_flag == 0 %i<=N_dx3_max
            i = i + 1;
            for j = 0 : N_shu  %iΪdx3��dy3�����Ӵ�������jΪdx2��dy2�����Ӵ���
                dx = j * dx2 + (i - 1) * dx3;
                dy = j * dy2;% + (i - 1) * dy3;
                %             dy1=-dy;
                %             disp([num2str(dx),',',num2str(dy1)]);
                if (abs(dx) >= (numel(I0(1,:))-bujin-abs(dx1)-lx_single)-100) && (dx <= (numel(I0(1,:))) && dy <= abs(dy1)) && (dy >= -numel(I0(:,1)))%�ж�ͼ���ǲ������ڴ�ŵ���ȷ�ķ�Χ��
                    
                    if dy1 <= 0
                        if dx > (rightdoor - bujin) && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         rightdoor - bujin
                            %                         rightdoor
                            %                     if dx > (rightdoor - bujin - (lx_single / 2)) && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m | I_mix;
                            if (single == 1 || (p == N_part && n_w == 1)) %|| dx > l_x_initial  %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m | I_ref;
                                    %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                    if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        mm=mm+1;
                                        coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        
                        if dx > (rightdoor_r - bujin) && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         rightdoor_r - bujin
                            %                         rightdoor_r
                            %                     if dx > (rightdoor_r - bujin - (lx_single / 2)) && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m_r | I_mix;
                            if (single == 1 || (p == N_part && n_w == 1)) %|| dx > l_x_initial %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m_r | I_ref;
                                    %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                    if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        nn=nn+1;
                                        coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m_r;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        
                    else
                        
                        if dx > (rightdoor_r - bujin) && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         rightdoor_r - bujin
                            %                         rightdoor_r
                            %                     if dx > (rightdoor_r - bujin - (lx_single / 2)) && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                            %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m_r | I_mix;
                            if (single == 1 || (p == N_part && n_w == 1)) %|| dx > l_x_initial %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m_r | I_ref;
                                    %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                    if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        nn=nn+1;
                                        coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m_r;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                        
                        if dx > (rightdoor - bujin) && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         rightdoor - bujin
                            %                         rightdoor
                            %                     if dx > (rightdoor - bujin - (lx_single / 2)) && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                            %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            %                         I_tt=I_s_m | I_mix;
                            if (single == 1 || (p == N_part && n_w == 1)) %|| dx > l_x_initial  %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                %I_mix=I_tt;
                            else
                                I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                                if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                else
                                    I_tt=I_s_m | I_ref;
                                    %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                    if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                        %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                        N=N+1;
                                        if eval(['N>Nmax.' key])
                                            N=N-1;
                                            N_flag=1;
                                            break
                                        end
                                        mm=mm+1;
                                        coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                        I_mix_tt=I_mix_tt | I_s_m;
                                        %imshow(I_mix_tt)
                                        %I_mix=I_tt;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
else %ȡ����after_try��ֵ����after_try��ֵΪ1�Ļ���������ߵĻ��Ų��������������ʼִ���ұߵĶ�����Ĳ��֣�����ӵ�lx_add�ĺ��򲿷�
    
    if dx1 <= 0 && dy1 <= 0
        leftdoor = - abs(dx1) + numel(I0(1,:)) - lx_add - lx_single;
        leftdoor_r = 0 + numel(I0(1,:)) - lx_add - lx_single;
        rightdoor = numel(I0(1,:)) - lx_single - abs(dx1);
        rightdoor_r = numel(I0(1,:)) - lx_single;
        bottomdoor = 0;
        bottomdoor_r = abs(dy1);
        topdoor = - numel(I0(:,1)) + ly_single;
        topdoor_r = - numel(I0(:,1)) + ly_single + abs(dy1);
    elseif dx1 <= 0 && dy1 > 0
        leftdoor = - abs(dx1) + numel(I0(1,:)) - lx_add - lx_single;
        leftdoor_r = 0 + numel(I0(1,:)) - lx_add - lx_single;
        rightdoor = numel(I0(1,:)) - lx_single - abs(dx1);
        rightdoor_r = numel(I0(1,:)) - lx_single;
        bottomdoor = abs(dy1);
        bottomdoor_r = 0;
        topdoor = - numel(I0(:,1)) + ly_single + abs(dy1);
        topdoor_r = - numel(I0(:,1)) + ly_single;
    elseif dx1 > 0 && dy1 <= 0
        leftdoor = 0 + numel(I0(1,:)) - lx_add - lx_single;
        leftdoor_r = - abs(dx1) + numel(I0(1,:)) - lx_add - lx_single;
        rightdoor = numel(I0(1,:)) - lx_single;
        rightdoor_r = numel(I0(1,:)) - lx_single - abs(dx1);
        bottomdoor = 0;
        bottomdoor_r = abs(dy1);
        topdoor = - numel(I0(:,1)) + ly_single;
        topdoor_r = - numel(I0(:,1)) + ly_single + abs(dy1);
    else % dx1>0 && dy1>0
        leftdoor = 0 + numel(I0(1,:)) - lx_add - lx_single;
        leftdoor_r = - abs(dx1) + numel(I0(1,:)) - lx_add - lx_single;
        rightdoor = numel(I0(1,:)) - lx_single;
        rightdoor_r = numel(I0(1,:)) - lx_single - abs(dx1);
        bottomdoor = abs(dy1);
        bottomdoor_r = 0;
        topdoor = - numel(I0(:,1)) + ly_single + abs(dy1);
        topdoor_r = - numel(I0(:,1)) + ly_single;
    end
    
    
    if ZP == 1
        i = 0;
        i_max = N_dx3_max + extend_max;
    else
        i = N_dx3 - N_dx3_max - extend_max;
        i_max = N_dx3 + extend_max;
    end
    
    while i <= i_max && N_flag == 0 %i<=N_dx3_max
        i = i + 1;
        for j = 0 : N_shu  %iΪdx3��dy3�����Ӵ�������jΪdx2��dy2�����Ӵ���
            dx = j * dx2 + (i - 1) * dx3;
            dy = j * dy2;% + (i - 1) * dy3;
            if (dx <= numel(I0(1,:))) && (dx >= -abs(dx1)) && dy >= -numel(I0(:,1)) && dy <= abs(dy1)
                
                if dy1 <= 0
                    
                    if dx > leftdoor && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                        %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                        %                         I_tt=I_s_m | I_mix;
                        if single == 1 || (p == N_part && n_w == 1) %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                            N=N+1;
                            if eval(['N>Nmax.' key])
                                N=N-1;
                                N_flag=1;
                                break
                            end
                            mm=mm+1;
                            coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)]; %#ok<*AGROW>
                            %I_mix=I_tt;
                        else
                            I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                I_mix_tt=I_mix_tt | I_s_m;
                                %imshow(I_mix_tt)
                                %I_mix_temp=I_tt;
                            else
                                I_tt=I_s_m | I_ref;%��һ����Ϊ�˺����ж��ǲ����ص���׼����
                                %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                    %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                end
                            end
                        end
                    end
                    if dx > leftdoor_r && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                        %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                        %                         I_tt=I_s_m_r | I_mix;
                        if single == 1 || (p == N_part && n_w == 1)   %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                            N=N+1;
                            if eval(['N>Nmax.' key])
                                N=N-1;
                                N_flag=1;
                                break
                            end
                            nn=nn+1;
                            coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                            %I_mix=I_tt;
                        else
                            I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            if n_w == 1%imrf_special
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                I_mix_tt=I_mix_tt | I_s_m_r;
                                %imshow(I_mix_tt)
                                %I_mix=I_tt;
                            else
                                I_tt=I_s_m_r | I_ref;
                                %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                    %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                end
                            end
                        end
                    end
                    
                else
                    
                    if dx > leftdoor_r && dx <= rightdoor_r && dy <= bottomdoor_r &&  dy >= topdoor_r
                        %                         I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                        %                         I_tt=I_s_m_r | I_mix;
                        if single == 1 || (p == N_part && n_w == 1)   %&& n_w>1 %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                            N=N+1;
                            if eval(['N>Nmax.' key])
                                N=N-1;
                                N_flag=1;
                                break
                            end
                            nn=nn+1;
                            coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                            %I_mix=I_tt;
                        else
                            I_s_m_r=imtranslate(I_s_m_r0,[round(dx),round(dy)]);
                            if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                nn=nn+1;
                                coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                I_mix_tt=I_mix_tt | I_s_m_r;
                                %imshow(I_mix_tt)
                                %I_mix=I_tt;
                            else
                                I_tt=I_s_m_r | I_ref;
                                %d_pix = abs(sum(I_s_m_r(:))-sum(I_s_m_r0(:)));
                                if abs(sum(I_tt(:)) - N_I_ref - N_single_r) <= 1 %|| abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m_r0(:)))==d_pix
                                    %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    nn=nn+1;
                                    coordinate2(nn,:)=[BS*dx,BS*dy,mod(180+sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m_r;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                end
                            end
                        end
                    end
                    
                    if dx > leftdoor && dx <= rightdoor && dy <= bottomdoor &&  dy >= topdoor
                        %                         I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                        %                         I_tt=I_s_m | I_mix;
                        if single == 1 || (p == N_part && n_w == 1) %ֻ����������Ų���Ҫ�����Ƿ�����������ص�������֤
                            N=N+1;
                            if eval(['N>Nmax.' key])
                                N=N-1;
                                N_flag=1;
                                break
                            end
                            mm=mm+1;
                            coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)]; %#ok<*AGROW>
                            %I_mix=I_tt;
                        else
                            I_s_m=imtranslate(I_s_m0,[round(dx),round(dy)]);
                            if n_w == 1%���ŵ�ʱ������������̵��ǵ�һ�������Ҳ����Ҫ�����Ƿ�����������ص�������֤
                                N=N+1;
                                if eval(['N>Nmax.' key])
                                    N=N-1;
                                    N_flag=1;
                                    break
                                end
                                mm=mm+1;
                                coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                I_mix_tt=I_mix_tt | I_s_m;
                                %imshow(I_mix_tt)
                                %I_mix_temp=I_tt;
                            else
                                I_tt=I_s_m | I_ref;%��һ����Ϊ�˺����ж��ǲ����ص���׼����
                                %d_pix = abs(sum(I_s_m(:))-sum(I_s_m0(:)));
                                if (abs(sum(I_tt(:)) - N_I_ref - N_single) <= 1) %|| (abs(sum(I_tt(:))-sum(I_ref(:)) - sum(I_s_m0(:)))==d_pix)
                                    %˵���ƶ����ͼ����ȫ���ڰ��ڣ����Ҳ�δ����������ص�
                                    N=N+1;
                                    if eval(['N>Nmax.' key])
                                        N=N-1;
                                        N_flag=1;
                                        break
                                    end
                                    mm=mm+1;
                                    coordinate1(mm,:)=[BS*dx,BS*dy,mod(sita1,360)];
                                    I_mix_tt=I_mix_tt | I_s_m;
                                    %imshow(I_mix_tt)
                                    %I_mix=I_tt;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% I_mix = I_mix | I_mix_tt;
% I_mix = imresize(I_mix,BS);
% I_mix_tt = imresize(I_mix_tt,BS);

if (eval(['Nmax.' key ' == N']) || flag_N_h == 1) && single == 0 %��������˻��߰��������,�����ǻ��ţ����Ų���
    I_mix=I_mix | I_mix_tt;
end

end
