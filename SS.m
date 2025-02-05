clear
clc
for ppp=1:1
    disp(['ppp=',num2str(ppp)]);
% %            json_path_input = 'D:\input\0331.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\axuediche-2.dwg_input.json';
            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part5.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\ceshi2.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part5.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part2.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part5.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\axuediche-2.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\ceshi2.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part1.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\ceshi1_25.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\part1.dwg_input.json';
% %            json_path_input = 'C:\Users\Administrator\AppData\Roaming\WITSOFT\WIT 3D\WIT 3D 2023\output\temp\Part4_1.dwg_input.json';
    
    
    total_time = tic;
    % feature('DefaultCharacterSet', 'UTF-8');
    delete(gcf)
    jsonStr = fileread(json_path_input);% 读取JSON文件
    data_list = jsondecode(jsonStr);% 解析JSON文件
    data_list_o = data_list;
    N_part = numel(fieldnames(data_list.part));
    time = data_list.param.time;
    enablemix = data_list.param.enablemix;
    toolszie = data_list.param.toolszie;
    ignorehole = data_list.param.ignorehole;
    partgap = data_list.param.partgap;
    leftmargin = data_list.param.leftmargin;
    rightmargin = data_list.param.rightmargin;
    topmargin = data_list.param.topmargin;
    downmargin = data_list.param.downmargin;
    sheetsize = data_list.param.specifysheet;
    % enablenesting = data_list.param.enablenesting;
    enablerevolve = data_list.param.enablerevolve;
    accuracy = data_list.param.accuracy;
    N_times = data_list.param.multiple_of_population;
    algorithmselect = data_list.param.algorithmselect;
    outputpath = data_list.param.resultjsonpath;
    outputpath = strsplit(outputpath, '\');
    output_json_name = outputpath{end};
    outputpath(end) = [];
    outputpath = strjoin(outputpath, '\');
    sheet_size = str2double(strsplit(sheetsize, 'x'));
    sheet_X = sheet_size(1);
    sheet_Y = sheet_size(2);
    num_forname = num2str((str2double(regexp(output_json_name, '\d+', 'match'))));
    num_forname = strrep(num_forname,' ','_');
    % 获取结构体变量中所有的键值
    keys = fieldnames(data_list.part);
    s_filepath = struct();
    L_sz = zeros(N_part,1);
    
    % 遍历所有的键值
    for i = 1:length(keys)
        key = keys{i};
        eval(['Path.' key ' =  data_list.part.' key '.filepath;'])
        eval(['X_length.' key ' =  data_list.part.' key '.length;'])
        eval(['Y_length.' key ' =  data_list.part.' key '.width;'])
        eval(['Nmax.' key ' =  data_list.part.' key '.quantity;'])
        eval(['img.' key ' = {rgb2gray(imread(Path.' key '))};'])
        eval(['L_total.' key ' = cell2mat({X_length.' key '+Y_length.' key '});'])
        eval(['L_sz(i) = L_total.' key ';'])
        %         eval(['alg1.' key ' =  data_list.part.' key '.alg1;'])
        %         eval(['alg2.' key ' =  data_list.part.' key '.alg2;'])
        eval(['kb.' key ' = data_list.part.' key '.kb;'])
    end
    
    [L_sz, so] = sort(L_sz, 'descend');
    key_seq = cell(N_part,1);
    keys_t = keys;
    for i = 1:N_part
        for j = 1:length(keys_t)
            key = keys_t{j};
            if eval(['L_sz(i) == L_total.' key])
                key_seq{i} = key;
                keys_t(j) = [];
                break
            end
        end
    end
    
    N_linjie = 21;
    n = 2;%图案的数量
    nVar = 25*n; % x的长度
    nPop0 = 20;
    nPop = 20; % 种群规模大小
    maxIt = 3; % 迭代次数
    maxIt_va = 1;%有的时候怕陷入局部最优，因此多试几次，能找到最优值
    nPc = 1; % 交叉的比例
    imrf = 1/accuracy;%缩放比例
    N_multiple = 1;
    n_w = 1;%n_w表示混排的时候，同一个板子上的第几种零件，将其初始化为1
    N_white_pix = 0;
    
    I0 = zeros(sheet_Y - leftmargin - rightmargin + partgap, sheet_X - topmargin - downmargin + partgap);%空白板
    heiban = logical(I0);
    I0 = logical(I0);
    mm_pix = (sheet_X - topmargin - downmargin) / numel(I0(1,:));%每个像素代表几mm
    I_mix = heiban;
    lx_heiban = numel(heiban(1,:));
    ly_heiban = numel(heiban(:,1));
    
    for p = 1 : N_part %所排的图的个数
        
        loop_time = tic;
        lx_add = 0;
        danpai_exist = 0;%单排是否存在
        wp_new = 1;%是不是新的图案 workpiece_new
        find_best_dxy5_used = 0;
        if exist('N_pre','var')%如果该值存在
            clear N_pre
        end
        
        if exist('coordinate_5','var')%如果该值存在
            clear coordinate_5
        end
        
        if exist('coordinate_sbs','var')%如果该值存在
            clear coordinate_sbs
        end
        
        if exist('coordinate_addi_top_s_b_s','var')%如果该值存在
            clear coordinate_addi_top_s_b_s
        end
        
        if exist('coordinate_sbs_ahead','var')%如果该值存在
            clear coordinate_sbs_ahead
        end
        
        key = key_seq{p};
        eval(['I2_tt = cell2mat(img.' key ');'])
        %         I2_0 = heiban;
        I2 = I2_tt;%这个纯属为了消除警告
        ave_I2 = sum(I2(:)) / (numel(I2(:,1)) * numel(I2(1,:)));
        I2(I2 > ave_I2) = 255;
        I2(I2 <= ave_I2) = 0;%对I2进行二值化
        
        I2 = ~I2;
        I2_x = l_x(I2);
        I2_y = l_y(I2);
        
        I2_guang = trans2lefttop(I2);
        I2_guang = imcrop(I2_guang,[0, 0, l_x(I2_guang), l_y(I2_guang)]);
        
        if (ceil(eval(['Y_length.' key])) + partgap > numel(heiban(:,1))) || (ceil(eval(['X_length.' key])) + partgap > numel(heiban(1,:)))
            temp_x = eval(['X_length.' key]);
            temp_y = eval(['Y_length.' key]);
            eval(['Y_length.' key '= temp_x;']);
            eval(['X_length.' key '= temp_y;']);
            I2_guang = imrotate(I2_guang,90,'bicubic');
            sita_xz = 90;
        else
            sita_xz = 0;
        end
        
        
        I2_t_max_or = imresize(I2_guang,[ceil(eval(['Y_length.' key])), ceil(eval(['X_length.' key]))],'bicubic');
        %         I2_t_max_or = imresize(I2_guang,[round(eval(['Y_length.' key])),NaN],'bicubic');
        I_base = zeros(l_y(I2_t_max_or) + round(2*partgap) + 10, l_x(I2_t_max_or) + round(2 * partgap) + 10);
        I2_t_max_or = paste_img(I_base,I2_t_max_or,1,1);
        I2_t_max_or = trans2mid(I2_t_max_or,0);
        I2_t_max_or_noexpand = I2_t_max_or;
        I2_t_max_or = expand(I2_t_max_or,partgap);
        
        I2_t_max = imresize(I2_guang,[ceil(eval(['Y_length.' key]) * imrf),ceil(eval(['X_length.' key]) * imrf)],'bicubic');
        I_base = zeros(l_y(I2_t_max)+round(2 * imrf * partgap) + 10,l_x(I2_t_max) + round(2 * imrf * partgap) + 10);
        
        I2_t_max = paste_img(I_base,I2_t_max,1,1);
        I2_t_max = trans2mid(I2_t_max,0);
        I2_t_max = expand(I2_t_max,partgap * imrf);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %以下求解是为了求出后期旋转后的中心点
        size_max = max(eval(['Y_length.' key]),eval(['X_length.' key]));
        
        if size_max <= 150
            imrf_special = 10;
        elseif size_max > 150 && size_max <=300
            imrf_special = 5;
        elseif size_max > 300 && size_max <=1000
            imrf_special = 2;
        else
            imrf_special = 1;
        end
        
        I2_t_max_t_10 = imresize(I2_guang,[round(eval(['Y_length.' key]) * imrf * imrf_special),round(eval(['X_length.' key]) * imrf * imrf_special)],'bicubic');
        %         I2_t_max_or_t_10 = imresize(I2_guang,[round(eval(['Y_length.' key])*10),NaN],'bicubic');
        l_max = max(numel(I2_t_max_t_10(:,1)),numel(I2_t_max_t_10(1,:)));
        I_base = zeros(round(l_max * 1.45 + imrf * imrf_special * partgap),round(l_max * 1.45 + imrf * imrf_special * partgap));
        I2_t_max_t_10 = paste_img(I_base,I2_t_max_t_10,1,1);
        I2_t_max_t_10 = trans2mid(I2_t_max_t_10,1);
        I2_t_max_t_10_noexpand = I2_t_max_t_10;
        
        %         I2_t_max_t_10 = trans2mid(I2_t_max_t_10,1);
        %         [nl,nr,nt,nb] = margin(I2_t_max_or_t_10);
        middle_coordinate_x = numel(I2_t_max_t_10_noexpand(1,:)) / (2 * imrf * imrf_special);
        middle_coordinate_y = numel(I2_t_max_t_10_noexpand(:,1)) / (2 * imrf * imrf_special);
        I2_t_max_t_10 = expand(I2_t_max_t_10, partgap * imrf * imrf_special);
        I2_t_max_t_10 = double(I2_t_max_t_10);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         I2_t_max_10 = imresize(I2_guang,[round(eval(['Y_length.' key]) * imrf * imrf_special),round(eval(['X_length.' key]) * imrf * imrf_special)],'bicubic');
        %         I2_t_max_10 = imresize(I2_guang,[round(eval(['Y_length.' key])*imrf*10),NaN],'bicubic');
        %         I2_t_max_10 = trans2lefttop(I2_t_max_t_10);
        %         I2_t_max_10 = imcrop(I2_t_max_10,[0,0,l_x(I2_t_max_10),l_y(I2_t_max_10)]);
        %         I_base = zeros(round(size(I2_t_max_10) * 1.2));
        %         I2_t_max_10 = paste_img(I_base,I2_t_max_10,1,1);
        %         I2_t_max_10 =trans2mid(I2_t_max_10,0);
        %         I2_t_max_10 = expand(I2_t_max_10,partgap * imrf * imrf_special);
        
        I2_t_max = trans2lefttop(I2_t_max);
        I_base = zeros(l_y(I2_t_max) * 2 + round(imrf * 2 * imrf_special),l_x(I2_t_max) * 2 + round(imrf * 2 * imrf_special));
        I2_t_max = paste_img(I_base,I2_t_max,1,1);
        I2_t_max = trans2mid(I2_t_max,0);
        
        I2_0 =  paste_img(heiban,I2_t_max_or,1,1);%这里的I2_0已经发生了膨胀
        I2 = I2_0;
        nC = round(nPop0 * nPc / 2) * 2; % 子代规模的大小
        %         nMu = 1; % 变异的概率
        [Parent0,template,I2_t,~] = Pretreat(I0,I2,nPop0);%图像预处理,这里的I2_t只是单纯的二值化和背景处理，Parent也是空的，没有赋值
        I2_t = trans2mid(I2_t,0);
        %         I2_t_ymin = I2_t;
        %         I2_t_ymin = trans2mid(I2_t_ymin,0);
        %         I2_t = I2_t_ymin;%初始赋值
        I2_t0 = I2_t;%这是为了保留原始图像的尺寸
        
        W0 = min(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的宽度
        L0 = max(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的长度
        
        %         I2_t=expand(I2_t,partgap);
        lx1 = l_x(I2_t);
        ly1 = l_y(I2_t);
        
        %[dxy_struct] = find_best_dxy(heiban,I0,I2_t0,I2_t,I2_t_max_or,I2_t_max_or_10,Parent0,...
        %template,Nmax,key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,nMu,maxIt,W0,L0,kb,partgap,N_part);
        S.alg1=[];S.alg2=[];
        S.dx1=[];S.dy1=[];S.dx2=[];S.dy2=[];S.dx3=[];S.dy3=[];
        S.sita10=[];S.sita1=[];S.beta=[];
        S.I1=[];S.I2=[];S.I3=[];
        S.SP=[];S.XP=[];S.ZP=[];S.YP=[];S.N=[];S.ii = [];S.coordinate = [];
        dxy_struct = repmat(S, 5, 1);
        
        Nmax = Nmax; %#ok<ASGSL>
        kb = kb; %#ok<ASGSL>
        
        isrectangle = IsRectangle(I2_t_max_or);
        
        if  eval(['Nmax.' key ]) >= N_linjie || (p == 1 && eval(['Nmax.' key ]) >= 2)%如果是第一个零件，则只要数目>=2就用遗传算法排图
            delete(gcf),delete(gcf),delete(gcf);
            
            if isrectangle == 1
                find_best_dxy5_used = 1;
                dxy_struct5 = find_best_dxy5(lx_heiban,ly_heiban,I2_t_max,imrf,sita_xz);
                coordinate_single_pre = dxy_struct5.coordinate;
                N_pre = length(coordinate_single_pre);
                
            else
                
                parfor mm = 1:5
                    if mm == 1
                        dxy_struct(mm) = find_best_dxy1(heiban,I2_t0,I2_t,I2_t_max_or,Parent0,template,...
                            key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,partgap,N_times,imrf_special,sita_xz);
                    elseif mm == 2
                        dxy_struct(mm) = find_best_dxy2(heiban,I2_t0,I2_t,I2_t_max_or,Parent0,template,...
                            key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,partgap,N_times,imrf_special,sita_xz);
                    elseif mm == 3
                        dxy_struct(mm) = find_best_dxy3(heiban,I2_t0,I2_t,I2_t_max_or,Parent0,template,...
                            key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,partgap,N_times,imrf_special,sita_xz);
                    elseif mm == 4
                        dxy_struct(mm) = find_best_dxy4(heiban,I2_t0,I2_t,I2_t_max_or,Parent0,template,...
                            key,maxIt_va,nPop,nVar,n,lx1,ly1,nC,maxIt,W0,L0,kb,partgap,N_times,imrf_special,sita_xz);
                    elseif mm == 5
                        if (max(lx_heiban,ly_heiban) / max(l_x(I2_t_max_or),l_y(I2_t_max_or))) < 100
                            dxy_struct(mm) = find_best_dxy5(lx_heiban,ly_heiban,I2_t_max,imrf,sita_xz);
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isempty(dxy_struct(5).N) == 0 %说明其进行了find_best_dxy5的计算
                    coordinate_5 = dxy_struct(5).coordinate;
                end
                
                dxy_struct_1_4 = [dxy_struct(1);dxy_struct(2);dxy_struct(3);dxy_struct(4)];
                [~, so] = sort([dxy_struct_1_4.N], 'descend'); %so为降序后的序号排序
                dxy_best = dxy_struct_1_4(so(1));
                
                if isempty(dxy_struct(5).N) == 1 || (isempty(dxy_struct(5).N) == 0 && dxy_struct(5).N < (dxy_best.N ))                              %说明其进行了find_best_dxy5的计算
                    %如果并未进行find_best_dxy5的计算，或者虽然计算了，但是其结果没有前面的四个结果好，则进行下面的计算
                    
                    if dxy_struct(1).N == dxy_struct(so(1)).N
                        dxy_best = dxy_struct(1);
                    end
                    
                    if dxy_struct(2).N == dxy_struct(so(1)).N
                        dxy_best = dxy_struct(2);
                    end
                    
                    if dxy_struct(3).N == dxy_struct(so(1)).N
                        dxy_best = dxy_struct(3);
                    end
                    
                    if (dxy_struct(4).N == dxy_struct(so(1)).N) || (dxy_struct(4).N >= dxy_struct(so(1)).N - 2)
                        dxy_best = dxy_struct(4);
                    end
                    
                    alg1_t=dxy_best.alg1;
                    alg2_t=dxy_best.alg2;
                    dx1 = dxy_best.dx1;
                    dx2 = dxy_best.dx2;
                    dx3 = dxy_best.dx3;
                    dy1 = dxy_best.dy1;
                    dy2 = dxy_best.dy2;
                    dy3 = dxy_best.dy3;
                    sita10 = dxy_best.sita10;
                    sita1 = dxy_best.sita1;
                    beta = dxy_best.beta;
                    I1 = dxy_best.I1;
                    I2 = dxy_best.I2;
                    I3 = dxy_best.I3;
                    SP = dxy_best.SP;
                    XP = dxy_best.ZP;
                    ZP = dxy_best.ZP;
                    YP = dxy_best.ZP;
                    ii = dxy_best.ii;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %         [dx1,dx2,dx3,dy1,dy2,dy3,sita10,sita1,beta,~,I1,I2,I3,SP,XP,ZP,YP] = ...%beta的值只有两种可能，要么是0，要么是90
                    %             find_position(I2_t,Parent0,template,key,alg1,alg2,maxIt_va,nPop,nVar,n,lx1,ly1,nC,nMu,maxIt,W0,L0,kb,N_times);
                    %上面的结果中，dy2的符号始终非正，dx3始终非负，而dx2和dy3的符号根据ZP和SP的值而定
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %以上完成了三次迭代，下面进行高精度寻优
                    %寻找第1次排图的高精度位置
                    I2_t_max = trans2lefttop(I2_t_max);
                    I2_t_max = imcrop(I2_t_max,[0,0,l_x(I2_t_max),l_y(I2_t_max)]);
                    I_base = zeros(ceil(l_y(I2_t_max)*(l_y(I1)/l_y(I2_t)))+round(8*4*imrf)+2,ceil(l_x(I2_t_max)*(l_x(I1)/l_x(I2_t))+round(8*4*imrf))+2);
                    I2_t_max = paste_img(I_base,I2_t_max,1,1);
                    I2_t_max = trans2mid(I2_t_max,0);
                    %         tx=1;%后期需要修改
                    [dx111,dy111,I1_hr] = high_accuracy1(I2_t_max,alg1_t,key,dx1,dy1,sita10,imrf,0,4,1);
                    
                    %寻找第2次排图的高精度位置
                    I1_hr_t = trans2lefttop(I1_hr);
                    I1_hr_t = imcrop(I1_hr_t,[0,0,l_x(I1_hr_t),l_y(I1_hr_t)]);
                    I1_hr_10 = imresize(I1_hr_t,imrf_special);%这里有问题
                    I1_hr_10(I1_hr_10 < 0.5) = 0;
                    I1_hr_10(I1_hr_10 ~= 0) = 1;
                    I_base = zeros(ceil(l_y(I1_hr_t) * (l_y(I2) / l_y(I1))) + round(10 * partgap * imrf) + 10, ceil(l_x(I1_hr_t) * (l_x(I2) / l_x(I1)) + round(10 * partgap * imrf)) + 10);
                    I1_hr_t = paste_img(I_base,I1_hr_t,1,1);
                    I1_hr_t = trans2mid(I1_hr_t,0);
                    [dx222,dy222,I2_hr] = high_accuracy2(I1_hr_t,alg2_t,key,dx2,dy2,imrf,5,1);%这里的dx2要换
                    
                    if abs(dx222) >= abs(dy222)
                        beta1 = 90;
                    else
                        beta1 = 0;
                    end
                    
                    if beta == beta1
                        ...
                    else %beta ~= beta1
                    if beta == 0 && beta1 == 90
                        beta = 90;
                        dx3_t = -dy3;
                        dy3_t = dx3;
                        dx3 = dx3_t;
                        dy3 = dy3_t;
                    else %beta == 90 && beta1 == 0
                        beta = 0;
                        dx3_t = dy3;
                        dy3_t = -dx3;
                        dx3 = dx3_t;
                        dy3 = dy3_t;
                    end
                    end
                    
                    sita1 = sita1 + beta;
                    if beta == 90
                        dx111_ro =  dy111;
                        dy111_ro = -dx111;
                    else
                        dx111_ro = dx111;
                        dy111_ro = dy111;
                    end
                    
                    if beta == 90
                        dx222_ro =  dy222;
                        dy222_ro = -dx222;
                    else
                        dx222_ro = dx222;
                        dy222_ro = dy222;
                    end
                    
                    %寻找第3次排图的高精度位置
                    sita0 = atan(- dx222_ro / dy222_ro) * 180 / pi;%最终结果180度数制
                    if dy222_ro == 0%目前看没必要了，因为dy222_ro不可能为0
                        kk_end = 10;
                    else
                        kk_end = 1 + floor(abs(40 / dy222_ro));
                    end
                    clear I2_hr_3
                    %         if dx222_ro ~= 0 && dy222_ro ~= 0
                    gaodu = zeros(1,2);
                    
                    I3_t = imrotate(I1_hr,beta,'bicubic');
                    I3_t = trans2lefttop(I3_t);
                    I3_t = imcrop(I3_t,[0,0,l_x(I3_t),l_y(I3_t)]);
                    I_base = zeros(l_y(I3_t) + round((2 * kk_end + 1) * abs(dy222_ro) * imrf),l_x(I3_t) + round((2 * kk_end + 1) * abs(dx222_ro) * imrf));
                    I3_t = paste_img(I_base,I3_t,1,1);
                    I3_t = trans2mid(I3_t,0);
                    I3_tt = I3_t;
                    for kk = 1:kk_end
                        I3_tt = imtranslate(I3_tt,[dx222_ro*imrf,dy222_ro*imrf]);
                        I3_t = I3_t | I3_tt;
                    end
                    I3_t = trans2lefttop(I3_t);
                    I3_t = imcrop(I3_t,[0,0,l_x(I3_t),l_y(I3_t)]);
                    
                    Ymax = round(1.1*(l_y(I3_t) + 2 * abs(dy3) * imrf + round(2 * ii * abs(dy222_ro)) * imrf));
                    Xmax = round(1.1*(l_x(I3_t) + 2 * abs(dx3) * imrf + round(2 * ii * abs(dx222_ro)) * imrf));
                    N_I3_t_pix = sum(I3_t(:));
                    
                    I_base = zeros(Ymax,Xmax);
                    I3_t = paste_img(I_base,I3_t,1,1);
                    I3_t = trans2mid(I3_t,0);
                    I3_t_r = imtranslate(I3_t,[round(dx3 * imrf - ii * dx222_ro * imrf), round(dy3 * imrf - ii * dy222_ro * imrf)]);
                    I3_t_r1 = imtranslate(I3_t,[round(dx3 * imrf + ii * dx222_ro * imrf), round(dy3 * imrf + ii * dy222_ro * imrf)]);
                    
                    while sum(I3_t_r(:)) ~= N_I3_t_pix || sum(I3_t_r1(:)) ~= N_I3_t_pix
                        Ymax = round(1.2 * Ymax);
                        Xmax = round(1.2 * Xmax);
                        I_base = zeros(Ymax,Xmax);
                        I3_t = paste_img(I_base,I3_t,1,1);
                        I3_t = trans2mid(I3_t,0);
                        I3_t_r = imtranslate(I3_t,[round(dx3 * imrf - ii * dx222_ro * imrf), round(dy3 * imrf - ii * dy222_ro * imrf)]);
                        I3_t_r1 = imtranslate(I3_t,[round(dx3 * imrf + ii * dx222_ro * imrf), round(dy3 * imrf + ii * dy222_ro * imrf)]);
                    end
                    
                    %         l_y_last = inf;
                    %         I3_t_r = imtranslate(I3_t,[dx3 * imrf - ii * dx222_ro * imrf, dy3 * imrf - ii * dy222_ro * imrf]);
                    I3_t_combine = I3_t | I3_t_r;
                    I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
                    x_del = min(round(9.8 * l_x(I3_t_combine_ro) / 20),round(l_x(I3_t_combine_ro)/2)-4);
                    m = 0;
                    for i = -ii : ii
                        m = m + 1;
                        I3_t_r = imtranslate(I3_t,[dx3 * imrf + i * dx222_ro * imrf, dy3 * imrf + i * dy222_ro * imrf]);
                        I3_t_combine = I3_t | I3_t_r;
                        I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
                        I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
                        I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
                        
                        while (sum(I3_t_combine_ro_t(:))) <= (sum(I3_t_combine_ro(:)) * 0.501)
                            x_del = min(round(x_del * 0.8), x_del-2);
                            I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
                            I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
                        end
                        
                        I3_t_combine_ro_t = trans2righttop(I3_t_combine_ro_t);
                        I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[ x_del,0]);
                        I3_t_combine_ro_t = trans2mid(I3_t_combine_ro_t,0);
                        gaodu(m,1) = i;
                        gaodu(m,2) = l_y(I3_t_combine_ro_t);
                    end
                    
                    [~, so] = sort(gaodu(:,2), 'ascend');
                    i_best = gaodu(so(1),1);
                    
                    while i_best == abs(ii)
                        ii = ii + 1;
                        Ymax = round(1.1*(l_y(I3_t) + 2 * abs(dy3) * imrf + round(2 * ii * abs(dy222_ro)) * imrf));
                        Xmax = round(1.1*(l_x(I3_t) + 2 * abs(dx3) * imrf + round(2 * ii * abs(dx222_ro)) * imrf));
                        I_base = zeros(Ymax,Xmax);
                        I3_t = paste_img(I_base,I3_t,1,1);
                        I3_t = trans2mid(I3_t,0);
                        I3_t_r = imtranslate(I3_t,[round(dx3 * imrf - ii * dx222_ro * imrf), round(dy3 * imrf - ii * dy222_ro * imrf)]);
                        I3_t_r1 = imtranslate(I3_t,[round(dx3 * imrf + ii * dx222_ro * imrf), round(dy3 * imrf + ii * dy222_ro * imrf)]);
                        
                        while sum(I3_t_r(:)) ~= N_I3_t_pix || sum(I3_t_r1(:)) ~= N_I3_t_pix
                            Ymax = round(1.2 * Ymax);
                            Xmax = round(1.2 * Xmax);
                            I_base = zeros(Ymax,Xmax);
                            I3_t = paste_img(I_base,I3_t,1,1);
                            I3_t = trans2mid(I3_t,0);
                            I3_t_r = imtranslate(I3_t,[round(dx3 * imrf - ii * dx222_ro * imrf), round(dy3 * imrf - ii * dy222_ro * imrf)]);
                            I3_t_r1 = imtranslate(I3_t,[round(dx3 * imrf + ii * dx222_ro * imrf), round(dy3 * imrf + ii * dy222_ro * imrf)]);
                        end
                        
                        m = 0;
                        
                        for i=-ii:ii
                            m = m + 1;
                            I3_t_r = imtranslate(I3_t,[dx3 * imrf + i * dx222_ro * imrf, dy3 * imrf + i * dy222_ro * imrf]);
                            %imshow(I3_t_r);
                            I3_t_combine = I3_t | I3_t_r;
                            I3_t_combine_ro = imrotate(I3_t_combine,sita0,'bicubic');
                            %x_del = round(9 * l_x(I3_t_combine_ro) / 20);
                            I3_t_combine_ro_t = trans2lefttop(I3_t_combine_ro);
                            I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
                            
                            while (sum(I3_t_combine_ro_t(:))) <= (sum(I3_t_combine_ro(:)) / 2)
                                x_del = min(round(x_del * 0.8), x_del-2);
                                I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[-x_del,0]);
                            end
                            
                            I3_t_combine_ro_t = trans2righttop(I3_t_combine_ro_t);
                            I3_t_combine_ro_t = imtranslate(I3_t_combine_ro_t,[ x_del,0]);
                            I3_t_combine_ro_t = trans2mid(I3_t_combine_ro_t,0);
                            %imshow(I3_t_combine_ro_t),pause(0.5)
                            gaodu(m,1) = i;
                            l_y_now = l_y(I3_t_combine_ro_t);
                            gaodu(m,2) = l_y_now;
                        end
                        
                        [~, so] = sort(gaodu(:,2), 'ascend');
                        i_best = gaodu(so(1),1);
                    end
                    
                    dx3 = dx3 + i_best * dx222_ro;
                    dy3 = dy3 + i_best * dy222_ro;
                    
                    I2_hr_3 = I3_t;
                    
                    I2_hr_3_guang = trans2lefttop(I2_hr_3);
                    I2_hr_3_guang = imcrop(I2_hr_3_guang,[0,0,l_x(I2_hr_3_guang),l_y(I2_hr_3_guang)]);
                    I2_hr_3_10 = imresize(I2_hr_3_guang,imrf_special);
                    I2_hr_3_10(I2_hr_3_10 < 0.5) = 0;
                    I2_hr_3_10(I2_hr_3_10 ~= 0) = 1;
                    %         I2_hr_3_10 = trans2lefttop(I2_hr_3_10);
                    %         I2_hr_3_10 = imcrop(I2_hr_3_10,[0,0,l_x(I2_hr_3_10),l_y(I2_hr_3_10)]);
                    
                    I2_hr = trans2lefttop(I2_hr);
                    I2_hr = imcrop(I2_hr,[0,0,l_x(I2_hr),l_y(I2_hr)]);
                    I2_hr_10 = imresize(I2_hr,imrf_special);
                    I2_hr_10(I2_hr_10 < 0.5) = 0;
                    I2_hr_10(I2_hr_10 ~= 0) = 1;
                    %         I2_hr_t = imrotate(I2_hr_3,beta);
                    I2_hr_t = I2_hr_3_guang;
                    sita0 = atan(-dx222_ro/dy222_ro)*180/pi;%最终结果180度数制
                    
                    I_base = zeros(1*l_y(I2_hr_t)+(round(abs(dy3))+60)*imrf,1*l_x(I2_hr_t)+(round(abs(dx3))+60)*imrf);
                    
                    I2_hr_t = paste_img(I_base,I2_hr_t,1,1);
                    I2_hr_t = trans2mid(I2_hr_t,0);
                    [dx333_ro,dy333_ro,I3_hr] = high_accuracy3(I2_hr_t,dx3,dy3,sita0,imrf,6,1);
                    I3_hr = trans2lefttop(I3_hr);
                    I3_hr = imcrop(I3_hr,[0,0,l_x(I3_hr),l_y(I3_hr)]);
                    I3_hr_10 = imresize(I3_hr,imrf_special);
                    I3_hr_10(I3_hr_10 < 0.5) = 0;
                    I3_hr_10(I3_hr_10 ~= 0) = 1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %以上已经找到了高精度的位置，下面进行排图
                    if beta == 90
                        ZP = 1 - ZP;
                        YP = 1 - YP;
                    end
                    dx1=dx111_ro;
                    dy1=dy111_ro;
                    
                    dx2=dx222_ro;
                    dy2=dy222_ro;
                    
                    if dy2 > 0
                        dx2=-dx2;
                        dy2=-dy2;
                    end
                    
                    dx3=dx333_ro;
                    dy3=dy333_ro;
                    %[dx1,dy1,dx2,dy2,dx3,dy3,sita1] = find_suit_dxy(ZP,dx1,dy1,dx2,dy2,dx3,dy3,sita1);%最终不能确定dx3的符号，但是可以确定dy3的符号为非负
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %接下来将图排成水平，即将dy3变为0
                    if exist('i_best','var')
                        gama = atan((dy333_ro - i_best * dy222_ro) / (dx333_ro - i_best * dx222_ro)) * 180 / pi;
                        %gama为为了将图像进行水平排布为旋转的角度,有问题，因为有可能dx3和dy3不是相邻的两列的差值
                    else
                        gama = atan(dy333_ro / dx333_ro) * 180 / pi;
                    end%经过验证，结果正确
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [dx1,dy1,dx2,dy2,dx3,dy3,sita1] = find_suit_dxy(ZP,dx1,dy1,dx2,dy2,dx3,dy3,sita1);%最终不能确定dx3的符号，但是可以确定dy3的符号为非负
                    %上面的表达式主要是调整dx3和dy3的值，并未对dx2和dy2做调整
                    I1_hr_or = I1_hr;
                    %下面将零件摆平，并寻找合适的间距值
                    [dx1,dy1,dx2,dy2,dx3,dy3,sita1,ZP,YP,~,~,~,~,~] = baiping(I1_hr,I2_hr,I3_hr,I2_hr_3,I2_hr_3_10,I2_t_max,...
                        I1_hr_10,I2_hr_10,I3_hr_10,I2_t_max_t_10,sita1,beta,gama,ZP,YP,SP,XP,dx1,dy1,dx2,dy2,dx3,dy3,imrf,imrf_special);
                    dx3 = abs(dx3);
                    %         imshow(I3_hr);
                    %以上所所得到的dx1,dx2,dx3,dy1,dy2,dy3均为摆平后的值
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %beta
                    sita2 = atan(-dx2/dy2);
                    sita0 = sita2*180/pi;%最终结果180度数制
                    
                    I2_shu = imrotate(I1,sita1,'bicubic');
                    lmax = round(2*(max(l_y(I1)*2+abs(dy2)+abs(dy3),l_x(I1)*2+abs(dx2)+abs(dx3))));
                    I_base = zeros(lmax,lmax);
                    I2_shu = trans2lefttop(I2_shu);
                    I2_shu = imcrop(I2_shu,[0,0,l_x(I2_shu),l_y(I2_shu)]);
                    I2_shu = paste_img(I_base,I2_shu,1,1);
                    I2_shu = trans2mid(I2_shu,0);
                    %     I2_shu = imrotate(I2_shu,beta);
                    %     I2_shu = trans2rightbottom(I2_shu);
                    I2_shu_r = I2_shu;
                    I2_shu_r = imtranslate(I2_shu_r,[dx2,dy2]);
                    I2_shu = I2_shu | I2_shu_r;%正确
                    
                    I3_shu = I2_shu;
                    I3_shu = trans2mid(I3_shu,0);
                    I3_shu_r = I3_shu;
                    I3_shu_r = imtranslate(I3_shu_r,[dx3,dy3]);
                    I3_shu = I3_shu | I3_shu_r;
                    
                    I2_shu = trans2mid(I2_shu,0);
                    I3_shu = trans2mid(I3_shu,0);
                    I2_shu = imrotate(I2_shu,sita0,'bicubic');
                    I3_shu = imrotate(I3_shu,sita0,'bicubic');
                    %         imshow(I3_shu)
                    
                    d3_shu = l_x(I3_shu)-l_x(I2_shu); %这里有问题
                    d3_shu = floor(d3_shu / cos(sita2));%/imrf;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %         bujin =round(0.8 * l_x(logical(imrotate(expand(I2_t0, partgap), sita1, 'bicubic'))));
                    N_shu = ceil(numel(heiban(:,1)) / abs(dy2));  %结果正确 其表示纵向最多能排多少个图案
                    delta_xx = abs(numel(heiban(:,1)) * tan(sita2));%斜铺后横向多出来得到部分像素（实际上不存在）
                    x_total = delta_xx + numel(heiban(1,:));
                    N_dx3 = round(ceil(numel(heiban(1,:)) / d3_shu) * 1.4);%仅仅在板子上，横向走能排多少个
                    N_dx3_max = round(ceil(x_total / d3_shu) * 1.4);%+3000;%加上扩展的区域，横向走能排多少个
                    
                    N_flag = 0;%其为排图的数目达到要求的标记
                    single = 1;%单排的标记
                    
                    d3_nitial = 0;%dx3和dy3是否被初始调整过的标记
                    coordinate_single = [];
                    coordinate_multiple = [];
                    
                    %         imshow(I2_t_max_or_t_10)
                    %I2_t_max_t_10_ro = imrotate(I2_t_max_t_10,sita1,'bicubic','crop');
                    I2_t_max_t_10_ro = imrotate(I2_t_max_t_10_noexpand,sita1,'bicubic','crop');
                    
                    I2_t_max_t_10_ro(I2_t_max_t_10_ro < 0.5) = 0;
                    I2_t_max_t_10_ro(I2_t_max_t_10_ro ~= 0) = 1;
                    
                    [nl,~,~,nb] = margin(I2_t_max_t_10_ro);
                    middle_coordinate_x_ro = (nl + l_x(I2_t_max_t_10_ro) / 2) / (imrf * imrf_special);
                    middle_coordinate_y_ro = (nb + l_y(I2_t_max_t_10_ro) / 2) / (imrf * imrf_special);
                    
                    middle_dx = middle_coordinate_x - middle_coordinate_x_ro;
                    middle_dy = middle_coordinate_y - middle_coordinate_y_ro;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %上面是调整板材尺寸
                    %下面是真正的着手排布零件
                    if wp_new == 1
                        %I2_t0_ex = I2_t_max_or_noexpand;
                        I2_t0_ex = I2_t_max_or;%已经膨胀过
                        %             I_base = zeros(size(I_mix_tt));
                        I2_t0_ex = trans2lefttop(I2_t0_ex);
                        I2_t0_ex = imcrop(I2_t0_ex,[0,0,l_x(I2_t0_ex),l_y(I2_t0_ex)]);
                        I2_t_temp0 = I2_t0_ex;
                        
                        I_s_m_guang = imrotate(I2_t_temp0,sita1,'bicubic');
                        I_s_m_guang (I_s_m_guang < 0.5) = 0;
                        I_s_m_guang = logical(I_s_m_guang);
                        I_s_m_guang = trans2lefttop(I_s_m_guang);
                        I_s_m_guang = imcrop(I_s_m_guang,[0,0,l_x(I_s_m_guang),l_y(I_s_m_guang)]);
                        lx_single = l_x(I_s_m_guang);
                        ly_single = l_y(I_s_m_guang);
                    end
                    
                    I_s_m = paste_img(heiban,I_s_m_guang,1,1);
                    I_s_m = trans2leftbottom(I_s_m);%单个移动图元至左下角，此时还不知道dy1是不是大于0，如果dy1大于0，说明其下偏，那么I_s_m_r显示不完全
                    
                    if wp_new  == 1
                        wp_new = 0;
                        I_s_m_r_guang = imrotate(I2_t_temp0,sita1+180,'bicubic');%单个移动图旋转180
                        I_s_m_r_guang (I_s_m_r_guang < 0.5) = 0;
                        I_s_m_r_guang = logical(I_s_m_r_guang);
                        I_s_m_r_guang = trans2lefttop(I_s_m_r_guang);
                        %         I_base = zeros(numel(I0(:,1)),l_x(I_s_m_r));
                        I_s_m_r_guang = imcrop(I_s_m_r_guang,[0,0,l_x(I_s_m_r_guang),l_y(I_s_m_r_guang)]);
                    end
                    
                    I_s_m_r = paste_img(heiban,I_s_m_r_guang,1,1);
                    I_s_m_r = trans2leftbottom(I_s_m_r);%单个移动图元至左下角，此时还不知道dy1是不是大于0，如果dy1大于0，说明其下偏，那么I_s_m_r显示不完全
                    
                    if dy1 > 0
                        I_s_m = imtranslate(I_s_m,[0,round(-dy1)]);
                        I_s_m_r = imtranslate(I_s_m_r,[0,round(-dy1)]);
                    end
                    
                    if dx1 < 0
                        I_s_m = imtranslate(I_s_m,[-dx1,0]);
                        I_s_m_r = imtranslate(I_s_m_r,[-dx1,0]);
                    end
                    
                    I_s_m_r = imtranslate(I_s_m_r,[round(dx1),round(dy1)]);
                    I_s_m0 = I_s_m;
                    I_accm = I_s_m0;%初始化
                    I_s_m_r0 = I_s_m_r;
                    I_accm_r = I_s_m_r0;%初始化
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %在进行排图之前，先进行单个零件的预排
                    if exist('N_pre','var') == 0 %如果该值不存在
                        sn_template.N = [];
                        sn_template.C1 = [];
                        sn_template.C2 = [];
                        sn_template.XPY = [];
                        sn_template.YPY = [];
                        SN = repmat(sn_template, abs(round(dx3) * round(dy2)), 1);
                        
                        q=0;
                        for x_py = 0:-2:-abs(round(dx3))
                            for y_py = 0:2:abs(round(dy2))
                                q = q + 1;
                                [coordinate1_pre_t,coordinate2_pre_t,N_single] = paitu_pre(ZP,dx1,dx2,dx3,dy1,dy2,dy3,sita1,N_dx3,N_dx3_max,N_shu,heiban,lx_single,ly_single,x_py,y_py,sita_xz);
                                SN(q).N = -N_single;
                                SN(q).C1 = coordinate1_pre_t;
                                SN(q).C2 = coordinate2_pre_t;
                                SN(q).XPY = x_py;
                                SN(q).YPY = -y_py;
                            end
                        end
                        
                        N_array = [SN.N; SN.XPY; SN.YPY]';
                        [N_array_sorted,index] = sortrows(N_array, [1 2 3]);
                        index_best = index(1);
                        
                        coordinate1_pre = SN(index_best).C1;
                        coordinate2_pre = SN(index_best).C2;
                        clear SN
                        
                        %以上得到的坐标，要先进行上层的排图，下面的进行忽略。
                        [coordinate_addi] = addi_top(dx3,dy2,coordinate1_pre,coordinate2_pre,middle_coordinate_x,middle_coordinate_y,...
                            heiban,I_s_m0,I_s_m_r0,I2_t_temp0,I2_t_max_t_10,I1,lx_single,ly_single,imrf,imrf_special,I2_t_max_t_10,I1_hr_10,I2_t_max,I1_hr_or);
                        
                        if isempty(coordinate_addi)
                            clear coordinate_addi
                        end
                        
                        [x_ave1, y_ave1] = cor_average(I_s_m0);
                        [x_ave2, y_ave2] = cor_average(I_s_m_r0);
                        coordinate1_pre(:,4) = round((coordinate1_pre(:,4) + x_ave1) / (lx_single / 2));
                        coordinate1_pre(:,5) = -coordinate1_pre(:,5) + y_ave1;
                        coordinate2_pre(:,4) = round((coordinate2_pre(:,4) + x_ave2) / (lx_single / 2));
                        coordinate2_pre(:,5) = -coordinate2_pre(:,5) + y_ave2;
                        
                        if isempty(coordinate1_pre) == 0
                            coordinate1_pre(:,2) = -coordinate1_pre(:,2);
                            if dy1 > 0
                                coordinate1_pre(:,2) = coordinate1_pre(:,2) + dy1;
                            end
                            if dx1 < 0
                                coordinate1_pre(:,1) = coordinate1_pre(:,1) - dx1;
                            end
                            coordinate1_pre(:,1) = coordinate1_pre(:,1) + l_x(I_s_m) / 2 + middle_dx;
                            coordinate1_pre(:,2) = coordinate1_pre(:,2) + l_y(I_s_m) / 2 + middle_dy;
                        end
                        
                        if isempty(coordinate2_pre) == 0
                            coordinate2_pre(:,2) = -coordinate2_pre(:,2);
                            if dy1 > 0
                                coordinate2_pre(:,2) = coordinate2_pre(:,2) + dy1;
                            end
                            if dx1 < 0
                                coordinate2_pre(:,1) = coordinate2_pre(:,1) - dx1;
                            end
                            coordinate2_pre(:,1) = coordinate2_pre(:,1) + dx1 + l_x(I_s_m) / 2 - middle_dx;
                            coordinate2_pre(:,2) = coordinate2_pre(:,2) - dy1 + l_y(I_s_m) / 2 - middle_dy;
                        end
                        coordinate_single_pre = [coordinate1_pre; coordinate2_pre];
                        
                        if exist('coordinate_addi','var') %如果该值存在
                            coordinate_single_pre = [coordinate_single_pre;coordinate_addi]; %#ok<AGROW>
                            clear coordinate_addi
                        end
                        
                        if isempty(coordinate_single_pre) == 0%如果其存在
                            coordinate_single_pre(:,1) = coordinate_single_pre(:,1) + leftmargin - partgap / 2;
                            coordinate_single_pre(:,2) = coordinate_single_pre(:,2) + downmargin - partgap / 2;
                        end
                        coordinate_single_pre = sortrows(coordinate_single_pre,[4 5]);
                        
                        if exist('coordinate_5','var')
                            if length(coordinate_5) >= length(coordinate_single_pre)
                                N_pre = length(coordinate_5);
                            else
                                N_pre = length(coordinate_single_pre);
                            end
                        else
                            N_pre = length(coordinate_single_pre);
                        end
                    end
                    
                else
                    
                    find_best_dxy5_used = 1;
                    coordinate_single_pre = dxy_struct(5).coordinate;
                    N_pre = length(coordinate_single_pre);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if eval(['Nmax.' key ' >= N_pre'])%说明单排没有排完，即板子上排满了，还有剩余的，或者刚好排完，那么需要保存单排的结果
                %如果上述条件不成立，则说明一个板子能排完，则直接进行后面的混排
                danpai_exist = 1;%说明单排存在
                eval(['Nmax.' key ' = mod(Nmax.' key ',N_pre);'])
                
                if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                    eval(['data_list_o.result.single.' key '.partquantity = ' num2str(0) ';'])%正确
                end
                
                if ~isfield(data_list_o.result,'single')%如果data_list_o.result这个键不存在
                    eval(['data_list_o.result.single.' key '.partquantity = ' num2str(0) ';'])%正确
                end
                
                if ~isfield(data_list_o.result.single,key)
                    eval(['data_list_o.result.single.' key '.partquantity = ' num2str(0) ';'])%正确
                end
                eval(['data_list_o.result.single.' key '.partquantity = data_list_o.result.single.' key '.partquantity + N_pre;'])
                
                %                     if eval(['Nmax.' key]) == 0 || flag_N_h == 1 %零件排完了或者板子达到了最大
                single = 0;
                single_file_name = strcat('single_', key);
                %eval(['data_list_o.result.single.' key '.pngfilepath = strcat(outputpath,"\",single_file_name,".png");'])
                eval(['data_list_o.result.single.' key '.txtfilepath = strcat(outputpath,"\",single_file_name,".txt");'])
                %eval(['data_list_o.result.single.' key '.partquantity = ' num2str(N) ';'])%正确
                eval(['sheetquantity = floor((data_list_o.part.' key '.quantity) / N_pre);'])%正确
                eval(['data_list_o.result.single.' key '.sheetquantity = ' num2str(sheetquantity) ';'])%正确
                eval(['data_list_o.result.single.' key '.specifysheet = sheetsize;'])%正确
                ita = sum(I2_t0(:)) * eval(['data_list_o.result.single.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                eval(['data_list_o.result.single.' key '.materialutilizationratio = ' num2str(ita) ';'])%正确
                
                if exist('coordinate_5','var')
                    if length(coordinate_5) >= length(coordinate_single_pre)
                        coordinate_single_pre = coordinate_5;
                        coordinate_single = coordinate_single_pre(:,1:3);
                    else
                        coordinate_single = coordinate_single_pre(:,1:3);
                    end
                else
                    coordinate_single = coordinate_single_pre(:,1:3);
                end
                save(strcat(outputpath,'\',single_file_name,'.txt'), 'coordinate_single', '-ascii')
                %             single_file_name = strcat(single_file_name,'A');
                %             save(strcat(outputpath,'\',single_file_name,'.txt'), 'coordinate_5', '-ascii')
                coordinate_single = [];
            end
            
            if strcmp(enablemix,'yes')%表示混排
                
                while eval(['Nmax.' key ' > 0'])
                    
                    if sum(I_mix(:)) > 0 && eval(['Nmax.' key ' < N_linjie']) %如果零件不是板子上的第一个，且数目较少
                        %则进行single by single排布
                        N_left = eval(['Nmax.' key]);
                        if sum(I_mix(:)) > 0
                            %[I1_t,dx1_t,dy1_t,sita1_t] = find_I1_t(I2_t,I2_t_max_or,Parent0,4,template,key,maxIt_va,N_times,nPop,nVar,nC,n,maxIt);
                            [I_mix,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz);
                        else
                            [coordinate,I_mix,~] = paitu_mix(Nmax,key,sita1,double(I2_t),I_mix,coordinate_single_pre,0,find_best_dxy5_used,sita_xz);
                        end
                        
                        if isempty(coordinate_sbs)%如果是空的
                            N_pai = 0;%p排的零件的数目
                        else
                            N_pai = length(coordinate_sbs(:,1));
                            N_left = N_left - N_pai;
                        end
                        
                        eval(['Nmax.' key ' = Nmax.' key '- N_pai;'])
                        
                        if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result,'multiple')%如果multiple这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result.multiple', key )%如果这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                        end
                        
                        multiple_file_name = strcat('multiple_',nameahead,key,'R_',num_forname);
                        
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = data_list_o.result.multiple{N_multiple,1}.' key '.partquantity + N_pai;'])
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                        eval('data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;')
                        eval('data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;')%正确
                        ita = sum(I2_t0(:)) * eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                        eval(['data_list_o.result.multiple{N_multiple,1}.materialutilizationratio = ' num2str(ita) ';'])%正确
                        
                        save(strcat(outputpath,'\',multiple_file_name,'.txt'), 'coordinate_sbs', '-ascii')
                        
                        if board_is_full == 1 %说明板子排满了
                            N_multiple = N_multiple + 1;
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                            I_mix = heiban;
                        end
                        
                    else
                        
                        if find_best_dxy5_used == 1
                            sita1 = 0;
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if (max(l_x(I2_t_max_or),l_y(I2_t_max_or))) > (min(numel(heiban(:,1)),numel(heiban(1,:))) / 3)
                            c1 = coordinate_single_pre(:,1);
                            c2 = coordinate_single_pre(:,2);
                            ue1 = unique(c1);
                            ue2 = unique(c2);
                            if length(ue1) > length(ue2)
                                coordinate_single_pre = sortrows(coordinate_single_pre,[2 1]);
                            else
                                coordinate_single_pre = sortrows(coordinate_single_pre,[1 2]);
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [coordinate,I_mix,N] = paitu_mix(Nmax,key,sita1,double(I2_t),I_mix,coordinate_single_pre,0,find_best_dxy5_used,sita_xz);
                        eval(['Nmax.' key ' = Nmax.' key '- N;'])
                        
                        I_mix_show = imresize(I_mix,0.5);
                        figure,imshow(I_mix_show,'Border','tight')
                        set(gcf,'ToolBar','none','MenuBar','none','NumberTitle','off');
                        %set(gcf, 'CloseRequestFcn', @(src,event) myCloseFcn(src,event,h));
                        pause(0.0001)
                        
                        N_white_pix = N_white_pix + sum(I2_t_max_or_noexpand(:)) * N;
                        coordinate_multiple = coordinate;
                        I_mix_multiple = I_mix;
                        
                        if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result,'multiple')%如果multiple这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if isfield(data_list_o.result,'multiple')%如果这个键存在
                            if numel(data_list_o.result.multiple) < N_multiple
                                eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            end
                        end
                        
                        if ~isfield(data_list_o.result.multiple{N_multiple,1}', key )%如果这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                        end
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = data_list_o.result.multiple{N_multiple,1}.' key '.partquantity + N;'])
                        
                        %             if (eval(['Nmax.' key ' == 0'])) || (flag_N_h) %某个零件排完了，或者板子达到最大尺寸了,都满足保存的条件
                        %             if (eval(['Nmax.' key ' == 0'])) || (flag_N_h) %如果零件排完了，需要保存，如果没有排完，说明板子排满了，也需要保存
                        multiple_file_name = strcat('multiple_',nameahead,key,'R_',num_forname);
                        %multiple_file_name = strcat('multiple_',nameahead,key,'R_',num_forname,'_',num2str(ppp));
                        %eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                        %写json文件内容
                        
                        %eval(['data_list_o.result.multiple.' key '.sheetquantity = ' num2str(sheetquantity) ';'])%正确
                        %eval(['data_list_o.result.multiple.' key '.specifysheet = sheetsize;'])%正确
                        
                        if (eval(['Nmax.' key ' == 0'])) %&& (p ~= N_part) %板子没有排满，但是某个零件排完了，需要保存坐标信息，坐标信息而且保存在ZB结构体中
                            if isempty(coordinate_multiple) == 0
                                coordinate_multiple(:,1) = coordinate_multiple(:,1) + leftmargin - partgap / 2;
                                coordinate_multiple(:,2) = coordinate_multiple(:,2) + downmargin - partgap / 2;
                            end
                            
                            eval(['ZB.' multiple_file_name ' = coordinate_multiple;'])
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                            coordinate_multiple = coordinate_multiple(:,1:3);
                            data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;
                            data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;
                            save(strcat(outputpath,'\',multiple_file_name,'.txt'), 'coordinate_multiple', '-ascii')
                            coordinate_multiple = [];
                        end
                        
                        if (eval(['Nmax.' key ' > 0'])) % 说明板子排满了，就要输出图，而且要输出txt,不过这里先不急着输出txt，先放在ZB结构体中
                            if isempty(coordinate_multiple) == 0
                                coordinate_multiple(:,1) = coordinate_multiple(:,1) + leftmargin;
                                coordinate_multiple(:,2) = coordinate_multiple(:,2) + downmargin;
                            end
                            
                            eval(['ZB.' multiple_file_name ' = coordinate_multiple;'])
                            %imwrite(I_mix,strcat(outputpath,'\',multiple_file_name,'.png'),'png');
                            %multiple_file_name=strcat(num2str(ppp),multiple_file_name);
                            imwrite(I_mix,strcat(outputpath,'\',multiple_file_name,'.png'),'png');
                            coordinate_multiple = [];
                            data_list_o.result.multiple{N_multiple,1}.pngfilepath = strcat(outputpath,'\',multiple_file_name,'.png');
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                            data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;
                            data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;
                            ita = N_white_pix / (numel(heiban(:,1)) * numel(heiban(1,:)));
                            data_list_o.result.multiple{N_multiple,1}.materialutilizationratio = ita; %这里计算板子利用率
                            N_white_pix = 0;
                        end
                        
                        %if eval(['Nmax.' key_seq{length(keys)} ' > 0'])  &&  flag_N_h == 1 %虽然板子铺满了，但是零件还没有排完,这个时候需要重新用一个空板子
                        if eval(['Nmax.' key ' > 0'])%  &&  flag_N_h == 1 %虽然板子铺满了，但是零件还没有排完,这个时候需要重新用一个空板子
                            %                 flag_N_h = 0;
                            %                 if exist('I_mix_tt','var') %如果该值存在
                            %                     clear I_mix_tt
                            %                 end
                            n_w = 1;%n_w表示混排的时候，同一个板子上的第几种零件
                            N_multiple = N_multiple + 1;%N_multiple为混排的板子的个数
                            %                 clear I_mix_tt
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                            
                            key_ZB = fieldnames(ZB);
                            for i = 1:length(key_ZB)
                                coor_temp = eval(['ZB.' cell2mat(key_ZB(i))]);
                                if isempty(coor_temp) == 0
                                    coor_temp(:,1) = coor_temp(:,1) + leftmargin - partgap / 2;
                                    coor_temp(:,2) = coor_temp(:,2) + downmargin - partgap / 2;
                                    coor_temp = coor_temp(:,1:3);
                                end
                                save(strcat(outputpath,'\',cell2mat(key_ZB(i)),'.txt'), 'coor_temp', '-ascii')
                            end
                            clear ZB
                            I_mix = heiban;
                        end
                    end
                end
                
            else%单排
                
                N_left = eval(['Nmax.' key ]);
                coordinate_single_left = coordinate_single_pre(1 : N_left, :);
                coordinate_single_left = coordinate_single_left(:, 1 : 3);
                single_left_file_name = strcat('multiple_', key);
                
                if ~isfield(data_list_o.result,'multiple')%如果multiple这个键不存在
                    eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                    nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                end
                
                if ~isfield(data_list_o.result.multiple', key )%如果这个键不存在
                    eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                end
                
                eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = data_list_o.result.multiple{N_multiple,1}.' key '.partquantity + N_left;'])
                eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",single_left_file_name,".txt");'])
                eval('data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;')
                eval('data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;')%正确
                ita = sum(I2_t0(:)) * eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                eval(['data_list_o.result.multiple{N_multiple,1}.materialutilizationratio = ' num2str(ita) ';'])%正确
                
                save(strcat(outputpath,'\',single_left_file_name,'.txt'), 'coordinate_single_left', '-ascii')
                N_multiple = N_multiple + 1;
            end
            
            fprintf('图形 %d 所用时间: %f \n',p, toc(loop_time));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        else %因为数目没有达到临界值，所以进行single by single排图
            
            N_left = eval(['Nmax.' key ]);
            
            if strcmp(enablemix,'yes')%表示混排
                
                while N_left > 0
                    
                    if sum(I_mix(:)) > 0
                        
                        [~,~,nt,~] = margin(I_mix);
                        if ((nt > l_x(I2_t_max_or)) || (nt > l_y(I2_t_max_or))) && N_left > 1
                            [I1_t,dx1_t,dy1_t,sita1_t] = find_I1_t(I2_t,I2_t_max_or,Parent0,4,template,key,maxIt_va,N_times,nPop,nVar,nC,n,maxIt);
                            
                            if l_y(I1_t) <= nt
                                [I_mix,coordinate_addi_top_s_b_s,N_left] = addi_top_single_by_single(I_mix,I1_t,I2_t_max_or,...
                                    heiban,dx1_t,dy1_t,sita1_t,N_left);
                                board_is_full = 0;
                            elseif l_x(I1_t) <= nt %旋转后铺图
                                I1_t_ro = imrotate(I1_t,90,'bicubic');
                                I2_t_max_or_ro = imrotate(I2_t_max_or,90,'bicubic');
                                dx1_t_ro = dy1_t;
                                dy1_t_ro = -dx1_t;
                                [I_mix,coordinate_addi_top_s_b_s,N_left] = addi_top_single_by_single(I_mix,I1_t_ro,I2_t_max_or_ro,...
                                    heiban,dx1_t_ro,dy1_t_ro,sita1_t + 90,N_left);
                            else %如果上面的间隙无法容纳下I1_t,则直接进行paitu_single_by_single
                                [I_mix,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz);
                            end
                        else
                            [I_mix,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz);
                        end
                        
                        if exist('coordinate_sbs','var') == 0%如果该值不存在
                            coordinate_sbs = [];
                        end
                        if exist('coordinate_addi_top_s_b_s','var') == 0%如果该值不存在
                            coordinate_addi_top_s_b_s = [];
                        end
                        coordinate_sbs_ahead = [coordinate_sbs; coordinate_addi_top_s_b_s];
                    end
                    
                    %                     if isempty(coordinate_addi_top_s_b_s) ==  0%说明其不为空
                    %                         if length(coordinate_addi_top_s_b_s) >= N_left%可以通过这个方式把零件排完
                    %                             N_left = 0;
                    %                             coordinate_sbs = coordinate_addi_top_s_b_s(1:N_left,:);
                    %                         else %不能通过这个方式把零件排完
                    %                             N_left = N_left - length(coordinate_addi_top_s_b_s);
                    %                             coordinate_sbs = coordinate_addi_top_s_b_s;
                    %                         end
                    %                     end
                    %                     [I_mix,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz);
                    
                    if N_left > 0
                        
                        [I_mix,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(I_mix,I2_t_max_or,N_left,sita_xz);
                        
                        if exist('coordinate_sbs_ahead','var')%如果该值存在
                            coordinate_sbs = [coordinate_sbs; coordinate_sbs_ahead]; %#ok<AGROW>
                        end
                        
                    else
                        coordinate_sbs = coordinate_sbs_ahead;
                    end
                    
                    if isempty(coordinate_sbs)%如果是空的
                        N_pai = 0;%p排的零件的数目
                    else
                        N_pai = length(coordinate_sbs(:,1));%
                    end
                    
                    if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                        nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                    end
                    
                    if ~isfield(data_list_o.result,'multiple')%如果multiple这个键不存在
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                        nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                    end
                    
                    if ~isfield(data_list_o.result.multiple', key )%如果这个键不存在
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                    end
                    
                    multiple_file_name = strcat('multiple_',nameahead,key,'R_',num_forname);
                    
                    eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = data_list_o.result.multiple{N_multiple,1}.' key '.partquantity + N_pai;'])
                    eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                    eval('data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;')
                    eval('data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;')%正确
                    ita = sum(I2_t0(:)) * eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                    eval(['data_list_o.result.multiple{N_multiple,1}.materialutilizationratio = ' num2str(ita) ';'])%正确
                    
                    save(strcat(outputpath,'\',multiple_file_name,'.txt'), 'coordinate_sbs', '-ascii')
                    
                    if board_is_full == 1 %说明板子排满了
                        N_multiple = N_multiple + 1;
                        nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        I_mix = heiban;
                    end
                end
                
            else%单排
                
                while N_left > 0
                    
                    %[I1_t,dx1_t,dy1_t,sita1_t] = find_I1_t(I2_t,I2_t_max_or,Parent0,4,template,key,maxIt_va,N_times,nPop,nVar,nC,n,maxIt);
                    [~,coordinate_sbs,board_is_full,N_left] = paitu_single_by_single(heiban,I2_t_max_or,N_left,sita_xz);
                    
                    if isempty(coordinate_sbs)%如果是空的
                        N_pai = 0;%p排的零件的数目
                    else
                        %N_left = N_left - length(coordinate_sbs(:,1));
                        N_left = mod(N_left,length(coordinate_sbs(:,1)));
                        N_pai = length(coordinate_sbs(:,1));
                    end
                    
                    if board_is_full == 1%板子满了，则认为是单排
                        
                        single_file_name = strcat('single_', key);
                        
                        if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                            eval(['data_list_o.result.single.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result,'single')%如果multiple这个键不存在
                            eval(['data_list_o.result.single.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result.single', key )%如果这个键不存在
                            eval(['data_list_o.result.single.' key '.partquantity = 0;'])
                        end
                        
                        eval(['data_list_o.result.single.' key '.partquantity = data_list_o.result.single.' key '.partquantity + N_pai;'])
                        eval(['data_list_o.result.single.' key '.txtfilepath = strcat(outputpath,"\",single_file_name,".txt");'])
                        eval(['data_list_o.result.single.' key '.sheetquantity = 1;'])
                        eval(['data_list_o.result.single.' key '.specifysheet = sheetsize;'])%正确
                        
                        ita = sum(I2_t0(:)) * eval(['data_list_o.result.single.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                        eval(['data_list_o.result.single.' key '.materialutilizationratio = ' num2str(ita) ';'])%正确
                        
                        save(strcat(outputpath,'\',single_file_name,'.txt'), 'coordinate_sbs', '-ascii')
                        %                     if board_is_full == 1 %说明板子排满了
                        %                         N_multiple = N_multiple + 1;
                        %                         nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        %                         I_mix = heiban;
                        
                    else%板子没有排满，则认为是混排
                        
                        if ~isfield(data_list_o,'result')%如果data_list_o.result这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result,'multiple')%如果multiple这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                            nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                        end
                        
                        if ~isfield(data_list_o.result.multiple', key )%如果这个键不存在
                            eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = 0;'])
                        end
                        
                        multiple_file_name = strcat('multiple_',nameahead,key,'R_',num_forname);
                        
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity = data_list_o.result.multiple{N_multiple,1}.' key '.partquantity + N_pai;'])
                        eval(['data_list_o.result.multiple{N_multiple,1}.' key '.txtfilepath = strcat(outputpath,"\",multiple_file_name,".txt");'])
                        eval('data_list_o.result.multiple{N_multiple,1}.sheetquantity = 1;')
                        eval('data_list_o.result.multiple{N_multiple,1}.specifysheet = sheetsize;')%正确
                        ita = sum(I2_t0(:)) * eval(['data_list_o.result.multiple{N_multiple,1}.' key '.partquantity']) / (numel(heiban(:,1)) * numel(heiban(1,:)));%正确
                        eval(['data_list_o.result.multiple{N_multiple,1}.materialutilizationratio = ' num2str(ita) ';'])%正确
                        
                        save(strcat(outputpath,'\',multiple_file_name,'.txt'), 'coordinate_sbs', '-ascii')
                        
                        N_multiple = N_multiple + 1;
                        nameahead = strcat(key,'_');%文件名称的前缀，为混排时候的第一个零件的名称
                    end
                end
            end
        end
    end
    
    data = [];
    finishpath = strcat(outputpath,'\','finish.txt');
    file = fopen(finishpath, 'w');
    fprintf(file, '%s\n', data);
    fclose(file);
    delete(gcf);delete(gcf),delete(gcf);
    output_args = 1;
    
    jsonStr = jsonencode(data_list_o);% 将 MATLAB 结构体编码为 JSON 字符串
    json_path_output = strcat(outputpath,'\',output_json_name);
    fid = fopen(json_path_output, 'w');% 将 JSON 字符串写入文件
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
    delete(gcf);delete(gcf),delete(gcf);
    fprintf('total time used: %f \n',toc(total_time));
    
end