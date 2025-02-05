clear
clc
tic
json_path_input='D:\input\0321.json';
%LAYOUTDESIGN 此处显示有关此函数的摘要
%   此处显示详细说明

jsonStr = fileread(json_path_input);% 读取JSON文件
data_list = jsondecode(jsonStr);% 解析JSON文件
data_list_o = jsondecode(jsonStr);%保持者
N_part = numel(fieldnames(data_list.part));
time=data_list.param.time;
toolszie=data_list.param.toolszie;
ignorehole=data_list.param.ignorehole;
partgap=data_list.param.partgap;
leftmargin=data_list.param.leftmargin;
rightmargin=data_list.param.rightmargin;
topmargin=data_list.param.topmargin;
downmargin=data_list.param.downmargin;
sheetsize=data_list.param.specifysheet;
enablenesting=data_list.param.enablenesting;
enablerevolve=data_list.param.enablerevolve;
accuracy=data_list.param.accuracy;
algorithmselect=data_list.param.algorithmselect;
outputpath=data_list.param.resultjsonpath;
outputpath = strsplit(outputpath, '\');
output_json_name=outputpath{end};
outputpath(end)=[];
outputpath=strjoin(outputpath, '\');
sheet_size = sort(str2double(strsplit(sheetsize, 'x')));
sheet_X=sheet_size(3);
sheet_Y=sheet_size(2);
% multiple_pngfile_name=data_list.result2.multiple.pngfilepath;
num_forname=num2str((str2double(regexp(output_json_name, '\d+', 'match'))));

% 获取结构体变量中所有的键值
keys = fieldnames(data_list.part);
s_filepath=struct();
L_sz=zeros(N_part,1);

% 遍历所有的键值
for i = 1:length(keys)
    key = keys{i};
    eval(['Path.' key ' =  data_list.part.' key '.filepath;']);
    eval(['X_length.' key ' =  data_list.part.' key '.length;']);
    eval(['Y_length.' key ' =  data_list.part.' key '.width;']);
    eval(['Nmax.' key ' =  data_list.part.' key '.quantity;']);
    eval(['img.' key '= {rgb2gray(imread(Path.' key '))};']);
    eval(['L_total.' key '= cell2mat({X_length.' key '+Y_length.' key '});']);
    eval(['L_sz(i) = L_total.' key ';']);
    eval(['alg1.' key '=  data_list.part.' key '.alg1;']);
    eval(['alg2.' key '=  data_list.part.' key '.alg2;']);
    eval(['kb.' key '=data_list.part.' key '.kb;']);
end

[L_sz, so]=sort(L_sz, 'descend');
key_seq=cell(N_part,1);
keys_t=keys;
for i=1:N_part
    for j=1:length(keys_t)
        key = keys_t{j};
        if eval(['L_sz(i) == L_total.' key])
            key_seq{i}=key;
            keys_t(j)=[];
            break
        end
    end
end

I0=zeros(sheet_Y-leftmargin-rightmargin,sheet_X-topmargin-downmargin);%空白板
heiban=logical(I0);

n=2;%图案的数量
nVar = 25*n; % x的长度
nPop0 = 20;
nPop = 20; % 种群规模大小
maxIt = 5; % 迭代次数
maxIt_va=1;%有的时候怕陷入局部最优，因此多试几次，能找到最优值
nPc = 1; % 交叉的比例
imrf=1;%缩放比例
N_times=8;%子代是父代的N_times倍

I0=imresize(I0,imrf);
I0=logical(I0);
mm_pix=(sheet_X-topmargin-downmargin)/numel(I0(1,:));%每个像素代表几mm
heiban=imresize(heiban,imrf);
I_finally=heiban;
coordinate_all=[];

for p=1:N_part %所排的图的个数
    key = key_seq{p};
    eval(['I2_tt=cell2mat(img.' key ');']);
    I2_0=heiban;
    I2=I2_tt;%这个纯属为了消除警告
    ave_I2=sum(I2(:))/(numel(I2(:,1))*numel(I2(1,:)));
    for i=1:numel(I2(:,1))%对I2进行二值化
        for j=1:numel(I2(1,:))
            if I2(i,j)>ave_I2
                I2(i,j)=255;
            else
                I2(i,j)=0;
            end
        end
    end
    
    I2=~I2;
    I2_x=l_x(I2);
    eval(['beishu=(X_length.' key '/I2_x)/mm_pix;']);
    
    I2=imresize(I2,beishu);
    
    for i=1:numel(I2(:,1))
        for j=1:numel(I2(1,:))
            I2_0(i,j)=I2(i,j);
        end
    end
    
    I2=I2_0;
    
    nC = round(nPop0 * nPc / 2) * 2; % 子代规模的大小
    nMu = 0.99; % 变异的概率
    [Parent0,template,I2_t,~] = Pretreat(I0,I2,nPop0);%图像预处理,这里的I2_t只是单纯的二值化和背景处理，Parent也是空的，没有赋值
    I2_t=~I2_t;
    I2_t=trans2mid(I2_t);
    [sita_best,~] = find_I2_t_xmin(I2_t);
    I2_t_ymin=imrotate(I2_t,sita_best-90,'crop','bicubic');
    I2_t_ymin=trans2mid(I2_t_ymin);
    
    lx1=l_x(I2_t);
    ly1=l_y(I2_t);
    
    I2_t=I2_t_ymin;%初始赋值
    I2_t0=I2_t_ymin;%这是为了保留原始图像的尺寸
    
    W0=min(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的宽度
    L0=max(numel(I2_t(:,1)),numel(I2_t(1,:)));%空白板的长度
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %第一次迭代寻优
    newRc=zeros(round(2.2*l_y(I2_t)),round(2.2*l_x(I2_t)));
    nl=round((numel(I2_t(1,:))-numel(newRc(1,:)))/2);
    nt=round((numel(I2_t(:,1))-numel(newRc(:,1)))/2);
    
    I2_t_temp=zeros(1,1);%初始化，不然会有警告
    for i=1:numel(newRc(:,1))
        for j=1:numel(newRc(1,:))
            I2_t_temp(i,j)=I2_t(i+nt,j+nl);
        end
    end
    
    I2_t=logical(I2_t_temp);%裁剪成newRc尺寸，newRc=zeros(min(2*lw,numel(I2_t(:,1))),min(2*lw,numel(I2_t(1,:))));
    imshow(I2_t);
    Parent=Parent0;
    flag=0;
    Alg1=eval(['alg1.' key]);
    while flag==0
        Grandma=template;
        for It = 1 : maxIt_va %迭代次数
            disp(['外层迭代次数：',num2str(It)]);
            Mother=repmat(template, 1, 1);
            Parent=Initial_pop(newRc,I2_t,Parent,nPop,nVar,n,lx1,ly1);%
            Offspring = repmat(Parent, N_times*3, 1);
            
            Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
            [Parent,Mother] = Up_pop_Mo1(Mother,Parent,Offspring,nPop,I2_t,newRc,lx1,ly1,n,template,Alg1);%0.5s 更新种群,或者说更新父辈
            
            Mother_last.x=0;
            Mother_last.y=0;
            Mother_last.z=0;
            
            It_Mo=0;%Mother迭代次数
            flag_y1=0;
            while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
                if isempty(Mother(1).x)==0%Mother被赋值
                    Offspring = repmat(Parent, N_times, 1);
                    Offspring(:)=Mother(end);
                    Parent(:)=Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                    Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                    [Parent,Mother] = Up_pop_Mo1(Mother,Parent(1),Offspring,nPop,I2_t,newRc,lx1,ly1,n,template,Alg1);%2s多 更新种群,或者说更新父辈
                else
                    while isempty(Mother(1).x)%循环，直至Mother被赋值，即出现完全分离的情况
                        Offspring = repmat(Parent, N_times*3, 1);
                        Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                        [Parent,Mother] = Up_pop_Mo1(Mother,Parent,Offspring,nPop,I2_t,newRc,lx1,ly1,n,template,Alg1);%2s多 更新种群,或者说更新父辈
                    end
                end
                
                if Mother(end).z~=Mother_last.z%Mother有更新
                    Mother_last=Mother(end);
                    I=Mother_last.tu;
                    imshow(I);
                    disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                    It_Mo=0;
                else%如果Mother没有更新
                    It_Mo=It_Mo+1;
                    disp(['Mother未更新次数：',num2str(It_Mo)]);
                end
            end
            
            if isempty(Grandma.x)==1 %如果Grandma是空的
                Grandma=Mother(end);
            else
                if Grandma.z<Mother(end).z
                    Grandma=Mother(end);
                end
            end
        end
        
        if It==maxIt_va
            flag=1;
        end
        
        dx1=Grandma.ddx;
        dy1=Grandma.ddy;
        sita1=Grandma.sita;
        I=Grandma.tu;
        I=trans2mid(I);
        imshow(I);
        disp('奶奶出现,第一次寻优完毕！');%得到目前最好的排布
        I1=I;%I1为第一次组合成功后的结果
    end
    
    I2_t_temp=zeros(W0,L0);%W0和L0为空白板的尺寸
    
    for i=1:numel(I1(:,1))
        for j=1:numel(I1(1,:))
            I2_t_temp(i,j)=I1(i,j);
        end
    end
    
    I2_t_temp=logical(I2_t_temp);
    
    I2_t=I2_t_temp;%这里的I2_t是已经排好的两个图元的组合
    I2_t=trans2leftbottom(I2_t);%将其移动到左下角
    I2_t_double=I2_t;%I2_t_double为两个69式排好的图元，为后期的X方向的寻优做准备。
    
    I2_t2=I2_t;
    I2_t2=trans2mid(I2_t2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %下面开始高度方向上的自循环，这是第二次迭代寻优
    eval(['newRc=zeros(round(2.5*l_y(I2_t2)),round(kb.' key '*l_x(I2_t2)));']);
    nt=round((numel(I2_t2(:,1))-numel(newRc(:,1)))/2);
    nl=round((numel(I2_t2(1,:))-numel(newRc(1,:)))/2);
    
    I2_t2_temp=zeros(1,1);%初始化，不然会有警告
    
    for i=1:numel(newRc(:,1))
        for j=1:numel(newRc(1,:))
            I2_t2_temp(i,j)=I2_t2(nt+i,j+nl);
        end
    end
    I2_t2_temp=logical(I2_t2_temp);
    
    lx2=l_x(I2_t2_temp);
    ly2=l_y(I2_t2_temp);
    
    I2_t2=I2_t2_temp;%裁剪成newRc尺寸，newRc=zeros(min(2*lw,numel(I2_t(:,1))),min(2*lw,numel(I2_t(1,:))));
    imshow(I2_t2)
    
    Parent=Parent0;%Parent0为空的
    flag=0;
    Alg2=eval(['alg2.' key]);
    while flag==0
        Grandma=template;
        for It = 1 : maxIt_va %迭代次数
            disp(['外层迭代次数：',num2str(It)]);
            Mother=repmat(template, 1, 1);
            Parent=Initial_pop(newRc,I2_t2,Parent,nPop,nVar,n,lx2,ly2);%2s多   生成子代种群（包含交叉操作）
            Offspring = repmat(Parent, N_times, 1);
            Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
            [Parent,Mother] = Up_pop_Mo2(Mother,Parent,Offspring,nPop,I2_t2,newRc,lx2,ly2,n,template,Alg2);%0.5s 更新种群,或者说更新父辈
            
            Mother_last.x=0;
            Mother_last.y=0;
            Mother_last.z=0;
            
            It_Mo=0;%Mother迭代次数
            flag_y1=0;
            while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
                if isempty(Mother(1).x)==0%Mother被赋值
                    Offspring = repmat(Parent, N_times, 1);
                    Offspring(:)=Mother(end);
                    Parent(:)=Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                    Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                    [Parent,Mother] = Up_pop_Mo2(Mother,Parent(1),Offspring,nPop,I2_t2,newRc,lx2,ly2,n,template,Alg2);%2s多 更新种群,或者说更新父辈
                else
                    while isempty(Mother(end).x)%循环，直至Mother被赋值，即出现完全分离的情况
                        Offspring = repmat(Parent, N_times, 1);
                        Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                        [Parent,Mother] = Up_pop_Mo2(Mother,Parent,Offspring,nPop,I2_t2,newRc,lx2,ly2,n,template,Alg2);%2s多 更新种群,或者说更新父辈
                    end
                end
                
                if Mother(end).z~=Mother_last.z%Mother有更新
                    Mother_last=Mother(end);
                    I=Mother_last.tu;
                    imshow(I);
                    disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                    It_Mo=0;
                else%如果Mother没有更新
                    It_Mo=It_Mo+1;
                    disp(['Mother未更新次数：',num2str(It_Mo)]);
                end
            end
            
            if isempty(Grandma(1).x)==1 %如果Grandma是空的
                Grandma=Mother(end);
            else
                if Grandma(1).z<Mother(end).z
                    Grandma=Mother(end);
                end
            end
        end
        
        if It==maxIt_va
            flag=1;
        end
        
        I2=Grandma.tu;
        dx2=Grandma.ddx;
        dy2=Grandma.ddy;
        
        I2=trans2mid(I2);
        imshow(I2);
        disp('奶奶出现,第二次寻优完毕！');%I2为第二次组合成功后的结果
    end
    
    ZP=0;%左偏
    YP=0;%右偏
    
    if dy2>=0 && dx2>=0
        ZP=1;
    elseif dy2>=0 && dx2<0
        YP=1;
    end
    
    if dy2<0 && dx2>=0
        YP=1;
    elseif dy2<0 && dx2<0
        ZP=1;
    end
    
    if ZP==1
        dx2=-abs(dx2);
        dy2=-abs(dy2);
    else
        dx2=abs(dx2);
        dy2=-abs(dy2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %以上完成了Y方向的排布，下面需要进行X方向的排布,这是第三次迭代
    if ZP==1
        sita2=-atan(abs(dx2/dy2));
        I_y5=trans2rightbottom(I2_t_double);
    else
        sita2=atan(abs(dx2/dy2));
        I_y5=trans2leftbottom(I2_t_double);
    end
    
    sita0=atan(-dx2/dy2)*180/3.1415926535;%最终结果180度数制
    
    I_y0=I_y5;
    for i=1:4
        I_y5=I_y5 | imtranslate(I_y0,[i*dx2, i*dy2]);
    end
    
    I_y5=trans2mid(I_y5);
    LY_I_y5=l_y(I_y5);
    LX_I_y5=l_x(I_y5);
    I_y5_s=imrotate(I_y5,sita0,'crop','bicubic');%I_y5_temp为将I_y5搬竖直了
    lx_I_y5_s=l_x(I_y5_s);%搬竖直后横向尺寸
    lx_I_y5=lx_I_y5_s/cos(sita2);
    
    ls=min(numel(I_y5(:,1)),round(LY_I_y5*1.25));
    lh=min(numel(I_y5(1,:)),round(LX_I_y5+lx_I_y5*1.25));
    
    newRc=zeros(ls,lh);%横竖正确
    nt=round((numel(I_y5(:,1))-numel(newRc(:,1)))/2);
    nl=round((numel(I_y5(1,:))-numel(newRc(1,:)))/2);
    
    I_y5_temp=zeros(1,1);%初始化，不然会有警告
    
    for i=1:numel(newRc(:,1))
        for j=1:numel(newRc(1,:))
            I_y5_temp(i,j)=I_y5(nt+i,j+nl);%到了这里，I_y5_temp依然是倾斜的角度，只是将其裁剪成newRc的大小
        end
    end
    
    I_y5_temp=logical(I_y5_temp);
    lx3=l_x(I_y5_temp);
    ly3=l_y(I_y5_temp);
    I_y5=I_y5_temp;%裁剪成newRc尺寸，newRc=zeros(min(2*lw,numel(I2_t(:,1))),min(2*lw,numel(I2_t(1,:))));
    imshow(I_y5)
    
    Parent=Parent0;
    flag=0;
    while flag==0
        Grandma=repmat(template, 1, 1);
        for It = 1 : maxIt_va %迭代次数
            disp(['外层迭代次数：',num2str(It)]);
            Mother=template;
            Parent=Initial_pop(newRc,I_y5,Parent,nPop,nVar,n,lx3,ly3);%2s多   生成子代种群（包含交叉操作）
            Offspring = repmat(Parent, N_times, 1);
            Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
            [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I_y5,newRc,lx3,ly3,n,template,sita0);%0.5s 更新种群,或者说更新父辈
            
            Mother_last.x=0;
            Mother_last.y=0;
            Mother_last.z=0;
            
            It_Mo=0;%Mother迭代次数
            flag_y1=0;
            while It_Mo<maxIt%如果Mother超过maxIt次没有找到更好的解，则认为达到最优
                if isempty(Mother(1).x)==0%Mother被赋值
                    Offspring = repmat(Parent, N_times, 1);
                    Offspring(:)=Mother(end);
                    Parent(:)=Mother(end);%因为目前看的话Mother中的元素最优，所以直接Mother中的元素给Parent和Offspring进行赋值
                    Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                    [Parent,Mother] = Up_pop_Mo3(Mother,Parent(1),Offspring,nPop,I_y5,newRc,lx3,ly3,n,template,sita0);%2s多 更新种群,或者说更新父辈
                else
                    while isempty(Mother(1).x)%循环，直至Mother被赋值，即出现完全分离的情况
                        Offspring = repmat(Parent, N_times, 1);
                        Offspring = Variation(Offspring,nC,n,newRc,nMu);%变异操作
                        [Parent,Mother] = Up_pop_Mo3(Mother,Parent,Offspring,nPop,I_y5,newRc,lx3,ly3,n,template,sita0);%2s多 更新种群,或者说更新父辈
                    end
                end
                
                if Mother(end).z~=Mother_last.z%Mother有更新
                    Mother_last=Mother(end);
                    %                     [I,~,~] = Nw_cal_xy(newRc,Mother,I_y5,n,lx3,ly3);%0.07s 选取方差最小的元素，将独立的经过解码后的图像合并在一张图I上
                    I=Mother_last.tu;
                    imshow(I);
                    disp(['Mother中的元素数目：',num2str(numel(Mother))]);
                    It_Mo=0;
                else%如果Mother没有更新
                    It_Mo=It_Mo+1;
                    disp(['Mother未更新次数：',num2str(It_Mo)]);
                end
            end
            
            if isempty(Grandma(1).x)==1 %如果Grandma是空的
                Grandma=Mother(end);
            else
                if Grandma(1).z<Mother(end).z
                    Grandma=Mother(end);
                end
            end
        end
        
        if It==maxIt_va
            flag=1;
        end
        
        I3=Grandma.tu;
        dx3=Grandma.ddx;
        dy3=Grandma.ddy;
        I3=trans2mid(I3);
        imshow(I3);
        disp('奶奶出现,第三次寻优完毕！');%第三次寻优结果
    end
    
    SP=0;%上偏
    XP=0;%下偏
    
    if dx3>=0 && dy3>=0
        XP=1;
    elseif dx3>=0 && dy3<0
        SP=1;
    end
    
    if dx3<0 && dy3>=0
        SP=1;
    elseif dx3<0 && dy3<0
        XP=1;
    end
    
    if SP==1
        dx3=abs(dx3);
        dy3=-abs(dy3);
    else
        dx3=abs(dx3);
        dy3=abs(dy3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %以上完成了第三次迭代，下面进行排图
    bujin=10;
    lx_single=l_x(I2_t0);
    
    N_shu=floor(numel(heiban(:,1))/abs(dy2));  %结果正确 其表示纵向最多能排多少个图案
    eval(['N_H_initial=floor(Nmax.' key '/(2*N_shu))*l_x(I1);'])%有问题    初始板的横向像素数目
    N_H_initial=min(N_H_initial,numel(heiban(1,:)));
    
    I2_s=imrotate(I2,sita2*180/3.1415926535,'crop','bicubic');%将I2搬竖直了
    I3_s=imrotate(I3,sita2*180/3.1415926535,'crop','bicubic');%将I3搬竖直了
    ddx3=(l_x(I3_s)-l_x(I2_s))/cos(sita2);%ddx3为第三次排图过程中两列之间的x方向的距离
    delta_xx=abs(numel(heiban(:,1))*tan(sita2));%斜铺后横向多出来得到部分像素（实际上不存在）
    x_total=delta_xx+numel(heiban(1,:));
    N_dx3_max=ceil(x_total/ddx3);%横向做能排多少个
    
    N_flag=0;%其为排图的数目达到要求的标记
    single=1;%单排的标记
    initial_flag1=0;
    initial_flag2=0;
    initial_flag3=0;
    d3_nitial=0;%dx3和dy3是否被初始调整过的标记
    while eval(['Nmax.' key '>0' ])% && flag_N_h==0'])
        %while eval(['Nmax.' key '>0 && (N<Nmax.' key ') && flag_N_h==0'])
        if single == 1 %用initial_flag1作为标记
            if initial_flag1==0
                initial_flag1=1;
                n_p=1;%铺图次数
                eval(['N_H_initial=floor(Nmax.' key '/(2*N_shu))*l_x(I1);'])%有问题  初始板的横向像素数目
                N_H_initial=min(N_H_initial,numel(heiban(1,:)));
                N_H_initial=max(N_H_initial,l_x(I1));
                I0=zeros(numel(heiban(:,1)),N_H_initial);
                I_mix_single=I0;%单个图形排布的I_mix
                I_mix_temp=I0;
                N_h=N_H_initial;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            else
                n_p=n_p+1;
                N_h=N_h+bujin;
                N_h=min(N_h,numel(heiban(1,:)));
                I0=zeros(numel(heiban(:,1)),N_h);
                I_mix_temp=I0;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            end
        end
        
        if single==0 && p==1 %用initial_flag2作为标记
            if initial_flag2==0
                initial_flag2=1;
                n_p=1;%铺图次数
                eval(['N_H_initial=floor(Nmax.' key '/(2*N_shu))*l_x(I1);'])%有问题  初始板的横向像素数目
                N_H_initial=min(N_H_initial,numel(heiban(1,:)));
                N_H_initial=max(N_H_initial,l_x(I1));
                I0=zeros(numel(heiban(:,1)),N_H_initial);
                I_mix_multiple=I0;
                I_mix_temp=I0;
                N_h=N_H_initial;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            else
                n_p=n_p+1;
                N_h=N_h+bujin;
                N_h=min(numel(heiban(1,:)),N_h);
                I0=zeros(numel(heiban(:,1)),N_h);
                I_mix_temp=I0;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            end
        end
        
        if single== 0 && p>1 %用initial_flag3作为标记
            if initial_flag3==0
                initial_flag3=1;
                n_p=1;%铺图次数
                eval(['N_H_initial=floor(Nmax.' key '/(2*N_shu))*l_x(I1);'])
                N_h=min(numel(heiban(1,:)),(numel(I_mix_multiple(1,:))+N_H_initial));%混排，且不是第一个零件，因此，板子的初始宽度因为前期混排好的板子
                I0=zeros(numel(heiban(:,1)),N_h);
                I_mix_temp=I0;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            else
                n_p=n_p+1;
                N_h=N_h+bujin;
                N_h=min(numel(heiban(1,:)),N_h);
                I0=zeros(numel(heiban(:,1)),N_h);
                I_mix_temp=I0;
                if N_h<numel(heiban(1,:))
                    flag_N_h=0;
                else
                    flag_N_h=1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %上面是调整板材尺寸
        %下面是真正的着手排布零件
        if single==1
            I_mix=I_mix_single;
        else
            I_mix=I_mix_multiple;
        end
        
        if ZP==1%如果是左偏，则往右边加空白
            for k=1:numel(I_mix(1,:))
                I_mix_temp(:,k)=I_mix(:,k);
            end
        else%如果是右偏，则往左边加空白
            delta_x=numel(I0(1,:))-numel(I_mix(1,:));
            for k=1:numel(I_mix(1,:))
                I_mix_temp(:,delta_x+k)=I_mix(:,k);
            end
        end
        
        I_mix=logical(I_mix_temp);
        
        if single==1
            I_mix_single=I_mix;
        else
            I_mix_multiple=I_mix;
        end
        
        I2_t_temp0=getmid(I0,I2_t0); %取其中的图案,是单个零件的图案,最终的大小为I0的大小
        I2_t_temp0=trans2mid(I2_t_temp0);
        I2_t_temp0=logical(I2_t_temp0);
        
        if ZP==1 %左偏
            if n_p==1 && d3_nitial==0
                d3_nitial=1;
                if dx1<0
                    dx1=-dx1;
                    dy1=-dy1;
                    sita1=sita1+180;
                end
                
                if dy3>0%说明其下偏，但是不知道偏的是不是合适
                    while dy3>0
                        dy3=dy3+dy2;
                        dx3=dx3+dx2;
                    end
                    dy3=dy3-dy2;%保证右边是下偏最小的距离
                    dx3=dx3-dx2;
                elseif dy3<0
                    while dy3<0%说明其上偏，必须让其下偏，这么做是为了在向右铺的时候，右边的最下边在分界线一下，不至于浪费空间
                        dy3=dy3-dy2;
                        dx3=dx3-dx2;
                    end
                end
            end
            I_s_m=imrotate(I2_t_temp0,sita1,'crop','bicubic');
            I_s_m=trans2leftbottom(I_s_m);%单个移动图元至左下角，此时还不知道dy1是不是大于0，如果dy1大于0，说明其下偏，那么I_s_m_r显示不完全
            I_s_m_r=imrotate(I2_t_temp0,sita1+180,'crop','bicubic');%单个移动图旋转180
            I_s_m_r=trans2leftbottom(I_s_m_r);
            
            if dy1>0
                I_s_m=imtranslate(I_s_m,[0,-dy1]);
                I_s_m_r=imtranslate(I_s_m_r,[0,-dy1]);
            end
            I_s_m_r=imtranslate(I_s_m_r,[dx1,dy1]);
            I_s_m0=I_s_m;
            I_accm=I_s_m0;%初始化
            I_s_m_r0=I_s_m_r;
            I_accm_r=I_s_m_r0;%初始化
            
        else%右偏
            
            if n_p==1 && d3_nitial==0
                d3_nitial=1;
                if dx1>0
                    dx1=-dx1;
                    dy1=-dy1;
                    sita1=sita1+180;
                end
                %dx3和dy3分别为右边的块相对于左边的块的相对位移，现在是从右往左铺，需要知道的是左边的快相对于右边的块的位移
                %因此需要进行如下取相反数的操作
                dx3=-dx3;
                dy3=-dy3;
                
                if dy3>0%说明其下偏，但是不知道偏的是不是合适
                    while dy3>0
                        dy3=dy3+dy2;
                        dx3=dx3+dx2;
                    end
                    dy3=dy3-dy2;%保证右边是下偏最小的距离
                    dx3=dx3-dx2;
                elseif dy3<0
                    while dy3<0%说明其上偏，必须使得其下偏
                        dy3=dy3-dy2;
                        dx3=dx3-dx2;
                    end
                end
            end
            I_s_m=imrotate(I2_t_temp0,sita1,'crop','bicubic');
            I_s_m=trans2rightbottom(I_s_m);%单个移动图元至右下角
            I_s_m_r=imrotate(I2_t_temp0,sita1+180,'crop','bicubic');%单个移动图旋转180
            I_s_m_r=trans2rightbottom(I_s_m_r);
            
            if dy1>0
                I_s_m=imtranslate(I_s_m,[0,-dy1]);
                I_s_m_r=imtranslate(I_s_m_r,[0,-dy1]);
            end
            I_s_m_r=imtranslate(I_s_m_r,[dx1,dy1]);
            I_s_m0=I_s_m;
            I_accm=I_s_m0;%初始化
            I_s_m_r0=I_s_m_r;
            I_accm_r=I_s_m_r0;%初始化
        end
        
        if single ==1 %单一的零件排布
            if flag_N_h==1 %横向尺寸达到最大了，下次就不是单一的排布了，而是跟其他零件混排
                single=0;%意为，下次再排的话就不是一种零件形状了
            end
            [coordinate1,coordinate2,I_mix,N_flag,N]=paitu(ZP,n_p,dx1,dx2,dx3,dy1,dy2,dy3,sita1,Nmax,key,...
                N_flag,N_dx3_max,N_shu,I0,bujin,I_mix,0,I_s_m0,I_s_m_r0,I_accm,I_accm_r,lx_single);
            I_mix_single=I_mix;
            eval(['Nmax.' key '=Nmax.' key '-N;']);
            partquantity=N;%一个板子上排的零件的数目
            if eval(['Nmax.' key])==0 || flag_N_h==1
                
                single_file_name = strcat('single_', key);
                single_file_name=strcat(outputpath,'\',single_file_name,'R_',num_forname);
                eval(['data_list_o.result.single.' key '.locationfilepath=strcat(single_file_name,".txt");']);
                eval(['data_list_o.result.single.' key '.pngfilepath=strcat(single_file_name,".png");']);
                eval(['data_list_o.result.single.' key '.partquantity=' num2str(partquantity) ';']);%正确
                eval(['sheetquantity=floor((data_list_o.part.' key '.quantity)/N);']);%正确
                eval(['data_list_o.result.single.' key '.sheetquantity=' num2str(sheetquantity) ';']);%正确
                eval(['data_list_o.result.single.' key '.specifysheet=data_list.param.specifysheet;']);%正确
                ita=sum(I_mix(:))/(numel(I_mix(:,1))*numel(I_mix(1,:)));%正确
                eval(['data_list_o.result.single.' key '.materialutilizationratio=' num2str(ita) ';']);%正确
                %eval(['data_list_o.result.single.' key '.sheetcenter_X=num2str(1220);']);%正确
                %eval(['data_list_o.result.single.' key '.sheetcenter_Y=num2str(610);']);%正确
                
                I_last=heiban;
                for i=1:numel(I_mix(1,:))
                    I_last(:,i)=I_mix(:,i);
                end
                
                I_last=trans2mid(I_last);
                imwrite(I_last,strcat(single_file_name,'.png'),'png');
                
                coordinate=[coordinate1;coordinate2];
                [nl0,nr0,nt0,nb0] = margin(I_last);
                if ZP ==1
                    coordinate(:,1)=coordinate(:,1)+l_x(I_s_m0)/2+(nr0-nl0)/2;
                else
                    coordinate(:,1)=coordinate(:,1)+numel(I0(1,:))/2-l_x(I_s_m0)/2+(nr0-nl0)/2;
                end
                coordinate(1,:)=coordinate(1,:)+(nb0-nt0)/2;
                save(strcat(single_file_name,'.txt'), 'coordinate', '-ascii')
                if eval(['N<Nmax.' key])
                    eval(['Nmax.' key '=mod(Nmax.' key ',N);']);
                else
                    eval(['Nmax.' key '=0;']);
                end
            end
            
        else %混排
            
            [coordinate1,coordinate2,I_mix,N_flag,N]=paitu(ZP,n_p,dx1,dx2,dx3,dy1,dy2,dy3,sita1,Nmax,key,...
                N_flag,N_dx3_max,N_shu,I0,bujin,I_mix,0,I_s_m0,I_s_m_r0,I_accm,I_accm_r,lx_single);
            I_mix_multiple=I_mix;
            eval(['Nmax.' key '=Nmax.' key '-N;']);
            if ~isfield(data_list_o.result,'multiple')
                eval(['data_list_o.result.multiple.' key '.partquantity=0;']);
            end
            if isfield(data_list_o.result,'multiple')
                if ~isfield(data_list_o.result.multiple, key)
                    eval(['data_list_o.result.multiple.' key '.partquantity=0;']);
                end
            end
            if ~eval(['isfield(data_list_o.result.multiple.' key ',"partquantity")'])
                eval(['data_list_o.result.multiple.' key '.partquantity=0;'])
            end
            eval(['data_list_o.result.multiple.' key '.partquantity= data_list_o.result.multiple.' key '.partquantity+N;']);
            
            if eval(['Nmax.' key '== 0'])  || flag_N_h  %后期需要修改,因为并未考虑混排一个板子不够的情况
                
                multiple_file_name = strcat('multiple_', key,'R_',num_forname);
                multiple_file_name=strcat(outputpath,'\',multiple_file_name);
                eval(['data_list_o.result.multiple.' key '.locationfilepath=strcat(multiple_file_name,".txt");']);
                coordinate=[coordinate1;coordinate2];
                save(strcat(multiple_file_name,'.txt'), 'coordinate', '-ascii')
                
                break
            end
        end
    end
    
    
    
    %     I_last=heiban;
    %
    %     for i=1:numel(I_mix(1,:))
    %         I_last(:,i)=I_mix(:,i);
    %     end
    
    %     eval(['DX.' key '.dx1=dx1;']);
    %     eval(['DX.' key '.dx2=dx2;']);
    %     eval(['DX.' key '.dx3=dx3;']);
    %
    %     eval(['DY.' key '.dy1=dy1;']);
    %     eval(['DY.' key '.dy2=dy2;']);
    %     eval(['DY.' key '.dy3=dy3;']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %上面完成排图，下面开始保存结果
    %         ita=sum(I_mix(:))/(numel(I_mix(:,1))*numel(I_mix(1,:)));
    %         imshow(I_mix)
    %
    %         if ZP==1%左偏
    %             coordinate2(:,1)=coordinate2(:,1)+abs(dx1);
    %             coordinate2(:,2)=coordinate2(:,2)+abs(dy1);
    %         else%右偏
    %             coordinate2(:,1)=coordinate2(:,1)-abs(dx1);
    %             coordinate2(:,2)=coordinate2(:,2)+abs(dy1);
    %         end
    %
    %         coordinate=[coordinate1;coordinate2];
    %         coordinate_all=[coordinate_all;coordinate]; %#ok<*AGROW>
    %
    %         coordinate_all(:,1)=coordinate_all(:,1)-(numel(heiban(1,:))-numel(I_mix(1,:)));
    %
    %         coordinate(:,1)=coordinate(:,1)/imrf;
    %         coordinate(:,2)=coordinate(:,2)/imrf;
    %
    %         I_last=heiban;
    %
    %         for i=1:numel(I_mix(1,:))
    %             I_last(:,i)=I_mix(:,i);
    %         end
    %
    %         imshow(I_last)
    %
    %         if single==1
    %             img_name = strcat('single_', key, '.jpg');
    %             imwrite(I_last,img_name,'jpg');
    %         end
    %
    %         eval(['save ' file_name '.txt -ascii coordinate_all'])
end

I_last=I_mix;
if isfield(data_list_o.result,'multiple')
    imwrite(I_last,strcat(multiple_file_name,'.png'),'png');
    data_list_o.result.multiple.pngfilepath=multiple_file_name;
    data_list_o.result.multiple.sheetquantity=1;
    data_list_o.result.multiple.specifysheet=sheetsize;
    data_list_o.result.multiple.materialutilizationratio=sum(I_last(:))/(numel(I_last(:,1))*numel(I_last(1,:)));
end

jsonStr = jsonencode(data_list_o);% 将 MATLAB 结构体编码为 JSON 字符串
json_path_output=strcat(outputpath,'\',output_json_name);
fid = fopen(json_path_output, 'w');% 将 JSON 字符串写入文件
fwrite(fid, jsonStr, 'char');
fclose(fid);

data = [];
finishpath=strcat(outputpath,'\','finish.txt');
file = fopen(finishpath, 'w');
fprintf(file, '%s\n', data);
fclose(file);

close all
output_args=1;
toc
