I00=zeros(1200,2400);
I0=I00;
zb=textread('afile.txt');
% I2_0=imread('01.jpg');
% I3_0=imread('02.jpg');
% 
% I2_0=rgb2gray(I2_0);
% I3_0=rgb2gray(I3_0);
% 
% I2_0=im2bw(I2_0);
% I3_0=im2bw(I3_0);
% I2_0=~I2_0;
% I3_0=~I3_0;
% 
% I22=heiban;
% I33=heiban;
% 
% [m,n]=size(I2_0);
% for i=1:m
%     for j=1:n
%         I22(i,j)=I2_0(i,j);       
%         
%     end
% end

% [m,n]=size(I3_0);
% for i=1:m
%     for j=1:n
%         I33(i,j)=I3_0(i,j);       
%         
%     end
% end
I1_00=I2_0;
% I2_00=I2_0;
% I3_00=I2_0;
% I2_0=I22;
% I3_0=I33;

% for i=1:160   
%     I22=trans2mid(I2_0);
%     I22=imrotate(I22,zb2(i,3),'crop');
%     I22=trans2rightbuttom(I22);
%     I22=imtranslate(I22,[zb2(i,1), -zb2(i,2)]);
%     I0=I0 | I22;
% end
% 
% imshow(I0);

for i=1:572
    I11=trans2mid(I1_00);
    I11=imrotate(I11,zb(i,3),'crop','bicubic');
    I11=trans2rightbuttom(I11);
    I11=imtranslate(I11,[zb(i,1), -zb(i,2)]);
    I0=I0 | I11;
    imshow(I0);
end
% 
% for i=51:100 
%     I22=trans2mid(I2_00);
%     I22=imrotate(I22,zb(i,3),'crop');
%     I22=trans2rightbuttom(I22);
%     I22=imtranslate(I22,[zb(i,1), -zb(i,2)]);
%     I0=I0 | I22;
%     imshow(I0);
% end
% 
% 
% for i=101:150 
%     I33=trans2mid(I3_00);
%     I33=imrotate(I33,zb(i,3),'crop');
%     I33=trans2rightbuttom(I33);
%     I33=imtranslate(I33,[zb(i,1), -zb(i,2)]);
%     I0=I0 | I33;
%     imshow(I0);
% end
% 



