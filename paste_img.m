function I_base = paste_img(I_base,I_paste,start_row,start_column)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% I_paste=trans2lefttop(I_paste);
[rows, columns, ~] = size(I_paste);
end_row = start_row + rows - 1;
end_column = start_column + columns - 1;
I_base(start_row:end_row, start_column:end_column, :) = I_paste;
% I_base = trans2mid(I_base);
end

