function [y1, y2] = crossPop(x1, x2) %�������

    n = numel(x1);%n = numel(A) ��������A��Ԫ�ظ�����
    s = randi([2, n - 1]);% r = randi([imin,imax],...)  ����һ����[imin,imax]��Χ�ڵ�α�������

    y1 = [x1(1 : s) x2(s + 1: end)];
    y2 = [x2(1 : s) x1(s + 1: end)];

end