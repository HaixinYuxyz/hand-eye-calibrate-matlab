function X = tsai(A,B)
% Calculates the least squares solution of
% AX = XB
% 
% A New Technique for Fully Autonomous 
% and Efficient 3D Robotics Hand/Eye Calibration
% Lenz Tsai
%
% Mili Shah
% July 2014

[m,n] = size(A); n = n/4;
S = zeros(3*n,3);
v = zeros(3*n,1);
%Calculate best rotation R
for i = 1:n
    %取出矩阵A和矩阵B的前半部分3*3的矩阵
    A1 = logm(A(1:3,4*i-3:4*i-1)); 
    B1 = logm(B(1:3,4*i-3:4*i-1));
    %罗德里格斯变换并归一化
    a = [A1(3,2) A1(1,3) A1(2,1)]'; a = a/norm(a);
    b = [B1(3,2) B1(1,3) B1(2,1)]'; b = b/norm(b);
    %计算初始向量
    %由于Skew(Pgij + Pcij)总是奇异的，所以至少要两组数据，也就是测三次
    S(3*i-2:3*i,:) = skew(a+b);
    v(3*i-2:3*i,:) = a-b;
end
%矩阵左除，得到的结果是初始旋转向量
x = S\v;
theta = 2*atan(norm(x));
x = x/norm(x);
R = (eye(3)*cos(theta) + sin(theta)*skew(x) + (1-cos(theta))*x*x')';
%Calculate best translation t
C = zeros(3*n,3);
d = zeros(3*n,1);
I = eye(3);
for i = 1:n
    %这是平移向量的前面的系数,C的维度是（3n,3）
    C(3*i-2:3*i,:) = I - A(1:3,4*i-3:4*i-1);
    %这是平移向量和系数乘积的结果C的维度是（3n,1）
    d(3*i-2:3*i,:) = A(1:3,4*i)-R*B(1:3,4*i);
end
%t的维度是(3,1)
t = C\d;
%Put everything together to form X
X = [R t;0 0 0 1];