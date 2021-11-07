function [K,T] = KTMat(Ele,No)

%相对介电常数
Eps = 2.2;

%获取元素数量和节点数
E_num = size(Ele,1);
N_num = size(No,1);

%获取节点坐标
x = No(:,1);
y = No(:,2);

%创建节点及其abc值矩阵
a = zeros(N_num,3);
b = zeros(N_num,3);
c = zeros(N_num,3);

%创建三角元面积矩阵
Area = zeros(E_num,1);

%创建临时矩阵
temp = zeros(2,3);

%创建三角元系数矩阵
K_all = zeros(N_num);
T_all = zeros(N_num);

for num = 1:E_num
    %获取全局节点编号
    [i,j,m] = deal(Ele(num,1),Ele(num,2),Ele(num,3));
    %获取编号对应的坐标值
    [xi,xj,xm] = deal(x(i),x(j),x(m));
    [yi,yj,ym] = deal(y(i),y(j),y(m));
    %获取a(i,j,m)、b(i,j,m)、c(i,j,m)值
    [a(num,1),a(num,2),a(num,3)] = deal(xj*ym - xm*yj,xm*yi - xi*ym,xi*yj - xj*yi);
    [b(num,1),b(num,2),b(num,3)] = deal(yj-ym,ym-yi,yi-yj);
    [c(num,1),c(num,2),c(num,3)] = deal(xm-xj,xi-xm,xj-xi);
    %获取每个三角元的面积
    Area(num) = 0.5*(b(num,1)*c(num,2)- b(num,2)*c(num,1));
end

%清除变量，后面仍会使用
clear num i j m;

for num = 1:E_num
    %获取bi、bj、bm;ci、cj、cm
    temp(1,:) = b(num,:);
    temp(2,:) = c(num,:);
    %求每个三角元系数矩阵
    Ke = ((temp')*temp).*(0.25*Eps/Area(num));
    %获取全局节点编号
    [i,j,m] = deal(Ele(num,1),Ele(num,2),Ele(num,3));
    %拓展为全局系数矩阵
    K_all(i,i) = K_all(i,i) + Ke(1,1);
    K_all(i,j) = K_all(i,j) + Ke(1,2);
    K_all(j,i) = K_all(j,i) + Ke(2,1);
    K_all(i,m) = K_all(i,m) + Ke(1,3);
    K_all(m,i) = K_all(m,i) + Ke(3,1);
    K_all(j,j) = K_all(j,j) + Ke(2,2);
    K_all(j,m) = K_all(j,m) + Ke(2,3);
    K_all(m,j) = K_all(m,j) + Ke(3,2);
    K_all(m,m) = K_all(m,m) + Ke(3,3);
    %求解T矩阵，并拓展为全局编号
    Te = [1/6,1/12,1/12;1/12,1/6,1/12;1/12,1/12,1/6].*Area(num).*Eps;
    T_all(i,i) = T_all(i,i) + Te(1,1);
    T_all(i,j) = T_all(i,j) + Te(1,2);
    T_all(j,i) = T_all(j,i) + Te(2,1);
    T_all(i,m) = T_all(i,m) + Te(1,3);
    T_all(m,i) = T_all(m,i) + Te(3,1);
    T_all(j,j) = T_all(j,j) + Te(2,2);
    T_all(j,m) = T_all(j,m) + Te(2,3);
    T_all(m,j) = T_all(m,j) + Te(3,2);
    T_all(m,m) = T_all(m,m) + Te(3,3);
end

K = K_all;
T = T_all;

end