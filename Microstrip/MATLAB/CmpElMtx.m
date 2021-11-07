function Ke = CmpElMtx(xy)

eps0 = 8.854e-12;

if (sum(xy(:,1))/3<=2)&&(sum(xy(:,2))/3<=0.25)
    epsr = 2.2;
else
    epsr = 1;
end

% 三角形的三边矢量
s1 = xy(3,:)-xy(2,:);
s2 = xy(1,:)-xy(3,:);
s3 = xy(2,:)-xy(1,:);

Atot = 0.5*(s2(1)*s3(2)-s2(2)*s3(1)); %计算三角形单元面积，矢量叉乘
% 检查面积是否为负
if (Atot<0)
    error('Wrong Oder!');
end

% 计算中间矩阵Se
Se1 = [-s1(2);s1(1)];
Se2 = [-s2(2);s2(1)];
Se3 = [-s3(2);s3(1)];
Se = [Se1 Se2 Se3];

% 计算该三角形单元的矩阵系数
for i = 1:3
    for j = 1:3
        Ke(i,j) = epsr*eps0*Se(:,i)'*Se(:,j)/(4*Atot);
    end
end

