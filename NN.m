function Z1=netnetwon(X,Y,Z,X1,Y1)%X,Y,Z为待插值点坐标，均为格子点
M=Z';
[n,m]=size(M);%取得插值矩阵的大小
%首先对x方向newton递推
x=Y(:,1)';
M1=zeros(n,m);
for i1=1:n
    L=M(i1,:);
    S=netwon(x,L);%对每一行进行牛顿递推
    M1(i1,:)=S;
end
%然后对y方向进行netwon递推
y=X(1,:);
M2=zeros(n,m);
for j1=1:m
    L1=M1(:,j1);
    L1=L1';
    S1=netwon(y,L1);
    M2(:,j1)=S1;
end
%利用M2表中的值来计算y待插值点值
y1=X1(1,:);
n1=size(y1,2);
%获得y的坐标
A2=zeros(m,n1);
for j1=1:m
    L1=M2(:,j1);
    L1=L1';
    L2=netwondis(y1,L1,y);
    A2(j1,:)=L2;
end
%计算（x,y）待插值点;
x1=Y1(:,1)';
m1=size(x1,2);
%获取x的坐标
Z1=zeros(n1,m1);
for i1=1:n1
    l1=A2(:,i1);
    l1=l1';
    Z=netwondis(x1,l1,x);
    Z1(i1,:)=Z;
end
Z1=Z1';
end
function S=netwon(x,y)
x=x';y=y';%转致
hangx=size(x,1);%获取x的行数;
H=ones(hangx,hangx);%生成行数等于x的矩阵，用来存取系数
for j=1:hangx%这个循环是算法核心，即那个倒三角的表
    if(j==1)
        H(:,1)=y;
    else
    for i=j:hangx
        H(i,j)=(H(i,j-1)-H(i-1,j-1))/(x(i,1)-x(i-j+1,1));
    end
    end
end
S=zeros(1,hangx);
for i=1:hangx
    S(1,i)=H(i,i);
end
end
function Z=netwondis(z,y,x)%输出插商的值，z为待插值点，x为原插值点，y为插商值
hangz=size(z,1);
Z=ones(hangz,1);
x=x';
hangx=size(x,1);
for k=1:hangx
    if(k==1)
        Z(:,1)=y(1,1);
        S=z-x(1,1);
    else
        Z=Z+y(1,k)*S;
        S=S.*(z-x(k,1));
    end
end
 Z=Z';
end