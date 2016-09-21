function Z1=netnetwon(X,Y,Z,X1,Y1)%X,Y,ZΪ����ֵ�����꣬��Ϊ���ӵ�
M=Z';
[n,m]=size(M);%ȡ�ò�ֵ����Ĵ�С
%���ȶ�x����newton����
x=Y(:,1)';
M1=zeros(n,m);
for i1=1:n
    L=M(i1,:);
    S=netwon(x,L);%��ÿһ�н���ţ�ٵ���
    M1(i1,:)=S;
end
%Ȼ���y�������netwon����
y=X(1,:);
M2=zeros(n,m);
for j1=1:m
    L1=M1(:,j1);
    L1=L1';
    S1=netwon(y,L1);
    M2(:,j1)=S1;
end
%����M2���е�ֵ������y����ֵ��ֵ
y1=X1(1,:);
n1=size(y1,2);
%���y������
A2=zeros(m,n1);
for j1=1:m
    L1=M2(:,j1);
    L1=L1';
    L2=netwondis(y1,L1,y);
    A2(j1,:)=L2;
end
%���㣨x,y������ֵ��;
x1=Y1(:,1)';
m1=size(x1,2);
%��ȡx������
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
x=x';y=y';%ת��
hangx=size(x,1);%��ȡx������;
H=ones(hangx,hangx);%������������x�ľ���������ȡϵ��
for j=1:hangx%���ѭ�����㷨���ģ����Ǹ������ǵı�
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
function Z=netwondis(z,y,x)%������̵�ֵ��zΪ����ֵ�㣬xΪԭ��ֵ�㣬yΪ����ֵ
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