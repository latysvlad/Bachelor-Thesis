function Res = pic2val( P, tol, stepval, extr, Q)
%pic2val ������� ����������� ������� ������� � ������ ������
%   ������� ������������� ���� �������(��� ����). ������ - ������, 
%   ��� - �������. 
%   � ����� ����������� � ������� ������������ ���������(���, �����, 
%   ������ �������) ����������� ������������� ������������ ��������� 
%   ��� �������� ���� �������� (������ ����������� ��������� ���� 
%   ������� �� ��������� tol). 
%   ������ ����������� ������� � ������ ������� �������� ����������� ������
%   ���������� ������ � �������������� ������� �������.
%   ������� ����������:
%   P - �����-����� �����������.
%   tol - ������(%). ���������, �� ������� ������ ����� ����������
%   ��������� �������� �� �����������. �� ��������� -  2%.
%   stepval - ��������� �������� �������������(%). ���������, � ������
%   �������� ����� ������� ���� ������. �� ��������� -  50%.
%   extr - ������� �������������. ����� ���� �������. ���� ����� 
%   (n-round(extr)) � (n-1), ������ ������, �� ���� ������ ����� n-� �����.
%   ������� ����� �������� ������� ������������� ������������������
%   �������, ���������� ����, � ������� � ������� ������� �������� �
%   ���������� �����. �� ��������� -  1, �.�. ������� ������ �����������
%   ��������.
%   Q - ���������� ����� ������. �� ��������� - ������ ����������� ��
%   ����������� � ��������.
%   �������� ����������:
%   Res - �������������(0...1) �������� "y" ����� �������
%% �������� ����������, �������������� ������� P
if nargin==1
    tol=2;
    stepval=50;
    extr=1;
elseif nargin==2
    stepval=50;
    extr=1;
elseif nargin==3
    extr=1;
end
P=P/max(max(P));
[M,N]=size(P);
P=double((1-P)>(1-stepval/100));
% figure(3);
% imagesc(P);
Res=zeros(1,N);
%% ������ �� ��������
for n=1:N
    %% ������� ���������� �� ����� ���� ������� �� �������
    m=1;
    quanzeros=0;%������� �����
    mf=0;%������ ����� ������, ��� ����������� 1
    ml=0;%��������� ����� ������, ��� ����������� 1
    %% ������ �� �������
    while m<=M
        if P(m,n)==1
             if n==1
                if ~mf
                    mf=m;
                end
                ml=m;
            elseif m>=Res(n-1)-tol*N/100 && m<=Res(n-1)+tol*N/100
                %�������� ����� �� ��������� � ���� �������
                if ~mf
                    mf=m;
                end
                ml=m;
             end
        else
            if mf
            %���� ������ �������� �������, �� ��� 5 ����� ��� ����� ������              
                quanzeros=quanzeros+1;   
                if quanzeros==5
                    m=M;                   
                end
            end
        end
        m=m+1;
        %% ������ �������� ������� � ���� �����
        %���� ������ ������� ��� ������, �� ������ ������� ��������������
        %������� � ����������. ���� ��� - ���������� �������������.
        if m>M
            if mf
                Res(n)=(mf+ml)/2;
            elseif n>extr
                s=(1-2*abs(round(extr)-extr));
                Res(n)=s*(3*Res(n-1)-Res(n-round(extr)))/2+(1-s)*Res(n-1);
                if Res(n)<0
                    Res(n)=0;
                end
            end
        end
        
    end
end
%% ����� ������������ ���������� ��������
if nargin~=5
    Res=1-Res/M;
else
    interp=@(x0,y0,x1,y1,x)[y0 y1]/[x0 x1; 1 1]*[x;1];
    NewRes=zeros(1,Q);
    k=0;
    for q=linspace(1,N,Q)
        k=k+1;
        for n=1:N
            if q==n
                NewRes(k)=Res(n);
                break;
            elseif n>q
                NewRes(k)=interp(n-1,Res(n-1),n,Res(n),q);
                break;
            end
        end
    end
    Res=1-NewRes/M;
end
end

