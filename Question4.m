fixed = imread('frc1.tif');
moving = imread('frc2.tif');
cpselect(moving,fixed);
X1=fixedPoints;
X2=movingPoints;

F=estimateF(X1',X2');
displayEpipolarF(fixed,moving,F);

%%normalization
function [Xn,T]=normal(X)
    Xc=sum(X(:,1))/size(X,1);
    Yc=sum(X(:,2))/size(X,1);
    Sum_d=sqrt(sum((X(:,1)-Xc).*(X(:,1)-Xc))+sum((X(:,2)-Yc).*(X(:,2)-Xc)));
    k=sqrt(2)*size(X,1)/Sum_d;
    T=[k 0 -k*Xc; 0 k -k*Yc; 0 0 1];
    Xn=[k*(X(:,1)-Xc) k*(X(:,2)-Yc)];
end

function F = estimateF(x1, x2)
    x1=x1';
    x2=x2';
    [Xn1,T1]=normal(x1);
    [Xn2,T2]=normal(x2);
    A=[Xn2(:,1).*Xn1(:,1) Xn2(:,1).*Xn1(:,2) Xn2(:,1) Xn2(:,2).*Xn1(:,1) Xn2(:,2).*Xn1(:,2) Xn2(:,2) Xn1(:,1) Xn1(:,2) ones(size(x1,1),1)];
    [U,S,V]=svd(A);
    f=V(:,size(V,2));
    F=reshape(f,[3,3])';
    [u,s,v]=svd(F);
    s(3,3)=0;
    F=u*s*v';
    F=T2'*F*T1;
end



    
    
    

