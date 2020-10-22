function [LRM,Mask]=LRMatrix(n1,n2,r,p)

    M1=randn(n1,r);
    M2=randn(n2,r);
    LRM=M1*M2';
    
    ind_l=zeros(n1*n2,1);
    vis=randperm(n1*n2,floor(p*n1*n2));
    ind_l(vis)=1;
    Mask=reshape(ind_l,n1,n2);

end