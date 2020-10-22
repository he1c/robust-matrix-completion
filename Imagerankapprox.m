function output=Imagerankapprox(I,r)
    
    [S,V,D]=svd(I);
    V(r+1:end,r+1:end)=0;
    output=S*V*D';

end