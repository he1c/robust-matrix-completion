function a=HQ_IRLS_solver(B,b,sigma)

[M,N]=size(B);

a=zeros(N,1);

for i=1:1:100
    if i==1
        kernelwidth=1000000;
    else
        kernelwidth=sqrt(1/2/M*(norm(b-B*a,'fro'))^2);
    end
    W=diag(exp(-(b-B*a).^2/kernelwidth^2/2));
    a_pre=a;
    p=B'*W*B;
    q=B'*W*b;
    a=p\q;
    if norm(a-a_pre)<1e-1
        kernelwidth=sigma;
        for ii=1:1:100
            W=diag(exp(-(b-B*a).^2/kernelwidth^2/2));
            a_pre=a;
            p=B'*W*B;
            q=B'*W*b;
            a=p\q;
            if norm(a-a_pre)<1e-5
                break;
            end
        end
        break;
    end
    
end
    
end