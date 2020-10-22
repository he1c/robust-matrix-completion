function M=HQ_PF(MissM,Mask,option)

    U=option.U;
    V=option.V;
    e_1=option.stop_1;
    e_2=option.stop_2;
    yita=option.yita;
    sigmamin=option.sigmamin;
    maxitr=option.maxitr;
    
    [n1,n2]=size(MissM);
    [~,r]=size(U);
    J_pre=0;
    
    for k=1:1:maxitr
        
        for j=1:1:n2
            Ij=find(Mask(:,j)>0);
            Uj=U(Ij,:);
            bj=MissM(Ij,j);
            if sum(Ij>0)>r
                V(:,j)=(Uj'*Uj)^-1*Uj'*bj;
            end
        end
        
        for i=1:1:n1
            Ji=find(Mask(i,:)>0);
            Vi=V(:,Ji);
            bi=MissM(i,Ji);
            if sum(Ji>0)>r
                U(i,:)=(Vi*Vi')\(Vi*bi');
            end
        end
        
        M=U*V;
        J=norm(Mask.*(MissM-M),'fro');
        
        if k>0&&abs(J-J_pre)<e_1

            for kk=1:1:maxitr

                AAA=reshape(MissM-Mask.*(U*V),n1*n2,1);
                AAA(AAA==0)=[];
                sigma=max((quantile(AAA,0.75)-quantile(AAA,0.25))*yita,sigmamin);

                for j=1:1:n2
                    Ij=find(Mask(:,j)>0);
                    Uj=U(Ij,:);
                    bj=MissM(Ij,j);
                    if sum(Ij>0)>r
                        V(:,j)=HQ_IRLS_solver(Uj,bj,sigma);
                    elseif sum(Ij>0)==r
                        V(:,j)=(Uj'*Uj)^-1*Uj'*bj;
                    end
                end

                for i=1:1:n1
                    Ji=find(Mask(i,:)>0);
                    Vi=V(:,Ji);
                    bi=MissM(i,Ji);
                    if sum(Ji>0)>r
                        U(i,:)=HQ_IRLS_solver(Vi',bi',sigma);
                    elseif sum(Ji>0)==r
                        U(i,:)=(Vi*Vi')^-1*Vi*bi';
                    end
                end

                M=U*V;
                J=norm(Mask.*(MissM-M),'fro');
                if abs(J-J_pre)<e_2
                    break;
                end
                J_pre=J;            
            end
            
            M=U*V;
            
            break;
        end
        
        J_pre=J;
        
        
    end
            
end