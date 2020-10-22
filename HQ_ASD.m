function M=HQ_ASD(MissM,Mask,option)

    U=option.U;
    V=option.V;
    e_1=option.stop_1;
    e_2=option.stop_2;
    yita=option.yita;
    sigmamin=option.sigmamin;
    maxitr=option.maxitr;

    [n1,n2]=size(MissM);
    changeflag=0;
    J_pre=0;
    Maskv=Mask(:);
    
    for i=1:1:maxitr
       
        AAA=reshape(MissM-Mask.*(U*V),n1*n2,1);
        
        AAA(Maskv==0)=[];
        
        kernelwidth=10000;
        
        if changeflag==1
            kernelwidth=max((quantile(AAA,0.75)-quantile(AAA,0.25))*yita,sigmamin);
        end
        
        GU=Mask.*exp(-(MissM-Mask.*(U*V)).^2/2/kernelwidth^2);
        dU=-GU.*(MissM-Mask.*(U*V))*V';
        du=-dU*((V*V')^(-1));
        tu=-trace(dU'*du)/(norm(sqrt(GU).*Mask.*(du*V),'fro'))^2;
        U=U+tu*du;
        
        dV=-U'*(GU.*(MissM-Mask.*(U*V)));
        dv=(U'*U)^(-1)*dV;
        tv=-trace(dV'*dv)/(norm(sqrt(GU).*Mask.*(U*dv),'fro'))^2;
        V=V+tv*dv;
        
        M=U*V;
        J=(norm(Mask.*(MissM-M),'fro'));
        
        if changeflag==0&&abs(J-J_pre)<e_1
            changeflag=1;
        end
        if i>10&&abs(J-J_pre)<e_2
            break;
        end

        J_pre=J;
    
    end

end


