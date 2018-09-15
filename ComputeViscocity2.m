function [viscocity]=ComputeViscocity2(x,v,rho,i,j,h,E)
   
    r_ab(1,1)=x(1,1,i)-x(1,1,j);
    r_ab(1,2)=x(1,2,i)-x(1,2,j);
    v_ab(1,1)=v(1,1,i)-v(1,1,j);
    v_ab(1,2)=v(1,2,i)-v(1,2,j);
    
    rho_a=rho(i);
    rho_b=rho(j);
    
	if dot(v_ab,r_ab)>=0
		viscocity=0;
    else
        alpha=1;   
        
        ro_ab=(rho_a+rho_b)/2;
        cs_a=sqrt(E/rho_a); %скорость продольной волны
        cs_b=sqrt(E/rho_b);
        r_ab_norm(1,1)=r_ab(1,1)/norm(r_ab);
        r_ab_norm(1,2)=r_ab(1,2)/norm(r_ab);
        w=dot(v_ab,r_ab_norm);
        v_sig=cs_a+cs_b-3*w;
        viscocity=(-alpha/2)*v_sig*w/ ro_ab;
    end
    