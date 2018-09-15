function [viscocity]=ComputeViscocity(x,v,rho,i,j,h,E)
    
  
	
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
        betta=2;
        h_ab=h;
        
        ro_ab=(rho_a+rho_b)/2;
        cs_a=sqrt(E/rho_a); %скорость продольной волны
        cs_b=sqrt(E/rho_b);
        c_ab=(cs_a+ cs_b)/2;
        
        nu=h_ab*0.1;
        mu=(h*dot(v_ab,r_ab)) / (dot(r_ab,r_ab) + nu^2);
        
        viscocity=(-alpha*c_ab*mu + betta*mu*mu) / ro_ab;
    end
    