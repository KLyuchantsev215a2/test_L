function [init] =initialization_v(x,v_0,N,sqn)

v=zeros(1,2,N);
    for i=1:(sqn*sqn)
        v(1,1,i)= fix((i-1)/sqn);
        v(1,2,i)=0;
    end
    

init=v;