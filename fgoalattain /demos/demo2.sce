function f1 = objfun(x)
f1(1)=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
f1(2)=x(2)-x(1)*5+x(2)*x(2)
endfunction

goal=[5,-6]
weight=[8,2]
x0=[-1,2]
A=[1 2]
b=[3]
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b)
