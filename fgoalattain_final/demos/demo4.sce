function f1 = objfun(x)
f1(1)=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
f1(2)=x(2)-x(1)*5+x(2)*x(2)
f1(3)=x(2)-4*x(1)
f1(4)=x(2)^3-x(1)^2+4
endfunction

goal=[5,-6,7,2]
weight=[8,2,3,4]
x0=[-1,2]
A=[1,2]
b=[3]
Aeq=[-1,4]
beq=[5]
lb=[-1,-1]
ub=[10,10]
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b,Aeq,beq,lb,ub)
