function y=_f(x)
y=-5*x(1)-4*x(2)-6*x(3)
endfunction

x0=[0.1,15.1,3.1]
A=[1 -1 1;3 2 4;3 2 0]
b=[20;42;30]
Aeq=[]
beq=[]
lb=zeros(3,1)
ub=[]

[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub)
