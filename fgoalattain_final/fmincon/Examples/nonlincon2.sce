function y=_f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0=[0,0];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[-2 -2]
ub=[2 2]
no_nlic=1;

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");

[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub)
