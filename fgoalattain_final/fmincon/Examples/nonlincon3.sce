function y=_f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0=[0,0];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[]
ub=[]
no_nlic=1;

function [y]=_nlc(x)
	y(1)=x(1)^2 + x(2)^2 - 1;

endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");

[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc)
