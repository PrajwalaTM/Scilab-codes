function y=_f(x)
	y=-(x(1)+x(2));
endfunction

x0 = [0.1,0.1];
A = [];
b = [];
Aeq = [];
beq = [];
no_nlic=2;
lb=[0,0];
ub=[];

exec builder.sce
exec loader.sce
function [y]=_nlc(x)
	y(1)=x(2)-x(1)^1.5;
	y(2)=-(x(1)+x(2));
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");


[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc)

