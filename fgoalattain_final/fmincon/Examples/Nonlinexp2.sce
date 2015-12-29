function y=_f(x)
	y=x(1)^2+2*x(2);
endfunction

x0 = [0,0];
A = [];
b = [];
Aeq = [];
beq = [];
no_nlic=0;
lb=[];
ub=[];

exec builder.sce
exec loader.sce
function [y]=_nlc(x)
	y(1)=x(1)^2+x(2)^2-1
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "ON");

function [y]=_cg(x)
	y(1)=2*x(1);
	y(2)=2*x(2);
endfunction

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_cg)

