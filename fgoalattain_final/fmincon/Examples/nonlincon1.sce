function y=_f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0=[1/4,1/4];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0.2]
ub=[0.5 0.8]
no_nlic=1;

function [y]=_nlc(x)
	y(1)=(x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2

endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "HessObj", "ON","GradCon", "OFF");
function y= _funG(x)
	y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction

function y= _funH(x)
	y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction

[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc)
