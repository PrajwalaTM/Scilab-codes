function y=_f(x)
	y=x(1)^2+x(2)^2+2*x(3)^2+x(4)^2-5*x(1)-5*x(2)-21*x(3)+7*x(4)
endfunction

x0 = [0,0,0,0];
A = [];
b = [];
Aeq = [];
beq = [];
no_nlic=3;
lb=[	];
ub=[];

//exec builder.sce
//exec loader.sce
function [y]=_nlc(x)
	y(1)=x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(1)-x(2)+x(3)-x(4)-8;
	y(2)=x(1)^2+2*x(2)^2+x(3)^2+2*x(4)^2-x(1)-x(4)-10;
	y(3)=2*x(1)^2+x(2)^2+x(3)^2+2*x(1)-x(2)-x(4)-5;
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");


[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc)

