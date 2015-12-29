function y=_f(x)
	y=.01 * x(1) * x(1) + x(2) * x(2) - 100;
endfunction

x0 = [-1,-1];
A = [-10,1];
b = [-10];
Aeq = [];
beq = [];
lb=[2,-50];
ub=[50,50];

//exec builder.sce
//exec loader.sce

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");


[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(_f,x0,A,b,Aeq,beq,lb,ub)

