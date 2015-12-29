function y=_f(x)
	y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
endfunction

x0 = [0.5,0];
A = [1,2];
b = 1;
Aeq = [2,1];
beq = 1;
no_nlic=[];
lb=[];
ub=[];
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");

function y= _funG(x)
	y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction

function y= _funH(x)
	y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction

[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,A,b,Aeq,beq)
