function y=f(x)
	y=x(1)*x(2)+x(2)*x(3);
endfunction

x0=[1,1,1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0.2,-%inf];
ub=[0.5 0.8,%inf];

function [c,ceq]=nlc(x)
	c=[x(1)^2,x(1)^2+x(2)^2,x(3)^2];
	ceq=[x(1)^3,x(2)^2+x(3)^2];
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");
function y= funG(x)
	y= [x(2),x(1)+x(3),x(2)];
endfunction

function y= funH(x,obj,lambda)
	y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,0,0;0,0,0] + lambda(2)*[2,0,0;0,2,0;0,0,0] +lambda(3)*[0,0,0;0,0,0;0,0,2] + lambda(4)*[6*x(1),0,0;0,0,0;0,0,0] + lambda(5)*[0,0,0;0,2,0;0,0,2];
endfunction

function [cg,ceqg] = conG(x)
	cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
	ceqg = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
endfunction

[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,funG,funH,conG)
