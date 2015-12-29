function y=_f(x)
	y=-(x(1)+x(2));
endfunction

x0=[0.1,0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[0 0]
ub=[]
no_nlic=2;

function [y]=_nlc(x)
	y(1)=x(1)^2+x(2)^2-2
	y(2)=1-x(1)^2-x(2)^2

endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "HessObj", "ON","GradCon", "ON");
function y= _funG(x)
	y= [-1,-1];
endfunction

function y= _funH(x)
	y= [0, 0; 0, 0];
endfunction


function y= _funJ(x)
	y= [2*x(1),2*x(2),-2*x(1),-2*x(2)];
endfunction



[x,fval,exitflag,output,lambda,grad,hessian,zl] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_funG,_funH,_funJ)
