function y=_f(x)
	y=100*(x(2) - x(1)*x(1))*(x(2) - x(1)*x(1)) + (1 - x(1))*(1 - x(1));
endfunction

x0 = [0.5,0];
A = [1,2];
b = 1;
Aeq = [2,1];
beq = 1;
no_nlic=[];
lb=[];
ub=[];

//exec builder.sce
//exec loader.sce
function [y]=_nlc(x)
	y(1)=(x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2
	y(2)=-x(1)^2
	y(3)=x(2)
endfunction

options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "HessObj", "OFF","GradCon", "OFF");
function y= _funG(x)
	y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction

function y= _funH(x)
	y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction

function [y]= _conG(x)
	y(1)=2*(x(1)-1/3)
	y(2)=2*(x(1)-1/3)
	y(3)=-2*x(1)
	y(4)=x(2)
	y(5)=0
	y(6)=1
endfunction
x = fmincon(_f,x0,A,b,Aeq,beq)

