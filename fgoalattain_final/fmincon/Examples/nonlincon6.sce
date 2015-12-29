//Objective function
function y= _f(x)
	y= x(1)^2 + x(2)^2;
endfunction

//Starting point, linear constraints and variable bounds
x0=[1.1,1.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
no_nlic=1;//Number of nonlinear inequality constraints

//Nonlinear constraints
function y= _nlc(x)
	y(1)=x(1)^2 + x(2)^2 -2;//Inequality constraints mentioned before all the equality constraints	     
	y(2)= x(1)^2 + x(2)^2 -1;//Equality constraints
endfunction

//Options
options=list("MaxIter", [1000000], "CpuTime", [60], "GradObj", "ON", "HessObj", "OFF", "GradCon", "ON");

//Gradient of objective function
function y=_fg(x)
	y=[2*x(1),2*x(2)];
endfunction 

//Gradient of constraints
function [y]=_cg(x)
	y=[2*x(1),2*x(2),2*x(1),2*x(2)];	 

endfunction
      
//Calling fmincon to solve the given problem
[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg,_cg)
	
