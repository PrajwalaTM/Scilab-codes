//Objective function to be minimised
function y=_f(x)
	y=-(x(1)^2+x(2)^2+x(3)^2);
endfunction

//Starting point, linear constraints and variable bounds
x0=[0.1,0.1,0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[]
ub=[]

//exec builder.sce;
//exec loader.sce;

//Nonlinear constraints
function [c,ceq]=nlc(x)
	c = [x(1) , x(2), x(3)];
    	ceq = [];
endfunction

//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "Hessian", "OFF","GradCon", "OFF");

//Gradient of objective function
function y= funG(x)
	y= [-2*x(1),-2*x(2),-2*x(3)];
endfunction

//Hessian of lagrangian
function y= funH(x,obj,lambda)
	y= obj*[-2,0,0;0,-2,0;0,0,-2] + lambda(1)*[0,0,0;0,0,0;0,0,0] + lambda(2)*[0,0,0;0,0,0;0,0,0] + lambda(3)*[0,0,0;0,0,0;0,0,0]
endfunction

//Gradient of constraints
function [cg,ceqg]= cg(x)
	cg=[1 , 0 , 0 ; 0 , 1 , 0 ; 0 , 0 , 1];
	ceqg=[];
	//or, this is another way of writing it : y=[2*x(1),-2*x(2),2*x(3),2*x(1),2*x(2),2*x(3)];	
endfunction

//Calling fmincon to solve the given problem
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,nlc,options)

