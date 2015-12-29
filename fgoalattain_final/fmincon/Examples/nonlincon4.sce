//Objective function to be minimised
function y=f(x)
	y=x(1)*x(2)+x(2)*x(3);
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
	c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
    	ceq = [];
endfunction

//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");

//Gradient of objective function
function y= funG(x)
	y= [x(2),x(1)+x(3),x(2)];
endfunction

//Hessian of lagrangian
function y= funH(x,obj,lambda)
	y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
endfunction

//Gradient of constraints
function [cg,ceqg]= cg(x)
	cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
	ceqg=[];
	//or, this is another way of writing it : y=[2*x(1),-2*x(2),2*x(3),2*x(1),2*x(2),2*x(3)];	
endfunction

//Calling fmincon to solve the given problem
[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,funG,funH,cg)

