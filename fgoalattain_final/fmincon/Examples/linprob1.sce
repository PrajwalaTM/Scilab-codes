function y=f(x)
	y=-x(1)-x(2)/3;
endfunction

x0=[0,0];
A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];
b=[2 1 2 1 -1 2];
Aeq=[1 1];
beq=[2];
lb=[]
ub=[]
nlc=[]

//Options
options=list("GradObj", "ON", "Hessian", "ON","GradCon", "OFF");

//Gradient of objective function
function y= funG(x)
	y= [-1,-1/3];
endfunction

//Hessian of lagrangian
function y= funH(x,obj,lambda)
	y= obj*[0,0;0,0] 
endfunction

//Calling fmincon to solve the given problem
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,nlc,options,funG,funH)
