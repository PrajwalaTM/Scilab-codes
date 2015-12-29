// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: rvidhyadar@gmail.com & vignesh2496@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt



function [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl,zu] = fmincon (varargin)
  // Solves a Constrainted Optimization Problem
  //
  //   Calling Sequence
  //   xopt = fmincon(f,x0,A,b)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,lHess)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,cGrad) 
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,cGrad)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,lHess,cGrad)
  //   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess,cGrad)
  //   [xopt,fopt] = fmincon(.....)
  //   [xopt,fopt,exitflag]= fmincon(.....)
  //   [xopt,fopt,exitflag,output]= fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl,zu]=fmincon(.....)
  //
  //   Parameters
  //   f : a function, representing objective function of the problem 
  //   x0 : a vector of doubles, containing starting values of variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  //   A : a matrix of doubles, containing coefficients of Linear Inequality Constraints of size (m X n) where 'm' is the number of Linear Inequality Constraints
  //   b : a vector of doubles, related to 'A' and containing the Right hand side equation of the Linear Inequality Constraints of size (m X 1)
  //   Aeq : a matrix of doubles, containing coefficients of Linear Equality Constraints of size (m1 X n) where 'm1' is the number of Linear Equality Constraints
  //   beq : a vector of doubles, related to 'Aeq' and containing the Right hand side equation of the Linear Equality Constraints of size (m1 X 1)
  //   lb : a vector of doubles, containing lower bounds of the variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  //   ub : a vector of doubles, containing upper bounds of the variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  //   nlc : a function, representing Non-linear Constraints functions(both Equality and Inequality) of the problem. It is declared in such a way that non-linear Inequality constraints are defined first as a single row vector (c), followed by non-linear Equality constraints as another single row vector (ceq)
  //   options : a list, containing option for user to specify -Maximum iteration, Maximum CPU-time, GradObj, Hessian & GradCon.                                       Syntax for options- options= list("MaxIter", [---], "CpuTime", [---], "GradObj", "ON/OFF", "Hessian", "ON/OFF", "GradCon", "ON/OFF");                                             Default Values for Options==> ("MaxIter", [10000], "CpuTime", [600], "GradObj", "OFF", "Hessian", "OFF", "GradCon", "OFF");
  //   fGrad : a function, representing gradient function of the Objective in Vector Form
  //   lHess : a function, representing hessian function of the Lagrange in Symmetric Matrix Form
  //   cGrad : a function, representing gradient of the Non-Linear Constraints (both Equality and Inequality) of the problem. It is declared in such a way that gradient of non-linear Inequality constraints are defined first as a separate Matrix (cg of size m2 X n as empty), followed by gradient of non-linear Equality constraints as a separate Matrix (ceqg of size m2 X n or as empty) where m2 & m3 are number of Non-linear Inequality and Equality constraints respectively
  //   xopt : a vector of doubles, cointating the computed solution of the optimization problem
  //   fopt : a scalar of double, containing the function value at x
  //   exitflag : a scalar of integer, containing flag which denotes the reason for termination of algorithm
  //   output : a structure, containing information about the optimization
  //   lambda : a vector of doubles, containing Lagrange multipliers at the optimized point
  //   gradient : a vector of doubles, containing Objective's gradient of the optimized point
  //   hessian  : a matrix of doubles, containing Objective's hessian of the optimized point
  //   zl : a vector of doubles, containing lower bound multipliers
  //   zu : a vector of doubles, containing upper bound multipliers
  //
  //   Description
  //   Search the minimum of a Constrained optimization problem specified by :
  //   find the minimum of f(x) such that 
  //
  //   <latex>
  //    \begin{eqnarray}
  //    &\mbox{min}_{x}
  //    & f(x) \\
  //    & \text{subject to} & A*x \leq b \\
  //    & & Aeq*x \ = beq\\
  //	& & c(x) \leq  0\\
  //    & & ceq(x) \ =  0\\
  //    & & lb \leq x \leq ub \\
  //    \end{eqnarray}
  //   </latex>
  //
  //   We are calling IPOpt for solving the Constrained problem, IPOpt is a library written in C++.
  //
  // Examples
  //
  //	//Find x in R^2 such that it minimizes:
  //    //f(x)= -x1 -x2/3
  //    //x0=[0,0]
  //    //constraint-1 (c1): x1 + x2 <= 2
  //    //constraint-2 (c2): x1 + x2/4 <= 1 
  //    //constraint-3 (c3): x1 - x2 <= 2
  //    //constraint-4 (c4): -x1/4 - x2 <= 1
  //    //constraint-5 (c5): -x1 - x2 <= -1
  //    //constraint-6 (c6): -x1 + x2 <= 2
  //    //constraint-7 (c7): x1 + x2 = 2  
  //
  //    //Objective function to be minimised
  //    function y=f(x)
  //		y=-x(1)-x(2)/3;
  //	endfunction
  //
  //	//Starting point, linear constraints and variable bounds  
  //	x0=[0 , 0 , 0];
  //	A=[1,1 ; 1,1/4 ; 1,-1 ; -1/4,-1 ; -1,-1 ; -1,1];
  //	b=[2;1;2;1;-1;2];
  //	Aeq=[1,1];
  // 	beq=[2];
  //	lb=[];
  //	ub=[];
  //
  //	//Options
  //	options=list("GradObj", "ON", "Hessian", "ON","GradCon", "OFF");
  //
  //	//Gradient of objective function
  //	function y= fGrad(x)
  //		y= [-1,-1/3];
  // 	endfunction
  //
  //	//Hessian of lagrangian
  //	function y= lHess(x,obj,lambda)
  //		y= obj*[0,0;0,0] 
  //  	endfunction
  //
  //    //Calling IPopt
  //	[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess)
  //
  // Examples
  //
  //	//Find x in R^3 such that it minimizes:
  //    //f(x)= x1*x2 + x2*x3
  //    //x0=[0.1 , 0.1 , 0.1]
  //    //constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
  //    //constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10  
  //
  //    //Objective function to be minimised
  //    function y=f(x)
  //		y=x(1)*x(2)+x(2)*x(3);
  //	endfunction
  //
  //	//Starting point, linear constraints and variable bounds  
  //	x0=[0.1 , 0.1 , 0.1];
  //	A=[];
  //	b=[];
  //	Aeq=[];
  // 	beq=[];
  //	lb=[];
  //	ub=[];
  //
  //	//Nonlinear constraints  
  //	function [c,ceq]=nlc(x)
  //		c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
  //  		ceq = [];
  //	endfunction
  //
  //	//Options  
  //	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");
  //
  //	//Gradient of objective function
  //	function y= fGrad(x)
  //		y= [x(2),x(1)+x(3),x(2)];
  //	endfunction
  //
  //    //Hessian of the Lagrange Function
  //	function y= lHess(x,obj,lambda)
  //		y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
  //	endfunction
  //
  //    //Gradient of Non-Linear Constraints
  //	function [cg,ceqg] = cGrad(x)
  //		cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
  //		ceqg=[];
  //	endfunction
  //
  //    //Calling IPopt
  //	[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess,cGrad)
  //
  // Examples
  //
  //    The below Problem is an Unbounded problem:
  //	//Find x in R^3 such that it minimizes:
  //    //f(x)= -(x1^2 + x2^2 + x3^2)
  //    //x0=[0.1 , 0.1 , 0.1]
  //    //  x1 <= 0
  //    //  x2 <= 0
  //    //  x3 <= 0
  //
  //    //Objective function to be minimised
  //    function y=f(x)
  //		y=-(x(1)^2+x(2)^2+x(3)^2);
  //	endfunction
  //
  //	//Starting point, linear constraints and variable bounds  
  //	x0=[0.1 , 0.1 , 0.1];
  //	A=[];
  //	b=[];
  //	Aeq=[];
  // 	beq=[];
  //	lb=[];
  //	ub=[0,0,0];
  //
  //	//Options
  //	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF", "Hessian", "OFF","GradCon", "OFF");
  //
  //    //Calling IPopt
  //	[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,[],options)
  //
  // Examples
  //
  //    The below Problem is an Infeasible problem:
  //	//Find x in R^3 such that in minimizes:
  //    //f(x)=x1*x2 + x2*x3
  //    //x0=[1,1,1]
  //    //constraint-1 (c1): x1^2 <= 1
  //    //constraint-2 (c2): x1^2 + x2^2 <= 1    
  //    //constraint-3 (c3): x3^2 <= 1  
  //    //constraint-4 (c4): x1^3 = 0.5  
  //    //constraint-5 (c5): x2^2 + x3^2 = 0.75
  //    // 0 <= x1 <=0.6
  //    // 0.2 <= x2 <= inf
  //    // -inf <= x3 <= 1
  //
  //    //Objective function to be minimised
  //    function y=f(x)
  //		y=x(1)*x(2)+x(2)*x(3);
  //	endfunction
  //
  //	//Starting point, linear constraints and variable bounds  
  //	x0=[1,1,1];
  //	A=[];
  //	b=[];
  //	Aeq=[];
  // 	beq=[];
  //	lb=[0 0.2,-%inf];
  //	ub=[0.6 %inf,1];
  //
  //	//Nonlinear constraints  
  //	function [c,ceq]=nlc(x)
  //		c=[x(1)^2-1,x(1)^2+x(2)^2-1,x(3)^2-1];
  //		ceq=[x(1)^3-0.5,x(2)^2+x(3)^2-0.75];
  //	endfunction
  //
  //	//Options  
  //	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "Hessian", "ON","GradCon", "ON");
  //
  //	//Gradient of objective function
  //	function y= fGrad(x)
  //		y= [x(2),x(1)+x(3),x(2)];
  //	endfunction
  //
  //    //Hessian of the Lagrange Function
  //	function y= lHess(x,obj,lambda)
  //		y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,0,0;0,0,0] + lambda(2)*[2,0,0;0,2,0;0,0,0] +lambda(3)*[0,0,0;0,0,0;0,0,2] + lambda(4)*[6*x(1),0,0;0,0,0;0,0,0] + lambda(5)*[0,0,0;0,2,0;0,0,2];
  //	endfunction
  //
  //    //Gradient of Non-Linear Constraints
  //	function [cg,ceqg] = cGrad(x)
  //		cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
  //		ceqg = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
  //	endfunction
  //
  //    //Calling IPopt
  //	[x,fval,exitflag,output,lambda,grad,hessian,zl,zu] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options,fGrad,lHess,cGrad)
  //
  // Authors
  // R.Vidyadhar , Vignesh Kannan
 

	//To check the number of input and output argument
   	[lhs , rhs] = argn();
	
	//To check the number of argument given by user
   	if ( rhs<4 | rhs>13 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,9,10,11,12,13"), "fmincon", rhs);
    		error(errmsg)
   	end
    	
	if (rhs==5 | rhs==7) then
    	errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,9,10,11,12,13"), "fmincon", rhs);
    	error(errmsg)
   	end
 
	//Storing the Input Parameters  
   	f    	 = varargin(1);
   	x0   	 = varargin(2);
   	A    	 = varargin(3);
   	b    	 = varargin(4);
   	Aeq  	 = [];
   	beq  	 = [];
   	lb       = [];
   	ub       = [];
   	nlc     = [];
   	
   	if (rhs>4) then
   		Aeq  	 = varargin(5);
   		beq  	 = varargin(6);
   	end

   	if (rhs>6) then
   		lb       = varargin(7);
   		ub       = varargin(8);
   	end

   	if (rhs>8) then
   		nlc      = varargin(9);
	end
	 
	//To check whether the 1st Input argument (f) is a function or not
   	if (type(f) ~= 13 & type(f) ~= 11) then
   		errmsg = msprintf(gettext("%s: Expected function for Objective (1st Parameter) "), "fmincon");
   		error(errmsg);
   	end
   
	//To check whether the 2nd Input argument (x0) is a Vector/Scalar
   	if (type(x0) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Starting Point (2nd Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check and convert the 2nd Input argument (x0) to row Vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector or Column Vector for x0 (Starting Point) or Starting Point cannot be Empty"), "fmincon");
   		error(errmsg);
    end

   	if(size(x0,2)==1) then
   		x0=x0';		//Converting x0 to row vector, if it is column vector
   	else 
   	 	x0=x0;		//Retaining the same, if it is already row vector
   	end   	 	
    	s=size(x0);
  	
  	//To check the match between f (1st Parameter) & x0 (2nd Parameter)
   	if(execstr('init=f(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "fmincon");
   		error(errmsg);
	end
   	
  	//To check whether the 3rd Input argument (A) is a Matrix/Vector
   	if (type(A) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Matrix/Vector for Constraint Matrix A (3rd parameter)"), "fmincon");
   		error(errmsg);
  	end

	//To check for correct size of A(3rd paramter)
   	if(size(A,2)~=s(2) & size(A,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of Linear Inequality Constraints X No of Variables) or an Empty Matrix for Linear Inequality Constraint coefficient Matrix A"), "fmincon");
   		error(errmsg);
   	end

   	s1=size(A);
   	
	//To check whether the 4th Input argument (b) is a Vector/Scalar
   	if (type(b) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for b (4th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of b (4th paramter) and Converting into Column Vector which is required for Ipopt
    if(s1(2)==0) then
    	if(size(b,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Inequality Constraint coefficient Matrix A (3rd parameter) is empty, b (4th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
	else
   		if((size(b,1)~=1) & (size(b,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for b (4th Parameter) for your Inputs "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)~=s1(1) & size(b,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (number of Linear inequality Constraints X 1) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==s1(1) & size(b,2)==1) then 
   	 		b=b;
   		elseif(size(b,1)==1 & size(b,2)~=s1(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Linear inequality Constraints) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==1 & size(b,2)==s1(1)) then
   			b=b';
   		end 
   	end
  	
  	//To check whether the 5th Input argument (Aeq) is a Matrix/Vector
   	if (type(Aeq) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Matrix/Vector for Equality Constraint Matrix Aeq (5th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of Aeq (5th paramter)
   	if(size(Aeq,2)~=s(2) & size(Aeq,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of Linear Equality Constraints X No of Variables) or an Empty Matrix for Linear Equality Constraint coefficient Matrix Aeq"), "fmincon");
   		error(errmsg);
   	end

   	s2=size(Aeq);

	//To check whether the 6th Input argument(beq) is a Vector/Scalar
   	if (type(beq) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for beq (6th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of beq(6th paramter) and Converting into Column Vector which is required for Ipopt
    if(s2(2)==0) then
    	if(size(beq,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Equality Constraint coefficient Matrix Aeq (5th parameter) is empty, beq (6th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
   	else
   		if((size(beq,1)~=1) & (size(beq,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for beq (6th Parameter)"), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)~=s2(1) & size(beq,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (number of Linear Equality Constraints X 1) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==s2(1) & size(beq,2)==1) then 
   	 		beq=beq;
   		elseif(size(beq,1)==1 & size(beq,2)~=s2(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Linear Equality Constraints) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==1 & size(beq,2)==s2(1)) then
   			beq=beq';
   		end 
   	end
   	
  	
  	//To check whether the 7th Input argument (lb) is a Vector/Scalar
   	if (type(lb) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Lower Bound Vector (7th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size and data of lb (7th paramter) and Converting it to Column Vector as required by Ipopt
   	if (size(lb,2)==0) then
        	lb = repmat(-%inf,1,s(2));
    end
    
   	if (size(lb,1)~=1) & (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound (7th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(lb,1)~=s(2) & size(lb,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (number of Variables X 1) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==s(2) & size(lb,2)==1) then
   	 	lb=lb;
   	elseif(size(lb,1)==1 & size(lb,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Variables) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==1 & size(lb,2)==s(2)) then
   		lb=lb';
   	end 
   	
   	//To check whether the 8th Input argument (ub) is a Vector/Scalar
   	if (type(ub) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Upper Bound Vector (8th Parameter)"), "fmincon");
   		error(errmsg);
  	end
   	
   	//To check for correct size and data of ub (8th paramter) and Converting it to Column Vector as required by Ipopt
    if (size(ub,2)==0) then
        ub = repmat(%inf,1,s(2));
    end
    
    if (size(ub,1)~=1)& (size(ub,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound (8th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(ub,1)~=s(2) & size(ub,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (number of Variables X 1) for upper bound (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==s(2) & size(ub,2)==1) then
   	 	ub=ub;
   	elseif(size(ub,1)==1 & size(ub,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Variables) for beq (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==1 & size(ub,2)==s(2)) then
   		ub=ub';
   	end 
    
    //To check the contents of lb & ub (7th & 8th Parameter)
    for i = 1:s(2)
		if (lb(i) == %inf) then
		   	errmsg = msprintf(gettext("%s: Value of Lower Bound can not be infinity"), "fmincon");
    		error(errmsg); 
  		end	

		if (ub(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of Upper Bound can not be negative infinity"), "fmincon");
    		error(errmsg); 
		end	

		if(ub(i)-lb(i)<=1e-6) then
			errmsg = msprintf(gettext("%s: Difference between Upper Bound and Lower bound should be atleast > 10^6 for variable number= %d "), "fmincon", i);
    		error(errmsg)
    	end
	end
  	
	//To check whether the 10th Input argument (nlc) is a function or an empty Matrix
	if (type(nlc) == 1 & size(nlc,2)==0 ) then
  		addnlc=[];
  		no_nlc=0;
  		no_nlic=0;
		no_nlec=0;
  	elseif (type(nlc) == 13 | type(nlc) == 11) then
  		
  		if(execstr('[sample_c,sample_ceq] = nlc(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Non-Linear Constraint function(9th Parameter) and x0(2nd Parameter) did not match"), "fmincon");
   			error(errmsg);
		end
  		[sample_c,sample_ceq] = nlc(x0);
  		
  		if (size(sample_c,1)~=1 & size(sample_c,1)~=0) then
  			errmsg = msprintf(gettext("%s: Definition of c in Non-Linear Constraint function(9th Parameter) should be in the form of Row Vector or Empty Vector"), "fmincon");
    		error(errmsg)
    	end
  		
  		if (size(sample_ceq,1)~=1 & size(sample_ceq,1)~=0) then
  			errmsg = msprintf(gettext("%s: Definition of ceq in Non-Linear Constraint function(9th Parameter) should be in the form of Row Vector or Empty Vector"), "fmincon");
    		error(errmsg)
    	end
    	
  		no_nlic = size(sample_c,2);
  		no_nlec = size(sample_ceq,2);
  		no_nlc = no_nlic + no_nlec;
  		
  		function allcon = addnlc(x)
  			[c,ceq] = nlc(x);
  			allcon = [c';ceq'];
  		endfunction
  		
  		if(execstr('sample_allcon = addnlc(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Non-Linear Constraint function(9th Parameter) and x0(2nd Parameter) did not match"), "fmincon");
   			error(errmsg);
		end
  		sample_allcon = addnlc(x0);
  		
  		if (size(sample_allcon,1)~=no_nlc | size(sample_allcon,2)~=1) then
  			errmsg = msprintf(gettext("%s: Please check the Non-Linear Constraint function (9th Parameter) function"), "fmincon");
    		error(errmsg)
    	end
    else 
    	errmsg = msprintf(gettext("%s: Non Linear Constraint (9th Parameter) should be a function or an Empty Matrix"), "fmincon");
    		error(errmsg)
    end
  		
  
   	//To check, Whether Options is been entered by user   
   	if ( rhs<10 ) then
      		param = list();
       else
      		param =varargin(10); //Storing the 3rd Input Parameter in intermediate list named 'param'
    end
   
	//If Options is entered then checking its type for 'list'   
   	if (type(param) ~= 15) then
   		errmsg = msprintf(gettext("%s: Options (10th parameter) should be a list"), "fmincon");
   		error(errmsg);
   	end
   
	//If Options is entered then checking whether even number of entires are entered   
   	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of Options (list) should be even"), "fmincon");
		error(errmsg);
   	end

	//To set Default Value for Options, If User Doesn't enter Options
   	options = list(..
      		"MaxIter"     , [10000], ...
      		"CpuTime"   , [600] ...
		);

	//Flags to check whether Gradient is "ON"/"OFF" & Hessian is "ON"/"OFF" 
   	flag1=0;
   	flag2=0;
   	flag3=0;
   	fGrad=[];
   	lHess=[];
   	cGrad=[];
	addcGrad=[];
   	
 
	//To check the User Entry for Options and storing it
   	for i = 1:(size(param))/2
       	select param(2*i-1)
    		case "MaxIter" then
          			options(2*i) = param(2*i);    //Setting the Maximum Iteration as per user entry
       		case "CpuTime" then
          			options(2*i) = param(2*i);    //Setting the Maximum CPU Time as per user entry
        	case "GradObj" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Gradient function if Gradient Option is "ON"
        					if (rhs<11) then      
				     			errmsg = msprintf(gettext("%s: Gradient function of Objective is missing, but GradObj=ON"), "fmincon");
				    			error(errmsg);     			
        					else
        						//This flag1 is activated(ie. =1) if Gradient is supplied
        						flag1=1;
        						pos_fGrad=11;
        						fGrad=varargin(11);
        					end
        				//To check whether Wrong entry(other than ON/OFF) is entered
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for GradObj should be either ON or OFF"), "fmincon");
						error(errmsg);     	
        				end
        	case "Hessian" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Hessian function if Hessian Option is "ON"
							if (flag1==0) then
								if (rhs<11) then    
				     				errmsg = msprintf(gettext("%s: Hessian function of Objective is missing, but HessObj=ON"), "fmincon");
				     				error(errmsg);		
        						else
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag2=1;
        							pos_lHess=11;
        							lHess=varargin(11);
        			    		end         			
        					elseif (flag1==1) then
								if (rhs<12) then    
				     				errmsg = msprintf(gettext("%s: Hessian function of Objective is missing, but HessObj=ON"), "fmincon");
				     				error(errmsg);     			
        						else
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag2=1;
        							pos_lHess=12;
        							lHess=varargin(12);
        			    		end
        			    	end       
        				//To check whether Wrong entry(other than ON/OFF) is entered	            
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for HessObj should be either ON or OFF"), "fmincon");
							error(errmsg);   
        				end
        	case "GradCon" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Gradient function if Gradient Option is "ON"
        					if (flag1==0 & flag2==0) then
        						if (rhs<11) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Non-Linear Constraint is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);     			
        						else
        							pos_cGrad=11;
        							//This flag is activated(ie. =1) if Gradient is supplied
        							flag3=1;
        							cGrad=varargin(11);
        						end
        					elseif((flag1==1 & flag2==0) |(flag1==0 & flag2==1) ) then
        						if (rhs<12) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Constraints is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);    			
        						else
        							pos_cGrad=12;
        							//This flag is activated(ie. =1) if Gradient is supplied
        							flag3=1;
        							cGrad=varargin(12);
        						end
        					elseif(flag1==1 & flag2==1) then
        						if (rhs<13) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Constraints is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);      			
        						else
        							pos_cGrad=13;
        							//This flag is activated(ie. =1) if Gradient is supplied
        							flag3=1;
        							cGrad=varargin(13);
        						end        					
        					end 	       				      
        			//To check whether Wrong entry(other than ON/OFF) is entered
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for GradCon should be either ON or OFF, but GradCon=ON"), "fmincon");
							error(errmsg);     	
        				end
            else
    	      	errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "fmincon", param(2*i-1));
    	      	error(errmsg)
    		end
   	end
   
   
	//Defining a function to calculate Gradient or Hessian if the respective user entry is OFF 
   	function y=gradhess(x,t)
		if t==1 then	//To return Gradient
			y=numderivative(f,x)		
		elseif t==2 then		//To return Hessian]n
			[grad,y]=numderivative(f,x)
		elseif t==3 then	//To return Gradient
			y=numderivative(addnlc,x)		
		elseif t==4 then		//To return Hessian]n
			[grad,y]=numderivative(addnlc,x)
		end
   	endfunction
 
   //To check the correct number of inputs given by the user	
   if (flag1==0 & flag2==0 & flag3==0)
   		if(rhs>10) then
        	errmsg = msprintf(gettext("%s: Only 10 Inputs are Needed for this option(GradObj=OFF, HessObj=OFF, GradCon=OFF), but %d were recorded"), "fmincon",rhs);
			error(errmsg); 
		end
   elseif ((flag1==1 & flag2==0 & flag3==0) | (flag1==0 & flag2==1 & flag3==0) | (flag1==0 & flag2==0 & flag3==1)) then
  		if(rhs>11) then
        	errmsg = msprintf(gettext("%s: Only 11 Inputs were needed for this option, but %d were recorded"), "fmincon",rhs);
			error(errmsg);
		end
   elseif ((flag1==1 & flag2==1 & flag3==0) | (flag1==0 & flag2==1 & flag3==1) | (flag1==1 & flag2==0 & flag3==1)) then
   		if(rhs>12) then
        	errmsg = msprintf(gettext("%s: Only 12 Inputs were needed for this option, but %d were recorded"), "fmincon",rhs);
			error(errmsg);
		end
   elseif (flag1==1 & flag2==1 & flag3==1)
   		if(rhs>13) then
        	errmsg = msprintf(gettext("%s: Only 13 Inputs are Needed for this option(GradObj=ON, HessObj=ON, GradCon=ON), but %d were recorded"), "fmincon",rhs);
			error(errmsg); 
		end
	end
	
   //To check the correct input of Gradient and Hessian Functions from Users	     	
   if (flag1==1) then
   		if (type(fGrad) ~= 11 & type(fGrad) ~= 13) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective, since GradObj=ON"), "fmincon");
   			error(errmsg);
   		end
   		if(execstr('sample_fGrad=fGrad(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		sample_fGrad=fGrad(x0);
		if (size(sample_fGrad,1)==s(2) & size(sample_fGrad,2)==1) then
		elseif (size(sample_fGrad,1)==1 & size(sample_fGrad,2)==s(2)) then
		elseif (size(sample_fGrad,1)~=1 & size(sample_fGrad,2)~=1) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Gradient function(%dth Parameter)---->Vector function is Expected"), "fmincon",pos_fGrad);
   			error(errmsg);
   		end
   end
   if (flag2==1) then
   		if (type(lHess) ~= 11 & type(lHess) ~= 13) then
  			errmsg = msprintf(gettext("%s: Expected function for Hessian of Objective, since Hessian=ON"), "fmincon");
   			error(errmsg);
   		end
   		if(execstr('sample_lHess=lHess(x0,1,1:no_nlc)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Hessian function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		sample_lHess=lHess(x0,1,1:no_nlc);
   		if(size(sample_lHess,1)~=s(2) | size(sample_lHess,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Hessian function(%dth Parameter)---->Symmetric Matrix function is Expected "), "fmincon",pos_lHess);
   			error(errmsg);
   		end
   	end
   	if (flag3==1) then
   		if (type(cGrad) ~= 11 & type(cGrad) ~= 13) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Constraint function,since GradCon=ON"), "fmincon");
   			error(errmsg);
   		end
   		
   		if(execstr('[sample_cGrad,sample_ceqg]=cGrad(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Constraint and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		[sample_cGrad,sample_ceqg]=cGrad(x0);
		
		if  (size(sample_cGrad,2)==0) then 		
		elseif (size(sample_cGrad,1)~=no_nlic | size(sample_cGrad,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Definition of (cGrad) in Non-Linear Constraint function(%dth Parameter) should be in the form of (m X n) or Empty Matrix where m is number of Non- linear Inequality Constraints and n is number of Variables"), "fmincon",pos_cGrad);
   			error(errmsg);		
		end
		
		if (size(sample_ceqg,2)==0) then		
		elseif (size(sample_ceqg,1)~=no_nlec | size(sample_ceqg,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Definition of (ceqg) in Non-Linear Constraint function(%dth Parameter) should be in the form of (m X n) or Empty Matrix where m is number of Non- linear Equality Constraints and n is number of Variables"), "fmincon",pos_cGrad);
   			error(errmsg);
		end

		function allcongrad = addcGrad(x)
  			[sample_cGrad,sample_ceqg] = cGrad(x);
  			allcongrad=[sample_cGrad;sample_ceqg];
  		endfunction
   	
   		sample_addcGrad=addcGrad(x0);
		if(size(sample_addcGrad,1)~=no_nlc | size(sample_addcGrad,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s:  Wrong Input for Constraint Gradient function(%dth Parameter) (Refer Help)"), "fmincon",pos_cGrad);
   			error(errmsg);
   		end
   	end
   	
   	empty=0;
   	

   	//Calling the Ipopt Function for solving the above Problem
    [xopt,fopt,status,iter,cpu,obj_eval,dual,lambda,zl,zu,gradient,hessian1] = solveminconp 					(f,gradhess,A,b,Aeq,beq,lb,ub,no_nlc,no_nlic,addnlc,flag1,fGrad,flag2,lHess,flag3,addcGrad,x0,options,empty)	
   
	//Calculating the values for output   	
   	xopt = xopt';
    exitflag = status;
    output = struct("Iterations", [],"Cpu_Time",[],"Objective_Evaluation",[],"Dual_Infeasibility",[]);
   	output.Iterations = iter;
    output.Cpu_Time = cpu;
    output.Objective_Evaluation = obj_eval;
    output.Dual_Infeasibility = dual;

    //Converting hessian of order (1 x (numberOfVariables)^2) received from Ipopt to order (numberOfVariables x numberOfVariables)
    s=size(gradient)
    for i =1:s(2)
    	for j =1:s(2)
			hessian(i,j)= hessian1(j+((i-1)*s(2)))
		end
    end
    
    //In the cases of the problem not being solved return NULL to the output matrices
    if( status~=0 & status~=1 & status~=2 & status~=4 & status~=7 ) then
		xopt=[]
		fopt=[]
		output = struct("Iterations", [],"Cpu_Time",[]);
		output.Iterations = iter;
    	output.Cpu_Time = cpu;
		lambda=[]
		gradient=[]
		hessian=[]
		zl=[]
		zu=[]
    end
		

	//To print Output Message
    select status
    
    	case 0 then
        	printf("\nOptimal Solution Found.\n");
    	case 1 then
        	printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
    	case 2 then
       		printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
    	case 3 then
        	printf("\nStop at Tiny Step\n");
    	case 4 then
        	printf("\nSolved To Acceptable Level\n");
    	case 5 then
        	printf("\nConverged to a point of local infeasibility.\n");
    	case 6 then
        	printf("\nStopping optimization at current point as requested by user.\n");
    	case 7 then
        	printf("\nFeasible point for square problem found.\n");
    	case 8 then 
        	printf("\nIterates diverging; problem might be unbounded.\n");
    	case 9 then
        	printf("\nRestoration Failed!\n");
    	case 10 then
        	printf("\nError in step computation (regularization becomes too large?)!\n");
    	case 11 then
        	printf("\nProblem has too few degrees of freedom.\n");
    	case 12 then
        	printf("\nInvalid option thrown back by IPOpt\n");
    	case 13 then
        	printf("\nNot enough memory.\n");
    	case 15 then
        	printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
    	else
        	printf("\nInvalid status returned. Notify the Toolbox authors\n");
        	break;
        end
    
    //Remark for the user, If the gradient and hessian is send by the User
    if (no_nlc~=0) then
		disp("||||||Please Make sure you have entered Correct Non-linear Constraints Functions (9th Parameter) in proper order -->Scilab Will Calculate Based on your input only||||||");	
    end
    
    //Remark for the user, If the gradient and hessian is send by the User
    if (flag1==1 |flag2==1 |flag3==1) then
		disp("||||||Please Make sure you have entered Correct Functions for Gradient or Hessian as per Help -->Scilab Will Calculate Based on your input only||||||");
    end
    		
endfunction
