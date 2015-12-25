// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Authors: Prajwala TM,Sheetal Shalini 
// Organization: FOSSEE, IIT Bombay
// Email: prajwala.tm@gmail.com,sheetalsh456@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
function [x,fval,attainfactor,exitflag,output,lambda] = fgoalattain(varargin)
    // Solves a multiobjective goal attainment problem
    // Calling sequence
    //x = fgoalattain(fun,startpoint,goal,weight)
    //x = fgoalattain(fun,startpoint,goal,weight,A,b)
    //x = fgoalattain(fun,startpoint,goal,weight,A,b,Aeq,beq)
    //x = fgoalattain(fun,startpoint,goal,weight,A,b,Aeq,beq,lb,ub)
    //x = fgoalattain(fun,startpoint,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon)
    //x = fgoalattain(fun,startpoint,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon,options)
    //[x,fval] = fgoalattain(...)
    //[x,fval,attainfactor] = fgoalattain(...)
    //[x,fval,attainfactor,exitflag] = fgoalattain(...)
    //[x,fval,attainfactor,exitflag,output] = fgoalattain(...)
    //[x,fval,attainfactor,exitflag,output,lambda] = fgoalattain(...)
    //
    // Input Parameters
    //  fun: a function that accepts a vector x and returns a vector F
    //  startpoint: a nx1 or 1xn matrix of double, where n is the number of variables.
    //      The initial guess for the optimization algorithm
    //  A: a nil x n matrix of double, where n is the number of variables and
    //      nil is the number of linear inequalities. 
    //      
    //  b: a nil x 1 matrix of double, where nil is the number of linear
    //      inequalities
    //  Aeq: a nel x n matrix of double, where n is the number of variables
    //      and nel is the number of linear equalities. 
    //  beq: a nel x 1 matrix of double, where nel is the number of linear
    //      equalities
    //  lb: a nx1 or 1xn matrix of double, where n is the number of variables. 
    //      The lower bound for x. If lb==[], then the lower bound is 
    //      automatically set to -inf
    //  ub: a nx1 or 1xn matrix of double, where n is the number of variables. 
    //      The upper bound for x. If ub==[], then the upper bound is 
    //      automatically set to +inf
    //  nonlinfun: a function, the nonlinear constraints
    //
    //  Output Parameters
    //  x: a nx1 matrix of double, the computed solution of the optimization problem
    //  fval: a vector of double, the value of functions at x
    //  attainfactor: The amount of over- or underachievement of the goals,γ at the solution.
    //  exitflag: a 1x1 matrix of floating point integers, the exit status
    //  output: a struct, the details of the optimization process
    //  lambda: a struct, the Lagrange multipliers at optimum
    
  // Description
  // fgoalattain solves the goal attainment problem, which is one formulation for minimizing a multiobjective optimization problem.
  // Finds the minimum of a problem specified by:
  // Minimise Y such that
  //
  //<latex>
  //\begin{eqnarray}
  //\mbox{min}_{x,\gamma}  & f(x)-weight \ast \gamma \leq goal \\
  //\mbox{subject to} & c(x) \leq 0 \\
  //                  & c_{eq}(x) = 0 \\
  //                  & Ax \leq b \\
  //                  & A_{eq} x = b_{eq} \\
  //                  & lb \leq x \leq ub
  //\end{eqnarray}
  //</latex>
  //
  // The solver makes use of fmincon to find the minimum.
  // The objective function must have header :
  // <programlisting>
  //   f = objfun ( x )
  // </programlisting>
  // where x is a n x 1 matrix of doubles and f is a m x 1 matrix of doubles, where m is the number of objectives  
  // On input, the variable x contains the current point and, on output, 
  // the variable f must contain the objective function values.
  // By default, the fgoalattain uses the method of finite differences to comput the gradient of the objectives
  //
  // In order to use exact gradients, we must update the header of the 
  // objective function to :
  // <programlisting>
  //   [f,G] = objfungrad ( x )
  // </programlisting>
  // where x is a n x 1 matrix of doubles, f is a 1 x 1 matrix of doubles
  // and G is a n x 1 matrix of doubles.
  //
  // On input, the variable x contains the current point and, on output, 
  // the variable f must contain the objective function value and the variable 
  // G must contain the gradient of the objective function.
  // Furthermore, we must enable the "GradObj" option with the statement :
  // <programlisting>
  //   options = optimset("GradObj","on");
  // </programlisting>
  // This will let fmincon know that the exact gradient of the objective 
  // function is known, so that it can change the calling sequence to the 
  // objective function.
  //
  // The constraint function must have header :
  // <programlisting>
  //   [c, ceq] = confun(x)
  // </programlisting>
  // where x is a n x 1 matrix of doubles, c is a nni x 1 matrix of doubles and 
  // ceq is a nne x 1 matrix of doubles (nni : number of nonlinear inequality 
  // constraints, nne : number of nonlinear equality constraints).
  // On input, the variable x contains the current point and, on output, 
  // the variable c must contain the nonlinear inequality constraints and ceq must contain the 
  // nonlinear equality constraints.
  //
  // By default, the gradient of the constraints is computed using the method of finite differences
  // If exact gradients can be computed, we need to update the header to
  // <programlisting>
  //   [c,ceq,DC,DCeq] = confungrad(x)
  // </programlisting>
  // where x is a n x 1 matrix of doubles, c is a nni x 1 matrix of doubles, 
  // ceq is a nne x 1 matrix of doubles, DC is a n x nni matrix of doubles and 
  // DCeq is a n x nne matrix of doubles.
  // On input, the variable x contains the current point and, on output, 
  // the variable c must contain the nonlinear inequality constraint function value,
  // the variable ceq must contain the nonlinear equality constraint function value,
  // the variable DC must contain the Jacobian matrix of the nonlinear inequality constraints
  // and the variable DCeq must contain the Jacobian matrix of the nonlinear equality constraints.
  // The i-th nonlinear inequality constraint is associated to the i-th column of 
  // the matrix DC, i.e, it is stored in DC(:,i) (same for DCeq).
  // Furthermore, we must enable the "GradObj" option with the statement :
  // <programlisting>
  //   options = optimset("GradConstr","on");
  // </programlisting>
  //
  // The fgoalattain finds out the maximum value of Y for the objectives evaluated at the starting point and
  // adds that as another variable to the vector x
  // This is passed to the fmincon function to get the optimised value of Y
  // Hence, the algorithm used mainly is "ipopt" to obtain the optimum solution
  // The relations between f(x), Y, weights and goals are added as additional non-linear inequality constraints 
  // 
  // The exitflag variable allows to know the status of the optimization.
  // <itemizedlist>
  //   <listitem>exitflag=0 : Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.</listitem>
  //   <listitem>exitflag=1 : First-order optimality measure was less than options.TolFun, and maximum constraint violation was less than options.TolCon.</listitem>
  //   <listitem>exitflag=-1 : The output function terminated the algorithm.</listitem>
  //   <listitem>exitflag=-2 : No feasible point was found.</listitem>
  //   <listitem>exitflag=%nan : Other type of termination.</listitem>
  // </itemizedlist>
  //
  // The output data structure contains detailed informations about the 
  // optimization process. 
  // It has type "struct" and contains the following fields.
  // <itemizedlist>
  //   <listitem>output.iterations: the number of iterations performed during the search</listitem>
  //   <listitem>output.funcCount: the number of function evaluations during the search</listitem>
  //   <listitem>output.stepsize: an empty matrix</listitem>
  //   <listitem>output.algorithm : the string containing the name of the algorithm. In the current version, algorithm="ipopt".</listitem>
  //   <listitem>output.firstorderopt: the max-norm of the first-order KKT conditions.</listitem>
  //   <listitem>output.constrviolation: the max-norm of the constraint violation.</listitem>
  //   <listitem>output.cgiterations: the number of preconditionned conjugate gradient steps. In the current version, cgiterations=0.</listitem>
  //   <listitem>output.message: a string containing a message describing the status of the optimization.</listitem>
  // </itemizedlist>
  //
  // The lambda data structure contains the Lagrange multipliers at the 
  // end of optimization.
  // It has type "struct" and contains the following 
  // fields.
  // <itemizedlist>
  //   <listitem>lambda.lower: the Lagrange multipliers for the lower bound constraints. In the current version, an empty matrix.</listitem>
  //   <listitem>lambda.upper: the Lagrange multipliers for the upper bound constraints. In the current version, an empty matrix.</listitem>
  //   <listitem>lambda.eqlin: the Lagrange multipliers for the linear equality constraints.</listitem>
  //   <listitem>lambda.eqnonlin: the Lagrange multipliers for the nonlinear equality constraints.</listitem>
  //   <listitem>lambda.ineqlin: the Lagrange multipliers for the linear inequality constraints.</listitem>
  //   <listitem>lambda.ineqnonlin: the Lagrange multipliers for the nonlinear inequality constraints.</listitem>
  // </itemizedlist>
  // TODO : implement Display option
  // TODO : implement Algorithm option
  // TODO : implement GoalsExactAchieve option
  // TODO : implement TolFun option
  // TODO : implement TolCon option
  // TODO : implement TolX option
  // TODO : test all exitflag values
  //Examples
  //function f1 = objfun(x)
  //f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
  //f1(2)=-x(1)*x(1)-3*x(2)*x(2)
  //f1(3)=x(1)+3*x(2)-18
  //f1(4)=-x(1)-x(2)
  //f1(5)=x(1)+x(2)-8
  //endfunction
  //x0=[-1,1];
  //
  //goal=[-5,-3,-2,-1,-4];
  //weight=abs(goal)
  //gval  =
  //[- 1.229D-09  
  //- 64.        
  //- 2.         
  //- 8.         
  //  3.842D-11 ] 
  //z  =
 
  //  [4.    4.]  
  //
  //Run fgoalattain
  //[x,fval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight)
  //
  // Authors
  // Prajwala TM, Sheetal Shalini , 2015
    
    // Check number of input and output arguments
    [lhs,rhs] = argn()
    fgoalattain_checkrhs("fgoalattain", rhs, [4 6 8 10 11 12])
    fgoalattain_checklhs("fgoalattain", lhs, 1:6)
    
    // initialisation of fun
    objfun = varargin(1)
    fgoalattain_checktype("fgoalattain", objfun, "objfun", 1, "function")
    
    // initialisation of startpoint
    startpoint = varargin(2)
    //to convert startpoint to a vector if user-entered startpoint is a matrix, using linear indexing
    nrow=size(startpoint,'r')
    ncol=size(startpoint,'c')
    lin=matrix(startpoint,1,nrow*ncol)
    startpoint=lin(:)

    fgoalattain_checktype("fgoalattain",startpoint, "startpoint", 2, "constant")
    
    numvar = size(startpoint,"*")
    fgoalattain_checkvector("fgoalattain",startpoint, "startpoint", 2, numvar)
    startpoint = startpoint(:)
    
    // initialisation of goal
    goal=varargin(3)
    fgoalattain_checktype("fgoalattain",goal,"goal",3,"constant")
    
    // initialisation of weight
    weight=varargin(4)
    fgoalattain_checktype("fgoalattain",weight,"weight",4,"constant")
    
    //initialisation of A and b
    if(rhs < 5) then
        A = []
        b = []
    else
        A = varargin(5)
        b = varargin(6)
    end
    
    fgoalattain_checktype("fgoalattain", A, "A", 5, "constant")
    fgoalattain_checktype("fgoalattain", b, "b", 6, "constant")
    
    numrowA = size(A,"r")
    if(A <> []) then
        fgoalattain_checkdims("fgoalattain", A, "A", 5, [numrowA numvar])
        fgoalattain_checkvector("fgoalattain", b, "b", 6, numrowA)
        b = b(:)
    end
    
    //initialisation of Aeq and beq
    if(rhs < 7) then
        Aeq = []
        beq = []
    else
        Aeq = varargin(7)
        beq = varargin(8)
    end
    
    fgoalattain_checktype("fgoalattain", Aeq, "Aeq", 7, "constant")
    fgoalattain_checktype("fgoalattain", beq, "beq", 8, "constant")
    
    numrowAeq = size(Aeq,"r")
    if(Aeq <> []) then
        fgoalattain_checkdims("fgoalattain", Aeq, "Aeq", 7, [numrowAeq numvar])
        fgoalattain_checkvector("fgoalattain", beq, "beq", 8, numrowAeq)
        beq = beq(:)
    end
    
    // initialisation of lb and ub
    if(rhs < 9) then
        lb = []
        ub = []
    else
        lb = varargin(9)
        ub = varargin(10)
    end
    
    fgoalattain_checktype("fgoalattain", lb, "lb", 9, "constant")
    fgoalattain_checktype("fgoalattain", ub, "ub", 10, "constant")
    
    // Check dimensions of lb and ub
    if(lb <> []) then
        fgoalattain_checkvector("fgoalattain", lb, "lb", 9, numvar)
        lb = lb(:)
    end
    
    if(ub <> []) then
        fgoalattain_checkvector("fgoalattain", ub, "ub", 10, numvar)
        ub = ub(:)
    end
    
    // Initialisation of nonlinfun
     function [c,ceq] = constr(z)
            c = []
            ceq = []
     endfunction
    
    if(rhs < 11) then
        nonlinfun = constr
    else
        nonlinfun = varargin(11)
    end
    
    fgoalattain_checktype("fgoalattain", nonlinfun, "nonlinfun", 11, "function")
    
    // initialisation of constants
    defaultoptions =list("MaxIter", [1500], "CpuTime", [500], "GradObj", "OFF","GradCon", "OFF","GoalsExactAchieve",0);
    //disp(defaultoptions)
    if(rhs < 12) then
        g_attainoptions = defaultoptions
    else
        g_attainoptions = varargin(12)
    end
    
   //Flags to check whether Gradient is "ON"/"OFF" & GoalsExactAchieve is 0 or not 
   	flag1=0;
   	flag2=0;
   	flag3=0;
        gexactachieve=0	
        param=g_attainoptions
        pos_fGrad=0
        pos_cGrad=0
        pos1=6
        pos2=8
        pos3=10
        options = list(..
      		"MaxIter"     , [10000], ...
      		"CpuTime"   , [600] ...
		);
      //disp(g_attainoptions)
	//To check the User Entry for Options and storing it
   	for i = 1:(size(param))/2
       	select param(2*i-1)
    		case "MaxIter" then
          			options(2*i) = param(2*i);    //Setting the Maximum Iteration as per user entry
       		case "CpuTime" then
          			options(2*i) = param(2*i);    //Setting the Maximum CPU Time as per user entry
        	case "GradObj" then
        				if (param(2*i)=="ON") then
        					flag1=1;
        					pos_fGrad=2*i
						pos1=2*i
        					end
        				
      
        	case "GradCon" then
        				if (param(2*i)=="ON") then
        					 pos_cGrad=2*i
        				         flag2=1;
						 pos2=2*i
        					 end
        						       				      
        			
            case "GoalsExactAchieve" then 
						if(param(2*i)<>0) then
						 flag3=1 
                                                 gexactachieve=param(2*i)
						 pos3=2*i	
						 end
					         
            else
    	      	errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "fgoalattain", param(2*i-1));
    	      	error(errmsg)
    		end
   	end
   
      

    objfunval = objfun(startpoint)
    objfunval=objfunval(:)
    goal1=goal(:)
    weight1=weight(:)
     tempvar=[]
    // appending the gamma value as another variable 
    for i=1:size(objfunval,'c')
      if(weight(i)<>0) then
      tempvar(i)=((objfunval(i)-goal1(i))/weight1(i))
     // tempvar=[tempvar;((objfunval(i)-goal1(i))/weight1(i))]
      end
    end
    
    startpoint(numvar+1)=max(tempvar)
    
    if(A <> []) then 
        A = [A'; zeros(1,numrowA)]'
    end
    if(Aeq <> []) then
        Aeq = [Aeq'; zeros(1,numrowAeq)]'
    end
    if(lb <> []) then
        lb(numvar+1) = -%inf
    end
    if(ub <> []) then
        ub(numvar+1) = +%inf
    end
    
    add_ineq_with_wt = []  //to store additional constraints if not all weights are non-zero
    add_ineq_witout_wt = []
    add_equ = []
    
    // function handle defining the additional inequalities
    function temp = add_ineq(z)
        tmp = objfun(z)
        temp = []
        for i = 1:size(tmp,'r')
            if(weight(i) <> 0) then
                add_ineq_with_wt = [add_ineq_with_wt; ( (tmp(i)-goal(i))/weight(i) )]
            else
                add_ineq_witout_wt = [add_ineq_witout_wt; tmp(i)-goal(i)]
            end
        end
        temp = [temp; add_ineq_with_wt - ones(size(add_ineq_with_wt,'r'),1)*z(numvar+1)]
        temp = [temp; add_ineq_witout_wt]
    endfunction
 
 // to check if the "GoalsExactAchieve" option is not 0
    
 // function handle defining additional equalities
    function temp = add_eq(z)
        temp = []
        if(gexactachieve<> 0) then
            for i = 1:gexactachieve
                tmp = objfun(z)
                add_equ = [add_equ; tmp(i)-goal(i)]
            end
        end
        temp = [temp; add_equ]
    endfunction
    

    der_newfunc=[]
    derobjfun=[]
    der_newfunc = zeros(numvar,1)
    der_newfunc = [der_newfunc; 1]
    // function handle defining new objective function
    function newfunc = new_objfun(z)
        newfunc = z(numvar+1)
    endfunction
    
    // function handle defining derivative via finite differences
    function func = der_app(f,z)
        func = []
        nvar = size(z,'*')
        for i = 1:nvar
            t = z
            t(i) = z(i) + 10^-7
            tmp = ((f(t) - f(z))/10^-7)'
            func = [func;tmp]
        end
    endfunction
    
   // function handle to add the derivatives of the objectives
    function func = der_obj_app(z)
        func = der_app(add_ineq,z)
    endfunction
  
    // function handle to add the derivatives of additional equalities and inequalities
    function func = der_ineq_app(z)
        func = der_app(add_ineq,z)
    endfunction
    
    function func = der_eq_app(z)
        func = der_app(add_eq,z)
    endfunction
    
  // function handle to define the derivative via finite differences for non-linear constraints  
    function [dc,dceq] = der_nonlin_app(z)
        dc = []
        dceq = []
        nvar = size(z,'*')
        for i = 1:nvar
            t = z
            t(i) = z(i) + 10^-7
            [c1,ceq1] = nonlinfun(t)
            [c2,ceq2] = nonlinfun(z)
            tmpc = ((c1 - c2)/10^-7)'
            tmpceq = ((ceq1 - ceq2)/10^-7)'
            dc = [dc; tmpc]
            dceq = [dceq; tmpceq]
        end
    endfunction
    
 // function handle defining new nonlinfun function
    function [nc,nceq] = new_nonlinfun_witout_der(z)
     [tmpvar1,tmpvar2] = nonlinfun(z)
        tmpvar3 = add_ineq(z)
        tmpvar4 = add_eq(z)
        // add the objectives specified by "GoalsExactAchieve" as hard constraints to the non-linear equality constraints
        nc = [tmpvar1; tmpvar3]
        nceq = [tmpvar2; tmpvar4]
        nc=matrix(nc,size(nc,"c"),size(nc,"r"))
        nceq=matrix(nceq,size(nceq,"c"),size(nceq,"r"))
    endfunction
    
    dnc=[]
    dnceq=[]
    
    function [nc,nceq] = new_nonlinfun_with_der(z)
        [tmpvar1,tmpvar2] = nonlinfun(z)
        tmpvar3 = add_ineq(z)
        tmpvar4 = add_eq(z)
        
        nc = [tmpvar1; tmpvar3]
        nceq = [tmpvar2; tmpvar4]
        dnc = []
        dnceq = []
        
        // check if "GradCon" option is turned on
        // if "GradCon" is turned on, use it
        if(flag2==1) then
            [a,b,derc,derceq] = nonlinfun(z)
            dnc = [dnc; derc]
            dnceq = [dnceq; derceq]
            dnc = [dnc; zeros(1,size(dnc,'c'))]
            dnceq = [dnceq; zeros(1,size(dnceq,'c'))]
            // else, calculate it using finite differences
        else
            [t1,t2] = der_nonlin_app(z)
            dnc = [dnc; t1]
            dnceq = [dnceq; t2]
        end
        
        // check if "GradObj" option is turned on
        //objder = optimget(g_attainoptions,'GradObj')
        // if "GradObj" is turned on, use it
        if(flag1==1) then
            [a,derobjfun] = objfun(z)
            tempvar2der = []
            tempvar3der = []
            tempvar4der = []
            for i = 1:size(a,'r')
                if weight(i) <> 0 then
                    tmp = [derobjfun(:,i)/weight(i) ; -1]
                    tempvar3der = [ tempvar3der, tmp ] 
                else
                    tmp = [derobjfun(:,i) ; 0]
                    tempvar4der = [ tempvar4der, tmp ]
                end
            end
            
            if(gexactachieve<> 0) then
                for (i=1:gexactachieve)
                    tmp = [derobjfun(:,i) ; 0]
                    tempvar2der = [ tempvar2der, tmp ]
                end
            end
            
            dnc = [ dnc, tempvar3der, tempvar4der ]
            dnceq = [ dnceq, tempvar2der ]
            // else, calculate it using finite differences
        else
            deraddineq = der_ineq_app(z)
            deraddeq = der_eq_app(z)
            dnc = [dnc, deraddineq]
            dnceq = [dnceq, deraddeq]
            derobjfun=der_newfunc
        end
        nc=matrix(nc,size(nc,"c"),size(nc,"r"))
        nceq=matrix(nceq,size(nceq,"c"),size(nceq,"r"))
    endfunction
    
    options = g_attainoptions
   // if either of "GradObj" or "GradConstr" is turned on, turn both on 
   if(flag1==1 | flag2==1) then
       g_attainoptions(pos_fGrad)="ON"
       g_attainoptions(pos_cGrad)="ON"
       updated_nonlinfun = new_nonlinfun_with_der
    else
       g_attainoptions(pos_fGrad)="ON"
       g_attainoptions(pos_cGrad)="OFF"
       updated_nonlinfun = new_nonlinfun_witout_der
    end
    
    flag1=1
    function [c1,c2]=cGrad(x)
    c1=dnc
    c2=dnceq
    endfunction

    function fder=fGrad(x)
    fder=der_newfunc
    endfunction
    //disp(startpoint)
    //disp(objfun(startpoint))
    //disp(new_objfun(startpoint))
    //disp(goal1)
    //disp(weight1)
   // [f,g] = new_nonlinfun_witout_der(startpoint)
    //disp(f)
    //disp(g)
    g_attainoptions(1)=null()
    g_attainoptions(1)=null()
    g_attainoptions(9)=null()
    g_attainoptions(9)=null()
    g_attainoptions(pos1)="ON"
    disp(g_attainoptions)
    disp(fGrad(x0))
    if(flag2==0 & flag3==0) then
     [tmpx,tmpfval,tmpexitflag,tmpoutput,tmplambda] = fmincon  (new_objfun,startpoint,A,b,Aeq,beq,lb,ub,updated_nonlinfun,g_attainoptions,fGrad)
   else
   [tmpx,tmpfval,tmpexitflag,tmpoutput,tmplambda] = fmincon  (new_objfun,startpoint,A,b,Aeq,beq,lb,ub,updated_nonlinfun,g_attainoptions,fGrad,cGrad)
    end
    x= tmpx(1:numvar)
    attainfactor=tmpfval
    fval = objfun(tmpx)
    exitflag = tmpexitflag
    output = tmpoutput
    lambda = tmplambda
    
endfunction
