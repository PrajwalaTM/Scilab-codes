function [x,gval,fval,exitflag,output,lambda] = fgoalattain( varargin )


// Inputs
// x = fgoalattain(fun,x0,goal,weight) tries to make the objective functions supplied by fun attain the goals specified by goal by varying x, starting at x0, with weight specified by weight.
// x = fgoalattain(fun,x0,goal,weight,A,b) solves the goal attainment problem subject to the linear inequalities A*x ≤ b.
// x = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq) solves the goal attainment problem subject to the linear equalities Aeq*x = beq as well. Set A = [] and b = [] if no inequalities exist.
// x = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub) defines a set of lower and upper bounds on the design variables in x, so that the solution is always in the range lb ≤ x ≤ ub.
// x = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon) subjects the goal attainment problem to the nonlinear inequalities c(x) or nonlinear equality constraints ceq(x) defined in nonlcon. fgoalattain optimizes such that c(x) ≤ 0 and ceq(x) = 0. Set lb = [] and/or ub = [] if no bounds exist.
// x = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon,options) minimizes with the optimization options specified in options. Use optimoptions to set these options.


// Outputs
// [x,fval] = fgoalattain(...) returns the values of the objective functions computed in fun at the solution x.
// [x,fval,attainfactor] = fgoalattain(...) returns the attainment factor at the solution x.
// [x,fval,attainfactor,exitflag] = fgoalattain(...) returns a value exitflag that describes the exit condition of fgoalattain.
// [x,fval,attainfactor,exitflag,output] = fgoalattain(...) returns a structure output that contains information about the optimization.
// [x,fval,attainfactor,exitflag,output,lambda] = fgoalattain(...) returns a structure lambda whose fields contain the Lagrange multipliers at the solution x.


// Parameters
//   fun: a set of functions to be minimized. See below for the complete specifications.
//   x0: a nx1 or 1xn matrix of doubles, where n is the number of variables. The initial guess for the optimization algorithm.
//   gaol:Vector of values that the objectives attempt to attain. The vector is the same length as the number of objectives F returned by fun.fgoalattain attempts to minimize the values in the vector F to attain the goal values given by goal.
//   weight:A weighting vector to control the relative underattainment or overattainment of the objectives in fgoalattain. When the values of goal are all nonzero, to ensure the same percentage of under- or overattainment of the active objectives, set the weighting function to abs(goal). (The active objectives are the set of objectives that are barriers to further improvement of the goals at the solution.)
//   A: a nil x n matrix of doubles, where n is the number of variables and nil is the number of linear inequalities. If A==[] and b==[], it is assumed that there is no linear inequality constraints. If (A==[] & b<>[]), fmincon generates an error (the same happens if (A<>[] & b==[])). 
//   b: a nil x 1 matrix of doubles, where nil is the number of linear inequalities.
//   Aeq: a nel x n matrix of doubles, where n is the number of variables and nel is the number of linear equalities.  If A==[] and b==[], it is assumed that there is no linear equality constraints. If (Aeq==[] & beq<>[]), fmincon generates an error (the same happens if (Aeq<>[] & beq==[])). 
//   beq: a nel x 1 matrix of doubles, where nel is the number of linear inequalities.
//   lb: a nx1 or 1xn matrix of doubles, where n is the number of variables. The lower bound for x. If lb==[], then the lower bound is automatically set to -inf.
//   ub: a nx1 or 1xn matrix of doubles, where n is the number of variables. The upper bound for x. If lb==[], then the upper bound is automatically set to +inf.
//   nonlcon: a function, the nonlinear constraints. See below for the complete specifications.
//   x: a nx1 matrix of doubles, the computed solution of the optimization problem
//   fval: a 1x1 matrix of doubles, the function value at x
//   attainfactor: returns the attainment factor at the solution x ( the over or under-achievement of the goal ).
//   exitflag: a 1x1 matrix of floating point integers, the exit status. See below for details.
//   output: a struct, the details of the optimization process.  See below for details.
//   lambda: a struct, the Lagrange multipliers at optimum.  See below for details.
//   options: an optional struct, as provided by optimset

// Description
// Search the minimum of a constrained optimization problem specified by :
// find the minimum of Y such that 

//<latex>
//\begin{eqnarray}
//\mbox{min}_{x}    & f(x)-weight*Y \leq goal \\
//\mbox{subject to} & c(x) \leq 0 \\
//                  & c_{eq}(x) = 0 \\
//                  & Ax \leq b \\
//                  & A_{eq} x = b_{eq} \\
//                  & lb \leq x \leq ub
//\end{eqnarray}
//</latex>

// Currently, we use ipopt for the actual solver of fgoalattain.
// The objective function must have header :
// <programlisting>
//   f = objfun ( x )
// </programlisting>
// where x is a n x 1 matrix of doubles and f is a 1 x 1 matrix of doubles.
// On input, the variable x contains the current point and, on output, 
// the variable f must contain the objective function value.

  // fgoalattain uses finite differences with order 2 formulas and 
  // optimum step size in order to compute a numerical gradient of the 
  // objective function.
  // If we can provide exact gradients, we should do so since it improves 
  // the convergence speed of the optimization algorithm.
  // In order to use exact gradients, we must update the header of the 
  // objective function to :
  // <programlisting>
  //   [f,G] = objfungrad ( x )
  // </programlisting>
  // where x is a n x 1 matrix of doubles, f is a 1 x 1 matrix of doubles
  // and G is a n x 1 matrix of doubles.
  // On input, the variable x contains the current point and, on output, 
  // the variable f must contain the objective function value and the variable 
  // G must contain the gradient of the objective function.
  // Furthermore, we must enable the "GradObj" option with the statement :
  // <programlisting>
  //   options = optimset("GradObj","on");
  // </programlisting>
  // This will let fgoalattain know that the exact gradient of the objective function.

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
  
  // By default, fgoalattain uses finite differences with order 2 formulas and
  // optimum step size in order to compute a numerical gradient of the 
  // constraint function.
  // In order to use exact gradients, we must update the header of the 
  // constraint function to :
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
  // By default, fgoalattain uses a L-BFGS formula to compute an 
  // approximation of the Hessian of the Lagrangian.

  // The exitflag variable allows to know the status of the optimization.
  // <itemizedlist>
  //   <listitem>exitflag=0 : Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.</listitem>
  //   <listitem>exitflag=1 : First-order optimality measure was less than options.TolFun, and maximum constraint violation was less than options.TolCon.</listitem>
  //   <listitem>exitflag=-1 : The output function terminated the algorithm.</listitem>
  //   <listitem>exitflag=-2 : No feasible point was found.</listitem>
  //   <listitem>exitflag=%nan : Other type of termination.</listitem>
  // </itemizedlist>
  // TODO : 2 : Change in x was less than options.TolX and maximum constraint violation was less than options.TolCon. 
  // TODO : -3 : Current point x went below options.ObjectiveLimit and maximum constraint violation was less than options.TolCon. 
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
  //
  // TODO : exitflag=2 : Change in x was less than options.TolX and maximum constraint violation was less than options.TolCon. 
  // TODO : exitflag=-3 : Current point x went below options.ObjectiveLimit and maximum constraint violation was less than options.TolCon. 
  // TODO : fill lambda.lower and lambda.upper consistently. See ticket #111 : http://forge.scilab.org/index.php/p/sci-ipopt/issues/111/
  // TODO : test with A, b
  // TODO : test with Aeq, beq
  // TODO : test with ceq
  // TODO : avoid using global for ipopt_data
  // TODO : implement Display option
  // TODO : implement FinDiffType option
  // TODO : implement DerivativeCheck option
  // TODO : implement MaxIter option
  // TODO : implement OutputFcn option
  // TODO : implement PlotFcns option
  // TODO : implement TolFun option
  // TODO : implement TolCon option
  // TODO : implement TolX option
  // TODO : implement GradConstr option
  // TODO : implement GradObj option
  // TODO : implement Hessian option
  // TODO : check that the hessian output argument is Hessian of f only
  // TODO : test all exitflag values
  
  // Example
  //  we provide only the objective functions and the nonlinear constraint
  //  function : we let fmincon compute the gradients by numerical 
  //  derivatives and fgoalattain display the results. 
  // function f = objfun ( x )
  // f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
  // f1(2)=-x(1)*x(1)-3*x(2)*x(2)
  // f1(3)=x(1)+3*x(2)-18
  // f1(4)=-x(1)-x(2)
  // f1(5)=x(1)+x(2)-8
  // endfunction
  // x0 = [-1,1];
  // goal=[-5,-3,-2,-1,4]
  // weight=abs(goal)
  // [z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight)
  //Output
  //iterations: 11
  // funcCount: 12
  // stepsize: [0x0 constant]
  // algorithm: "ipopt"
  // firstorderopt: 2.633D-09
  // constrviolation: 0
  // cgiterations: 0
  // message: [3x1 string]
 //output  =
 // 1.  
 // exitflag  =
 // 1.0000000  
 // gval  =
 //  3.928D-11  
 //- 1.257D-09  
 // - 2.0000003  
 // - 8.         
 // - 63.999998  
 //z  =  4.0000001    3.9999999    1.0000000  
 
 // Authors
 // Sheetal Shalini, Prajwala TM



  
  [lhs,rhs]=argn()  //to count the no.of output & input arguments respectively
  fgoalattain_checkrhs ( "fgoalattain" , rhs , [4 6 8 10 11 12] )  //to check if number of input arguments is either 4,6,8,10,11 or 12
  fgoalattain_checklhs ( "fgoalattain" , lhs , 1 : 6 )  //to check if number of output arguments is anywhere between 1 to 6 (both inclusive)
  objfun = varargin(1)  //first input is the set of functions to be maximized
  fgoalattain_checktype("fgoalattain", objfun, "objfun", 1, "function")
  x0 = varargin(2)  //second input is the starting point
  fgoalattain_checktype("fgoalattain", x0, "x0", 2, "constant")
  nrow=size(x0,'r')  // to find no.of rows, if user enters x0 as a matrix
  ncol=size(x0,'c')  // to find no.of columns, if user enters x0 as a matrix
  sizevar=ncol*nrow;
  x0=matrix(x0,1,sizevar) // to convert it into a matrix of 1 row and no.of columns equal to the number of entries in the matrix
  nvar=size(x0,"*")
  fgoalattain_checkvector("fgoalattain", x0, "x0", 2, nvar)
  x0=x0(:)  // to convert the resulting matrix into a vector
  goal = varargin(3)  //third input is the goal vector
  weight = varargin(4)  //fourth input is the weight vector
  if ( rhs<5 ) then
    A = []
    b = []
  else
    A = varargin(5) // (A.x)<=b
    b = varargin(6)
  end
  fgoalattain_checktype("fgoalattain", A, "A", 5, "constant")
  fgoalattain_checktype("fgoalattain", b, "b", 6, "constant")

  // Check if A and b of proper dimensions
    if(A <> [] & b == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is empty, but the column vector b is not empty"), "fgoalattain", 5, 6)
        error(errmsg)
    end
    
    if(A == [] & b <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is not empty, but the column vector b is empty"), "fgoalattain", 5, 6)
        error(errmsg)
    end
    
    numrowA = size(A,"r")
    if(A <> []) then
        fgoalattain_checkdims("fgoalattain", A, "A", 5, [numrowA nvar])
        fgoalattain_checkvector("fgoalattain", b, "b", 6, numrowA)
        b = b(:)
    end
  if ( rhs<7 ) then
    Aeq = []
    beq = []
  else
    Aeq = varargin(7)  // (Aeq.x)=beq
    beq = varargin(8)
  end
  fgoalattain_checktype("fgoalattain", Aeq, "Aeq", 7, "constant")
  fgoalattain_checktype("fgoalattain", beq, "beq", 8, "constant")
    
    // Check if Aeq and beq of proper dimensions
    if(Aeq <> [] & beq == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is empty, but the column vector beq is not empty"), "fgoalattain", 7, 8)
        error(errmsg)
    end
    
    if(Aeq == [] & beq <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is not empty, but the column vector beq is empty"), "fgoalattain", 7, 8)
        error(errmsg)
    end
    
    numrowAeq = size(Aeq,"r")
    if(Aeq <> []) then
        fgoalattain_checkdims("fgoalattain", Aeq, "Aeq", 7, [numrowAeq nvar])
        fgoalattain_checkvector("fgoalattain", beq, "beq", 8, numrowAeq)
        beq = beq(:)
    end
  if(rhs<9) then // if lb and ub are not provided, declare as empty
   lb=[]
   ub=[]
  else
   lb=varargin(9)  // lb<=x<=ub
   ub=varargin(10)
  end
  fgoalattain_checktype("fgoalattain", lb, "lb", 9, "constant")
  fgoalattain_checktype("fgoalattain", ub, "ub", 10, "constant")
    
    // Check dimensions of lb and ub
    if(lb <> []) then
        fgoalattain_checkvector("fgoalattain", lb, "lb", 9, nvar)
        lb = lb(:)
    end
    
    if(ub <> []) then
        fgoalattain_checkvector("fgoalattain", ub, "ub", 10, nvar)
        ub = ub(:)
    end
  if ( rhs<11 ) then

    nonlcon = nonlinfun  // if nonlinfun is not provided, declare as empty
  else
    nonlcon = varargin(11)
  end
  fgoalattain_checktype("fgoalattain", nonlinfun, "nonlinfun", 11, "function")
  defaultoptions = optimset1("fgoalattain") // the default options 
  if ( rhs<12) then
    options = defaultoptions  // if options not provided, declare as default options
  else
    options=varargin(12)
  end
  temp=[]
  objfunval = objfun(x0)
  // converting objfunval, goal and weight into column vectors
  objfunval1=objfunval(:) 
  goal1=goal(:)
  weight1=weight(:)
  for i=1:size(objfunval1,'c')
  if(weight1(i) <> 0)
  temp(i)=(objfunval1(i)-goal1(i))/weight1(i)
  end
  end

    x0(nvar+1) = max(temp) // defining a new variable x(n+1),and initialising it to gamma (the function to be minimised)
    if(A <> []) then 
        A = [A'; zeros(1,numrowA)]'
    end
    if(Aeq <> []) then
        Aeq = [Aeq'; zeros(1,numrowAeq)]'
    end
    if(lb <> []) then
        lb(nvar+1) = -%inf
    end
    if(ub <> []) then
        ub(nvar+1) = +%inf
    end

   function [c,ceq] = nonlinfun(x)  
   // fmincon library of scilab gives error when 'c' component of nonlinfun empty
   // add a trivial case of -5 <= 0 to c to bypass this error
    c = [-5]
    ceq = []
  endfunction

  fgoalattain_checktype("fgoalattain", nonlinfun, "nonlinfun", 11, "function")

  tmp=[]
  add_ineq_with_wt = []  //to store additional constraints if not all weights are non-zero
  add_ineq_witout_wt = []
  add_equ = []
    
  // function handle defining the additional inequalities
  function temp = add_ineq(z)
        tmp = objfun(z)
        temp = []
        for i = 1:size(tmp,'r')
            if(weight1(i) <> 0) then
                add_ineq_with_wt = [add_ineq_with_wt; ( (tmp(i)-goal1(i))/weight1(i) )]
                //tmpvar(i)=((abc(i)-goal1(i))/weight1(i))
            else
                add_ineq_witout_wt = [add_ineq_witout_wt; tmp(i)-goal1(i)]
            end
        end
        temp = [temp; add_ineq_with_wt - ones(size(add_ineq_with_wt,'r'),1)*z(nvar+1)]
        temp = [temp; add_ineq_witout_wt]
    endfunction

function temp = add_eq(z)
        temp = []
        if(goals<> 0) then
            for i = 1:goals
                tmp = objfun(z)
                add_equ = [add_equ; tmp(i)-goal1(i)]
            end
        end
        temp = [temp; add_equ]
    endfunction
  
  function [newfunc,der_newfunc] = new_objfun(z) // function handle defining new objective function
        newfunc = z(nvar+1)
        der_newfunc = zeros(nvar,1)
        der_newfunc = [der_newfunc; 1]
  endfunction

function func = der_app(f,z) 
    // function handle defining derivative via finite differences
    // f'(x) = ( f(x+e) - f(x) )/e where e = 10^-7
        func = []
        nvar = size(z,'*')
        for i = 1:nvar
            t = z
            t(i) = z(i) + 10^-7
            tmp = ((f(t) - f(z))/10^-7)'
            func = [func;tmp]
        end
    endfunction
    
    function func = der_obj_app(z)
        func = der_app(add_ineq,z)
    endfunction

function [dc,dceq] = der_nonlin_app(z) // finds the derivative of non-linear function
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

 // goals=optimget(options,'GoalsExactAchieve')

  function [nc,nceq] = new_nonlinfun_witout_der(z) 
  // function handle defining new nonlinfun function
  // this nonlinfun function passed when the gradient feature is off
        [tmpvar1,tmpvar2] = nonlinfun(z)
        tmpvar3 = add_ineq(z)
	if(goals <> 0 )  then
	abc=objfun(z)
	for i=1:goals
	tmpvar2=[tmpvar2; abc(i)-goal1(i)]
	end
	end
        nc = [tmpvar1; tmpvar3;tmp]
        nceq = tmpvar2
    endfunction
  
  dnc=[]
  dnceq=[]
  function [nc,nceq] = new_nonlinfun_with_der(z)
        [tmpvar1,tmpvar2] = nonlinfun(z)
        tmpvar3 = add_ineq(z)
        tmpvar4 = add_eq(z)
        
        nc = [tmpvar1; tmpvar3]
        nceq = [tmpvar2; tmpvar4]
        //dnc = []
        //dnceq = []
        
        // check if "GradConstr" option is turned on
        constrder = optimget(options,'GradConstr')
        // if "GradConstr" is turned on, use it
        if(constrder == 'on') then
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
        objder = optimget(options,'GradObj')
        // if "GradObj" is turned on, use it
        if(objder == 'on') then
            [a,derobjfun] = objfun(z)
            tempvar2der = []
            tempvar3der = []
            tempvar4der = []
            for i = 1:size(a,'r')
                if weight1(i) <> 0 then
                    tmp = [derobjfun(:,i)/weight1(i); -1]
                    tempvar3der = [ tempvar3der, tmp ] 
                else
                    tmp = [derobjfun(:,i) ; 0]
                    tempvar4der = [ tempvar4der, tmp ]
                end
            end
            
            if(goals<> 0) then
                for (i=1:goals)
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
        end
    endfunction
	
	function [dnc1,dnceq1]=cGrad(x)
	dnc1=dnc
	dnceq1=dnceq
	endfunction

    if(optimget(options,'GradObj') == 'on' | optimget(options,'GradConstr') == 'on') then
        options = optimset1(options,'GradObj','on','GradConstr','on')
	// pass new_nonlinfun_with_der if any gradient option is on
        updated_nonlinfun = new_nonlinfun_with_der
    else
	// pass new_nonlinfun_witout_der if both gradient options are off
		options = optimset1(options,'GradObj','on','GradConstr','off')
        updated_nonlinfun = new_nonlinfun_witout_der
    end

// The non-linear equality ceq(x)=0 and inequality c(x)<=0 constraints
//    disp(x0)
//    disp(A)
//    disp(b)
//    disp(Aeq)
//    disp(beq)
//    disp(lb)
//    disp(ub)
//    disp(objfun(x0))
//    [q,w] = nonlinfun(x0)
//    disp(q)
//    disp(w)
//    [q,w] = new_objfun(x0)
//    disp(q)
//    disp(w)
//    [q,w] = new_nonlinfun_witout_der(x0)
//    disp(q)
//    disp(w)
//    [q,w,c,d] = new_nonlinfun_with_der(x0)
//    disp(q)
//    disp(w)
//    disp(c)
//    disp(d)
    
    [x,fval,exitflag,output,lambda] = fmincon(new_objfun,x0,A,b,Aeq,beq,lb,ub,updated_nonlinfun,options,fGrad,cGrad)
   x=x(1:nvar)
   gval=objfun(x)
  //if(goals <> 0) then
  //gval=gsort(gval,'g','d')
 // end
endfunction

