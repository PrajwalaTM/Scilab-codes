function [x,fval,maxobjfun,exitflag,output,lambda] = fminimax( varargin )
//x = fminimax(fun,x0)
//x = fminimax(fun,x0,A,b)
//x = fminimax(fun,x0,A,b,Aeq,beq)
//x = fminimax(fun,x0,A,b,Aeq,beq,lb,ub)
//x = fminimax(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
//x = fminimax(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
//[x,fval] = fminimax(...)
//[x,fval,maxfval] = fminimax(...)
//[x,fval,maxfval,exitflag] = fminimax(...)
//[x,fval,maxfval,exitflag,output] = fminimax(...)
//[x,fval,maxfval,exitflag,output,lambda] = fminimax(...)
//fun : The function to be minimized;accepts a vector x and returns a vector F, the objective functions evaluated at x.
//x0: a nx1 or 1xn matrix of doubles, where n is the number of variables. The initial guess for the optimization algorithm.
//A: a nil x n matrix of doubles, where n is the number of variables and nil is the number of linear inequalities.
//b: a nil x 1 matrix of doubles, where nil is the number of linear inequalities.
//Aeq: a nel x n matrix of doubles, where n is the number of variables and nel is the number of linear equalities.
//beq: a nel x 1 matrix of doubles, where nel is the number of linear inequalities.
//lb: a nx1 or 1xn matrix of doubles, where n is the number of variables. The lower bound for x.
//ub: a nx1 or 1xn matrix of doubles, where n is the number of variables. The upper bound for x.
//nonlcon: a function returning the nonlinear equality and inequality constraints
//options: an optional struct, as provided by optimset
  [lhs,rhs]=argn()
//to check the validity of the number of input and output arguments
  fminimax_checkrhs ( "fminimax" , rhs , [2 4 6 8 9 10] )
  fminimax_checklhs ( "fminimax" , lhs , 1 : 6 )
  objfun = varargin(1)
  x0 = varargin(2)
//to convert x0 to a vector if user-entered x0 is a matrix, using linear indexing
 nrow=size(x0,'r')
 ncol=size(x0,'c')
 lin=matrix(x0,1,nrow*ncol)
 x0=lin(:)

//to set either default values or input argument values to the variables depending on the number of input arguments
  if(rhs<3) then
    A = []
    b = []
  else
    A = varargin(3)
    b = varargin(4)
  end
  if ( rhs<5 ) then
    Aeq = []
    beq = []
  else
    Aeq = varargin(5)
    beq = varargin(6)
  end
  if ( rhs<7 ) then
    lb = []
    ub = []
  else
    lb = varargin(7)
    ub = varargin(8)
  end
  if ( rhs<9 ) then

    nonlcon = nonlinfun
  else
    nonlcon = varargin(9)
  end
  defaultoptions = optimset("fminimax")
  if ( rhs<10 ) then
    options = defaultoptions
  else
    options = varargin(10)
  end
    
// function that returns the maximum of the set of objectives in objfun(z)    
    function f2 = maxobjfun(z)
        t = objfun(z)
        f2 = max(t)
    endfunction
//Using fmincon to optimise the maximum of set of objectives
    [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(maxobjfun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
endfunction

