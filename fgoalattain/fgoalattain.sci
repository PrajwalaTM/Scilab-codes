function [x,fval,attainfactor,exitflag,output,lambda] = fgoalattain(varargin)
    // solves a multiobjective goal attainment problem
    //Calling sequence
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
    //  startpoint: a nx1 or 1xn matrix of doubles, where n is the number of variables.
    //      The initial guess for the optimization algorithm
    //  A: a nil x n matrix of doubles, where n is the number of variables and
    //      nil is the number of linear inequalities. 
    //      
    //  b: a nil x 1 matrix of doubles, where nil is the number of linear
    //      inequalities
    //  Aeq: a nel x n matrix of doubles, where n is the number of variables
    //      and nel is the number of linear equalities. 
    //  beq: a nel x 1 matrix of doubles, where nel is the number of linear
    //      equalities
    //  lb: a nx1 or 1xn matrix of doubles, where n is the number of variables. 
    //      The lower bound for x. If lb==[], then the lower bound is 
    //      automatically set to -inf
    //  ub: a nx1 or 1xn matrix of doubles, where n is the number of variables. 
    //      The upper bound for x. If ub==[], then the upper bound is 
    //      automatically set to +inf
    //  nonlinfun: a function, the nonlinear constraints
    //
    // Output Parameters
    //  x: a nx1 matrix of doubles, the computed solution of the optimization problem
    //  fval: a vector of doubles, the value of functions at x
    //  attainfactor: The amount of over- or underachievement of the goals,γ at the solution.
    //  exitflag: a 1x1 matrix of floating point integers, the exit status
    //  output: a struct, the details of the optimization process
    //  lambda: a struct, the Lagrange multipliers at optimum
    
    // Check number of input and output arguments
    [lhs,rhs] = argn()
    apifun_checkrhs("fgoalattain", rhs, [4 6 8 10 11 12])
    apifun_checklhs("fgoalattain", lhs, 1:6)
    
    // initialisation of fun
    objfun = varargin(1)
    apifun_checktype("fgoalattain", objfun, "objfun", 1, "function")
    
    // initialisation of startpoint
    startpoint = varargin(2)
    apifun_checktype("fgoalattain",startpoint, "startpoint", 2, "constant")
    
    numvar = size(startpoint,"*")
    apifun_checkvector("fgoalattain",startpoint, "startpoint", 2, numvar)
    startpoint = startpoint(:)
    
    // initialisation of goal
    goal=varargin(3)
    apifun_checktype("fgoalattain",goal,"goal",3,"constant")
    
    // initialisation of weight
    weight=varargin(4)
    apifun_checktype("fgoalattain",weight,"weight",4,"constant")
    
    //initialisation of A and b
    if(rhs < 5) then
        A = []
        b = []
    else
        A = varargin(5)
        b = varargin(6)
    end
    
    apifun_checktype("fgoalattain", A, "A", 5, "constant")
    apifun_checktype("fgoalattain", b, "b", 6, "constant")
    
    numrowA = size(A,"r")
    if(A <> []) then
        apifun_checkdims("fgoalattain", A, "A", 5, [numrowA numvar])
        apifun_checkvector("fgoalattain", b, "b", 6, numrowA)
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
    
    apifun_checktype("fgoalattain", Aeq, "Aeq", 7, "constant")
    apifun_checktype("fgoalattain", beq, "beq", 8, "constant")
    
    numrowAeq = size(Aeq,"r")
    if(Aeq <> []) then
        apifun_checkdims("fgoalattain", Aeq, "Aeq", 7, [numrowAeq numvar])
        apifun_checkvector("fgoalattain", beq, "beq", 8, numrowAeq)
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
    
    apifun_checktype("fgoalattain", lb, "lb", 9, "constant")
    apifun_checktype("fgoalattain", ub, "ub", 10, "constant")
    
    // Check dimensions of lb and ub
    if(lb <> []) then
        apifun_checkvector("fgoalattain", lb, "lb", 9, numvar)
        lb = lb(:)
    end
    
    if(ub <> []) then
        apifun_checkvector("fgoalattain", ub, "ub", 10, numvar)
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
    
    apifun_checktype("fgoalattain", nonlinfun, "nonlinfun", 11, "function")
    
    // initialisation of constants
    defaultoptions = optimset("fgoalattain") // CHANGE HERE
    if(rhs < 12) then
        g_attainoptions = defaultoptions
    else
        g_attainoptions = varargin(12)
    end
    
    objfunval = objfun(startpoint)
    objfunval=objfunval(:)
    goal1=goal(:)
    weight1=weight(:)
    
    for i=1:size(objfunval,'c')
      tempvar(i)=((objfunval(i)-goal(i))/weight(i))
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
    
    // function handle defining the additional inequalities
    function temp = add_ineq(z)
        abc=objfun(z)
        for i=1:size(abc,'r') // CHANGE HERE : earlier size(abc,'c')
         tmpvar(i)=((abc(i)-goal1(i))/weight1(i)) // CHANGE HERE
        end
        tmpvar=tmpvar(:)
        temp = tmpvar - ones(size(tmpvar,'r'),1)*z(numvar+1)
    endfunction
    
    // function handle defining new objective function
    function [newfunc,der_newfunc] = new_objfun(z)
        newfunc = z(numvar+1)
        der_newfunc = zeros(numvar,1)
        der_newfunc = [der_newfunc; 1]
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
    
    function func = der_obj_app(z)
        func = der_app(add_ineq,z)
    endfunction
    
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
        nc = [tmpvar1; tmpvar3]
        nceq = tmpvar2
    endfunction
    
    disp(startpoint)
    disp(objfun(startpoint))
    disp(new_objfun(startpoint))
    disp(goal1)
    disp(weight1)
    [f,g] = new_nonlinfun_witout_der(startpoint)
    disp(f)
    disp(g)
    
   [tmpx,tmpfval,tmpexitflag,tmpoutput,tmplambda] = fmincon(new_objfun,startpoint,A,b,Aeq,beq,lb,ub,new_nonlinfun_witout_der,g_attainoptions)
    x= tmpx(1:numvar)
    attainfactor=tmpfval
    fval = objfun(tmpx)
    exitflag = tmpexitflag
    output = tmpoutput
    lambda = tmplambda
    
endfunction
