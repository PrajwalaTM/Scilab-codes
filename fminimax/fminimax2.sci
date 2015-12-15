function errmsg = fminimax_checkrhs ( funname , rhs , rhsset )
 errmsg = []
  if ( and(rhs <> rhsset) ) then
    rhsstr = strcat(string(rhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while the number of expected input arguments should be in the set [%s]."), funname , rhs , rhsstr );
    error(errmsg)
  end
endfunction

function [c,ceq] = nonlinfun(x)
    c = [-5]
    
    ceq = []
endfunction

function errmsg = fminimax_checklhs ( funname , lhs , lhsset )
errmsg = []
  if ( and ( lhs <> lhsset ) ) then
    lhsstr = strcat(string(lhsset)," ")
    errmsg = msprintf(gettext("%s: Unexpected number of output arguments : %d provided while the expected number of output arguments should be in the set [%s]."), funname , lhs , lhsstr );
    error(errmsg)
  end
endfunction



