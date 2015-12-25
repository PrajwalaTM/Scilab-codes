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
function err = fgoalattain_checkvector ( funname , var , varname , ivar , nbval )
// This function checks if the inputs are vectors with the correct dimensions
  err = []
  nrows = size(var,"r")
  ncols = size(var,"c")
  if ( nrows <> 1 & ncols <> 1 ) then
    strcomp = strcat(string(size(var))," ")
    err = msprintf(gettext("%s: Expected a vector matrix for input argument %s at input #%d, but got [%s] instead."), funname, varname , ivar , strcomp );
    error(err)
  end
  if ( ( nrows == 1 & ncols <> nbval ) | ( ncols == 1 & nrows <> nbval ) ) then
    strcomp = strcat(string(size(var))," ")
    err = msprintf(gettext("%s: Expected %d entries for input argument %s at input #%d, but current dimensions are [%s] instead."), funname, nbval , varname , ivar , strcomp );
    error(err)
  end
endfunction

function err = fgoalattain_checktype ( funname , var , varname , ivar , expectedtype )
// This function checks if the input arguments are of the correct type
err = []
  if ( and ( typeof ( var ) <> expectedtype ) ) then
    strexp = """" + strcat(expectedtype,""" or """) + """"
    err = msprintf(gettext("%s: Expected type [%s] for input argument %s at input #%d, but got ""%s"" instead."),funname, strexp, varname , ivar , typeof(var) );
    error(err);
  end
endfunction

function err = fgoalattain_checkrhs ( funname , rhs , rhsset )
// This function checks if the number of input arguments is valid
err = []
  if ( and(rhs <> rhsset) ) then
    rhsstr = strcat(string(rhsset)," ")
    err = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while the number of expected input arguments should be in the set [%s]."), funname , rhs , rhsstr );
    error(err)
  end
endfunction

function err = fgoalattain_checklhs ( funname , lhs , lhsset )
// This function checks if the number of output arguments is valid
err = []
  if ( and ( lhs <> lhsset ) ) then
    lhsstr = strcat(string(lhsset)," ")
    err = msprintf(gettext("%s: Unexpected number of output arguments : %d provided while the expected number of output arguments should be in the set [%s]."), funname , lhs , lhsstr );
    error(err)
  end
endfunction

function err = fgoalattain_checkdims ( funname , var , varname , ivar , matdims )
// This function checks the dimensions of the input arguments
[lhs,rhs]=argn()
  fgoalattain_checkrhs ( "fgoalattain_checkdims" , rhs , 5 )
  fgoalattain_checklhs ( "fgoalattain_checkdims" , lhs , [0 1] )

  err = []
  if ( or ( size(var) <> matdims ) ) then
    strexp = strcat(string(matdims)," ")
    strcomp = strcat(string(size(var))," ")
    errmsg = msprintf(gettext("%s: Expected size [%s] for input argument %s at input #%d, but got [%s] instead."), funname, strexp, varname , ivar , strcomp );
    error(err)
  end
endfunction
