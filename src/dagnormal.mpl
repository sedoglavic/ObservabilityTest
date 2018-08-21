libname := cat(currentdir(),"/release"),libname :
#------------------------------------------------------------------------------
### Name                  : dagnormal
### Input                 : A rational expression f.
### Output                : A list L such that:
###			 	- L[1] is a dag coding a numerator of f;
###				- L[2] is a dag coding a denominator of f.
### Description           : 
### References            : 
### Date                  : Wed Mar 22 16:00:57 2000
### Last modified         : 
### Implementation        : recursive
### Error Conditions      : errors if expr is not rational.
### Environment Variables : 
### Side Effects          : 
### Other                 : use mul

dagnormal := proc(expr)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. See ? dagnormal" :

	local output, donormal;
	
    donormal := proc(expr)
        option remember;  local tmp, z;
            if       type(expr, symbol) or type(expr, numeric) then [expr, 1]
		elif type(expr,list) then map(procname,expr)
		elif type(expr, `^`) then  
		     if 0 < op(2, expr) then 
		             tmp := procname(op(1, expr));
			     [tmp[1]^op(2, expr), tmp[2]^op(2, expr)]
			elif op(2, expr) = 0 then [1, 1]
			else tmp := procname(op(1, expr));
			     [tmp[2]^(-op(2, expr)), tmp[1]^(-op(2, expr))]
		     fi
                elif type(expr, `*`)  then 
		     tmp := map(procname, [op(expr)]) :
		     [mul(z[1], z = tmp), mul(z[2], z = tmp)]
		elif type(expr, `+`)  then 
		     tmp :=  procname(op(1,expr)) 	: 
		     z   := procname(subsop(1=0,expr))	:
		     [tmp[1]*z[2]+tmp[2]*z[1],tmp[2]*z[2]] 
		else ERROR(`the expression is not rational`)
            fi
        end; # donormal 
	
    output 	:= donormal(expr)			:
    donormal 	:= subsop(4 = NULL, op(donormal))	:
    output

end : # dagnormal
#------------------------------------------------------------------------------

savelib(`dagnormal`) :
quit :

