libname := cat(currentdir(),"/release"),libname :
savelibname := convert(libname[1],name) :
#------------------------------------------------------------------------------
### Name                  : BaurStrassen
### Input                 : A maple expression and a list of indeterminates
### Output                : The list of the partial derivative of expression 
###			    w.r.t. the indeterminates.
### Description           : 
### References            : @Misc{Morgenstern:1984,
###  			    author = {Jacques Morgenstern},
###  			    title =  {How to compute fast a function and all 
###				      its derivatives: a variation on the 
###				      theorem of Baur-Strassen},
###			    month =  may,
###			    year =   1984}
###			    
### Date                  : Thu Mar 23 09:29:38 2000
### Last modified         : 
### Implementation        : recursive
### Error Conditions      : errors if expr is not rational. 
### Environment Variables : 
### Side Effects          : 
### Other                 : use add 

BaurStrassen := proc(expr, ind ::{symbol, list(symbol)})

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. See ? BaurStrassen" :

    local output, doderivation ;
	
    doderivation := proc(expr, ind::{symbol, list(symbol)})
        option remember;	local tmp, x, i, _stck;
	
            if not type(expr, numeric) then
                if   type(expr, symbol) then map2(diff, expr, ind)
		elif type(expr,list) then map(procname,expr,ind)
                elif type(expr, `^`) then
                     map(unapply(
                     op(2, expr)*op(1, expr)^(op(2, expr) - 1)*x, x),
                     procname(op(1, expr), ind))
                elif type(expr, `+`) then
                     tmp := map(procname, [op(expr)], ind);
		     [seq(add(x[i], x = tmp),i=1..nops(ind))]
                elif type(expr, `*`) then
                     tmp   := procname(op(1,expr), ind)			   :
		     i     := subsop(1=1,expr)				   :
		     _stck := map(unapply(i*x,x),tmp) 	  	           :
	             tmp   := map(unapply(op(1,expr)*x,x),procname(i,ind)) :
		     [seq(_stck[i]+tmp[i],i=1..nops(ind))]
                else ERROR(`the expression is not rational`)
                fi
            else [0 $ nops(ind)]
            fi
        end: # doderivation
	
    output 		:= doderivation(expr, ind)		:
    doderivation 	:= subsop(4 = NULL, op(doderivation))	:
    output

end : # BaurStrassen
#------------------------------------------------------------------------------
 
savelib(`BaurStrassen`) :
quit :

