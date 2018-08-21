libname := cat(currentdir(),"/release"),libname :
savelibname := convert(libname[1],name) :
#------------------------------------------------------------------------------
### Name                  : Some tools used for determining system construction
### Input                 :
### Output                :
### Description           :
### References            :
### Date                  : Tue Jul 11 15:18:07 CEST 2006
### Last modified         :
### Implementation        :
### Error Conditions      :
### Environment Variables :
### Side Effects          :
### Other                 :
### To Do		          :

Bracket := proc(  Der1, Der2, Coord :: list(name) )

option `Copyleft 2006 by Alexandre.Sedoglavic@lifl.fr. Used by observabilitySymmetries` :

local i :: name;

	[seq(Der1(Der2(i))-Der2(Der1(i)),i=Coord)] ;

end proc :

MakeDerivation := proc(Coefficients :: list, Coord :: list(name))

option `Copyleft 2006 by Alexandre.Sedoglavic@lifl.fr. Used by observabilitySymmetries` :

local X :: name :

    subs( Diff=diff,
          unapply(add(Coefficients[i]*Diff(X,Coord[i]),i=1..nops(Coord)),
                  X,
                  proc_options = {remember,operator,arrow}
                      )
        )  :

end proc :

#------------------------------------------------------------------------------
### Name                  : FindSymmetries
### Input                 :
### Output                :
### Description           : Compute a symmetries group of the system
###			     i.e. return a basis of infinitesimal generators'
###			     coefficients vector base
###			    After some choice of symmetries type (translation, scaling)
###			  	the determining system (PDE) is reduce to kernel computation
### References            :
### Date                  : Tue Jul 11 15:18:07 CEST 2006
### Last modified         : Wed Jul 19 16:18:09 CEST 2006
### Implementation        :
### Error Conditions      :
### Environment Variables :
### Side Effects          :
### Other                 : use RandomTools[Generate], zip
### To Do		          : + change the name of this procedure because
###                         it computes kernel in Q of a matrix whose
###                         coefficients are in Q(Variables)
###                         + suppress the second parameters
###                         (Variables) because a call to indets could
###                         be used insteed



FindSymmetries := proc(	DeterminingMatrix  :: matrix,
                        Variables 	   :: list(name)
                      )

option `Copyleft 2006 by Alexandre.Sedoglavic@lifl.fr. Used by observabilitySymmetries` :

	local foo, bar, Mat :: matrix, i :: integer, 
		NbNonObsVar :: posint, primenumber :: posint,
		intrange :: range  :

    NbNonObsVar := linalg[coldim](DeterminingMatrix) :

    intrange := 1..10^5 :
    Mat := eval(DeterminingMatrix,
                zip((foo,bar)-> foo=bar,
                    Variables,
       RandomTools[Generate](list( integer(range=intrange),nops(Variables) )))
               )  :

    # At least, considered linear system should be square
    for i from 1 to ceil(nops(Variables)/linalg[rowdim](DeterminingMatrix)) do
        Mat := linalg[stackmatrix](Mat,
                                   eval(DeterminingMatrix,
                                        zip((foo,bar)-> foo=bar,
                                            Variables,
                            RandomTools[Generate](list(integer(range=intrange),
							nops(Variables))) )
                                       )
                                  )  :
    end do :
    primenumber := nextprime(10^20) :  # reduce computation time 
  # primenumber := nextprime(ceil(evalf(NbNonObsVar^(NbNonObsVar/2))*10^5)) : 
  # this choice is arbitrary and inspired by Hadamard's inequality
  # but as we do not compute determinant, we should use a much smaller bound
	
    `mod`:=mods ; # in order to replace primenumber-1 by -1
    userinfo(2,observabilitySymmetries,
	    `The computation are done modulo the following prime number`,
	    primenumber,cat(currentdir(),"/release")):

    Mat := Gausselim(Mat,'foo') mod primenumber :
    # otherwise
    # As we are looking only for kernel in Q, some supplementary
    # specialization are maybe necessary

    do
	# if there is no kernel, we return an empty set
        if foo=NbNonObsVar then return {} end if :

        Mat := linalg[stackmatrix](Mat,
                                   eval(DeterminingMatrix,
                                        zip((foo,bar)-> foo=bar,
                                            Variables,
                            RandomTools[Generate](list(integer(range=intrange),
						nops(Variables))) )
                                       )
                                  )  :

	Mat := Gausselim(Mat,'bar') mod primenumber :
	if bar < linalg[rowdim](Mat) 
	then Mat := linalg[delrows](Mat, bar+1..linalg[rowdim](Mat) ) :
	end if:

	userinfo(2,observabilitySymmetries,
		"Considering a specialized determining ",
             	cat(linalg[rowdim](Mat),"x",NbNonObsVar),
             	" matrix"
             ) :

	if bar = foo then break end if :
	foo := bar :
    end do :

	userinfo(2,observabilitySymmetries,"Start kernel computation"):
    Mat := map(convert,linalg[kernel](Mat),list) :

    if Mat={} then ERROR("Something wrong with kernel computation")
    else Mat := map(`mods`,Mat,primenumber) :
	 Mat := convert(IntegerRelations:-LLL(convert(Mat,list)),set) :
	 return Mat 
    end if :
end proc :

#------------------------------------------------------------------------------
### Name                  : observabilitySymmetries
### Input                 :
### Output                : a table res with 4 indices
###			`coordinates`,`translationcat(currentdir(),"/release")scaling` and `other`.
###	       The entrie associated to the indice coordinates define the 
###            coordinates set in whish infinitesimal generators are defined
###   	       Associated entries are sets of infinitesimal generators
###	       of symmetries groups acting on considered systems and leaving
###	       inputs and outputs invariants.
### Description           :
### References            :
### Date                  : Tue Jul 11 15:18:07 CEST 2006
### Last modified         :
### Implementation        :
### Error Conditions      :
### Environment Variables :
### Side Effects          :
### Other                 : The time is always denoted by `t`
### To Do	          :
	# create a procedure PrintAutomorphismGroup(infgen, coordinates)
	# in order to factorise code


CheckInfGen := proc (  Vect :: list, Mat :: matrix )

option `Copyleft 2006 by Alexandre.Sedoglavic@lifl.fr` :
description "This function is Used by observabilitySymmetries and check if an infinitesimal generators computed by the function FindSymmetries is a computational artefact or not" :

local res :: set ;

    res := map(Testzero,convert(evalm(Mat&*linalg[vector](Vect)),set)) :

    if nops(res)=1 then return op(res)
    elif nops(res)=2 then return `and`(op(res))
    else ERROR(`Unable to test zero in expression`,res)
    end if
end :

observabilitySymmetries := proc (   VectorField		:: list,
                                    Variables 		:: list(name),
                                    OutputRelation	:: list,
                                    Parameters		:: list(name),
                                    Input		:: list(name),
                                    OutputObsTest	:: list
                                )

option `Copyleft 2006 by Alexandre.Sedoglavic@lifl.fr` :

description "Given observability analysis done by the function observabilityTest, observabilitySymmetries computes symmetries acting on considered model and letting its inputs and outputs invariants. See ? observabilitySymmetries" :

local 	Mat :: matrix, Matt :: matrix, foo, bar, InfGenerators :: list(list),
    Var :: list(name), i :: integer, j :: integer,  t :: symbol,
    foobar, res :: table, L,DD,alpha, beta:

    userinfo(1,observabilitySymmetries,"Version 0.0 "):
    #----------------------------------------------------------------------
    # if there is nothing to do
    if OutputObsTest[1] = 0 then
        userinfo(1,observabilitySymmetries,"No symmetries"):
        return NULL ;
    end if :
    #----------------------------------------------------------------------
    # avoid maple garbage collecting message
    kernelopts(printbytes=false):
    #----------------------------------------------------------------------
    # prepare the table that stores results
    res := table([	`translation`={},
			`scaling`={},
			`other`={},
			`coordinates`=[t,op(OutputObsTest[2])]
		]) :
    userinfo(1,observabilitySymmetries,"The independant variable t represents time"):
    userinfo(1,observabilitySymmetries,"It is not supposed to be observable"):
    #----------------------------------------------------------------------
    # translation computation
    #----------------------------------------------------------------------
    # First, we compute a Jacobian matrix w.r.t. unobservable variables
    # used in the sequel
    Mat := linalg[matrix](BaurStrassen( [1,op(VectorField),op(OutputRelation)],
                                        res[`coordinates`] )	) :
    # This branche of the code  suppose that time is not observable
    #----------------------------------------------------------------------
    # Observable variables could be directly specialized before calling
    # FindSymmetries
    Var := [op(Input),t,op(Variables),op(Parameters) ] :

    #--------------------------------------------------------------------------
    userinfo(1,observabilitySymmetries,"Searching translation"):
    res[`translation`] := FindSymmetries( Mat, Var ) :
    #----------------------------------------------------------------------
    # As computation are based on specialization, spurious
    # infinitesimal generators could be computed. Thus, we must
    # check these symmetries.
    res[`translation`] := select( CheckInfGen, res[`translation`], Mat ) :

    #----------------------------------------------------------------------
    # Printing the result
    if evalb(res[`translation`]={}) then
        userinfo(1,observabilitySymmetries,"no translation found"):
    else
        for foo in res[`translation`] do
            userinfo(1,observabilitySymmetries,
         "I found a one-parameters (_lambda) group of symmetries defined by") :
            for i from 1 to nops(res[`coordinates`]) do
                if not Testzero(foo[i])  then
                    userinfo(1,observabilitySymmetries,
                             cat(
                                 res[`coordinates`][i],
                                 " ----> ",
                                 foo[i],
                                 " * ",
                                 " _lambda ",
                                 " + ",
                                 res[`coordinates`][i])
                            ):
                end if :
            end do :
            userinfo( 1,
                      observabilitySymmetries,
                      "others variables and parameters are unchanged"
                    ) :
        end do :
    end if :

    if nops(res[`translation`])=OutputObsTest[1]+1 then return res end if :
    # +1 because if system is autonomous, time translations always occurs
    #----------------------------------------------------------------------
    # scaling computation
    #--------------------------------------------------------------------------
    userinfo(1,observabilitySymmetries,"Searching dilatation"):
    # This time, we directly used definition of these systems
    L := MakeDerivation([1,op(VectorField),0$nops(Parameters)],
			Var[nops(Input)+1..-1]) :

    DD := MakeDerivation( [seq(	cat(currentdir(),"/release")(	member(foo,res[`coordinates`]),
                    			foo*cat('_',foo),
					0
				),foo=Var[nops(Input)+1..-1]
			  )],
                          Var[nops(Input)+1..-1]) :

    foobar :=  convert(Bracket(DD,L,Var[nops(Input)+1..-1]),set) ;
    foobar := foobar union convert(map(DD,OutputRelation),set) :
    Mat := linalg[matrix](BaurStrassen(convert(foobar,list), map2(cat,'_',res[`coordinates`]))):
    res[`scaling`] := FindSymmetries(Mat, Var) :

    #----------------------------------------------------------------------
    # As computation are based on specialization, spurious
    # infinitesimal generators could be computed. Thus, we must
    # check these symmetries.
    res[`scaling`] := select(CheckInfGen, res[`scaling`], Mat) :

    #----------------------------------------------------------------------
    # Printing the result
    if evalb(res[`scaling`] = {}) then
        userinfo(1,observabilitySymmetries,"no scaling found") :
    else
        for foo in res[`scaling`] do
            userinfo(1,observabilitySymmetries,
         "I found a one-parameters (_lambda) group of symmetries defined by") :
            for i from 1 to nops(res[`coordinates`]) do
                if not Testzero(foo[i])  then
                    userinfo(1,observabilitySymmetries,
                             cat(
                                 res[`coordinates`][i],
                                 " ----> ",
                                 res[`coordinates`][i],
                                 " * ",
                                 " _lambda^",
                                 foo[i])
                            ):
                end if :
            end do :
            userinfo( 1,
                      observabilitySymmetries,
                      "others variables and parameters are unchanged") :
        end do :

        res[`scaling`] := {seq( zip( (foo,bar)->foo*bar, foobar, 
						res[`coordinates`] ), 
				foobar=res[`scaling`])} :
    end if :

    if nops(res[`translation`])+nops(res[`scaling`])=OutputObsTest[1]+1
	then return res end if :

    #----------------------------------------------------------------------
    # Others symmetries computation
    #----------------------------------------------------------------------

    userinfo(1,observabilitySymmetries,"Computing more general symmetries") ;

    # Construction of determining system

    DD := MakeDerivation( [seq(	cat(currentdir(),"/release")(
					member(foo,res[`coordinates`]),
                    			add(j*alpha[cat('_',foo),cat('_',j)],
						j=res[`coordinates`])+
                         				beta[cat('_',foo)],
					0
				),foo=Var[nops(Input)+1..-1]
			  )],
                          Var[nops(Input)+1..-1]) :
	
    foobar :=  convert(Bracket(DD,L,Var[nops(Input)+1..-1]),set) ;
    foobar := foobar union convert(map(DD,OutputRelation),set) :

    foo := [seq(seq(alpha[cat('_',foo),cat('_',bar)],bar = res[`coordinates`]),
			foo = res[`coordinates`] ),
            seq(beta[cat('_',foo)],foo = res[`coordinates`])] :

    foobar := convert(remove(Testzero,foobar),list) :
    Matt := matrix(	nops(foobar),	nops(foo),
			[seq([seq(coeff(i,j),j=foo)],i=foobar)]
		):

    Matt := convert(FindSymmetries(Matt, Var),list) :

    # construct computed infinitesimal generators
    res[`other`] := {seq(
	[seq(add(foo[k]* res[`coordinates`][k-(j-1)*nops(res[`coordinates`])],
	k=1+(j-1)*nops(res[`coordinates`])..j*nops(res[`coordinates`])) 
			+ foo[nops(res[`coordinates`])^2+j],
		j=1..nops(res[`coordinates`]))]
		,foo=Matt
	      )}:

    foo := [seq(linalg[vector](foo[nops(res[`coordinates`])^2+1..-1]),
		foo=Matt)
	];
    foobar := [seq(linalg[matrix](nops(res[`coordinates`]),
				  nops(res[`coordinates`]),
				  foo[1..nops(res[`coordinates`])^2]),
		foo=Matt)
	]; 

    for i from 1 to nops(foo) do
       foobar[i] := linalg[jordan](foobar[i],'Matt') :
       Mat := matrix( nops(res[`coordinates`]),
			   nops(res[`coordinates`]),
		[
		[0$nops(res[`coordinates`])]$nops(res[`coordinates`])
		]
		) :
	for j from 1 to nops(res[`coordinates`])-1 do
		Mat[j,j] := _lambda^foobar[i][j,j] :
		if foobar[i][j,j+1]<>0 then
		Mat[j,j+1] := _lambda^Mat[j,j] : 
		end if ;
	end do ;
	Mat[j,j] := _lambda^foobar[i][j,j] :
	Mat := evalm(Matt&*Mat&*linalg[inverse](Matt)) :
	foobar[i] := evalm(Mat&*linalg[vector](res[`coordinates`]) +
			   (Mat - eval(Mat,_lambda=1))&*foo[i]
			)	:
	foobar[i] := convert(evalm(map(evalc,foobar[i]) assuming _lambda :: posint),list) :

      userinfo(1,observabilitySymmetries,
         "I found a one-parameters (_lambda) group of symmetries defined by") :
	for j from 1 to nops(res[`coordinates`]) do
             userinfo(1,observabilitySymmetries,
                             cat(
                                 res[`coordinates`][j],
                                 " ----> ",
                                 convert(foobar[i][j],string))
                            ):
	end do ;
    end do ;
	
    #----------------------------------------------------------------------
    if nops(res[`translation`])
	+ nops(res[`scaling`])
	+ nops(res[`other`]) 
	< OutputObsTest[1]
    then userinfo(1,observabilitySymmetries,
	"Sorry unable to compute all symmetries group") ;
    end if :

    #----------------------------------------------------------------------
    # set maple garbage collecting message
    kernelopts(printbytes=true):
    #----------------------------------------------------------------------

    return res :
end: # observabilitySymmetries
#------------------------------------------------------------------------------

savelib(`Bracket`,
        `MakeDerivation`,
        `FindSymmetries`,
        `CheckInfGen`,
        `observabilitySymmetries`) ;
quit :

#------------------------------------------------------------------------------
# That's all folks
