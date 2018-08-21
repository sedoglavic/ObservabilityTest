## ObservabilityTest - version 0.0
##
## Title	:	Implementation of a probabilistic algorithm to test 
##			local algebraic	observability in polynomial time.
##
## Created	:	Fri Jun  9 09:32:33 2000
## Author	: 	Alexandre Sedoglavic  <sedoglav@gage.polytechnique.fr>
##				<http://medicis.polytechnique.fr/~sedoglav>
## Description: 
##
##			       	diff(Theta,t) =  0	
##				  diff(X,t)   =  F(X,Theta,U),
##				       Y      =  G(X,Theta,U).
##
## References: 
##		@TechReport{Sedoglavic:2000:oct,
##  		author =	 {Alexandre Sedoglavic},
##	  	title =	         {A probabilistic algorithm to test local 
##				  algebraic observability in polynomial time},
##  		institution =	 {GAGE laboratory},
##  		year =	 	 2000,
##  		type =	 	 {Manuscript},
##  		address =	 {Available at {\sf
##                  		 http://www.gage.polytechnique.fr/notes/}},
##  		number =	 {2000-13},
##  		month =	 oct
##		}
##
## This package provides the internal functions :	
##
## Preprocessing	- ConstructLinearVariationalSystem	;
##			- ConstructOutput			;
##			- System2SLP				;
##
## Power Series Expansion 
##			- InverseMatrixSeries 			;
##			- ExpMatrixSeries			;
##			- HomogeneousResolution			;
##			- PreprocessingForVariationOfConstants 	;
##			- SecondMembersEvaluation		;
##			- VariationOfConstants			;
##
## Postprocessing	- ObservabilityAnalysis			;
##
## This package provides the external function :
##	
## Main procedure	- ObservabilityTest			.	
##
## Used procedure :
## 			for the construction of slp : 
##				BaurStrassen, dagnormal,
##				optimize, makeproc,
## 				ranpoly (for inputs).
##			for power series expansion : linear algebra (linalg);
##				mul, add, diff, series. 
## 			for modular computation:  modp, Gausselim, 
##				nextprime, randomize, rand.
##
## To do	
##		- change t in _t and verify if _t is affected.
##		- verify that the separant does not vanish
##		- a procedure for power series expansion of the
##		  solution of a system of ODE
##		- exhibit a group of transformation.
##		- the code is designed for one output. For many output,
##		  the output evaluation and the rank computation have to be 
##		  done in the main loop of the procedure ObservabilityTest.
##		- make a type for output of ObservabilityTest	
##
##
#______________________________________________________________________________

libname := cat(currentdir(),"/release"),libname :

		#---------------------------------------#	
		# 		Preprocessing		#
		#---------------------------------------#	

#-----------------------------------------------------------------------------#
### Name                  : ConstructLinearVariationalSystem
### Input                 : System::list, variables::list, parameters:: list
### Output 		  : a list such that:
###				- list[1] is a matrix composed of :
###					- System;
###					- Inverse of A;
###					- Inverse of A &* B (logarithmic 
###					derivative);
###					- partial derivatives of System w.r.t.
###					initial	conditions and parameters 
###					(also called sensitivity matrix).
###				- list[2] is a list of arguments used 
###				in list[1].
###
### Description           : 
### References            : 
### Date                  : Mon Aug 14 15:12:59 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : 
### Used procedure	  : BaurStrassen, dagnormal, linalg.
### Side Effects          : 
### Other                 : 
				    
ConstructLinearVariationalSystem := proc (	System     :: list,
						variables  :: list(name),
						parameters :: list(name)	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local Output, Stck,z,i, NbVariables,
	  dSysOndVd, InvdSysOndVd, dSysOndV, LogarithmicDerivative,
	  dSysOndP,  SLPArguments, VarSensitivityMatrix, VarSensitivityMatrixd,
	  variablesd:

    #--------------------------------------------------------------------------
    # Some shortcuts
      NbVariables := nops(variables) :
      variablesd := map(proc(i) parse(cat(i,`d`)) end,variables) ;
      VarSensitivityMatrixd := array(	1..NbVariables,
					1..NbVariables+nops(parameters)) :
      VarSensitivityMatrix  := array(	1..NbVariables,
					1..NbVariables+nops(parameters)) :
    #--------------------------------------------------------------------------
    # The arguments of the slp
      SLPArguments 	:= [op(variablesd),op(variables)]	: 
    #--------------------------------------------------------------------------
    # Convert System in a vector field
      Stck := [seq(variablesd[z] - System[z],z = 1..NbVariables)] 	:
    # Convert System in a polynomial formulation
      Stck := dagnormal(Stck) 						:
    # and keep just the numerators
      Stck := map(unapply(i[1],i),Stck) 				:
    #--------------------------------------------------------------------------
    # Compute the partial derivative w.r.t. variablesd, variables, parameters.
    z := linalg[matrix](BaurStrassen(Stck,[op(SLPArguments),op(parameters)])) :
    #--------------------------------------------------------------------------
    # Two submatrix of these partial derivatives.
      dSysOndV := linalg[submatrix](	z,	1..NbVariables,
				        NbVariables+1..2*NbVariables	) :
      dSysOndVd := linalg[submatrix](z,1..NbVariables,1..NbVariables) :
    #--------------------------------------------------------------------------
    # Construction of the inverse of dSysOndVd
      InvdSysOndVd := linalg[matrix](NbVariables,NbVariables,0) :
    # Remember that this is a diagonal matrix by construction.
      for i from 1 to NbVariables do InvdSysOndVd[i,i] := 1/dSysOndVd[i,i] od:
    #--------------------------------------------------------------------------
    # Construction of the logarithmic derivative
      LogarithmicDerivative := &*(InvdSysOndVd,dSysOndV) :
    #--------------------------------------------------------------------------
    # Construction of the output
      Output := linalg[matrix](nops(Stck),1,Stck) 			 :
      Output := linalg[augment](Output,InvdSysOndVd,LogarithmicDerivative) :
    #--------------------------------------------------------------------------
      if parameters<>[] then
	 dSysOndP := linalg[submatrix](	z,	1..NbVariables,
					2*NbVariables+1..linalg[coldim](z)  ) :
    # We use the << sensitivity >> matrix w.r.t. parameters
    # the << sensitivity >> matrix w.r.t. parameters is a solution of the 
    # following system:
	 Output	  := linalg[augment](	Output,
		    	dSysOndVd&*linalg[submatrix](   VarSensitivityMatrixd,
						    1..NbVariables,
						    1..nops(parameters)
						) + 
			dSysOndV&*linalg[submatrix](   VarSensitivityMatrix,
						    1..NbVariables,
						    1..nops(parameters)
						) + dSysOndP    ) :
      fi:
    #--------------------------------------------------------------------------
    # We use the << sensitivity >> matrix w.r.t. the initial condition
    # the << sensitivity >> matrix w.r.t. initial conditions is a solution 
    # of the following system:
      Output	 := linalg[augment](Output,
			dSysOndVd&*linalg[submatrix](   VarSensitivityMatrixd,
					    1..NbVariables,
			    nops(parameters)+1..NbVariables+nops(parameters)
						) + 
			dSysOndV&*linalg[submatrix](   VarSensitivityMatrix,
					    1..NbVariables,
			    nops(parameters)+1..NbVariables+nops(parameters)
						)) :
      SLPArguments := [	op(SLPArguments),
			VarSensitivityMatrixd,VarSensitivityMatrix] 	:
    #--------------------------------------------------------------------------
      Output 	   := evalm(Output) 			:
    # Output: the system and the argument necessary to evaluate it.

      [Output,SLPArguments]
    
end: # ConstructLinearVariationalSystem		    
#______________________________________________________________________________



#-----------------------------------------------------------------------------#
### Name                  : ConstructOutput
### Input                 : System::list, variables::list, parameters::list,
### Output                : a list such that:
###				- list[1] is a matrix composed of
###					- partial derivatives of System w.r.t.
###					initial	conditions and parameters 
###					(also called sensitivity matrix).
###				- list[2] is a list of arguments used 
###				in list[1].
### Description           : 
### References            : 
### Date                  : Mon Aug 14 15:12:59 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : 
### Used procedure	  : BaurStrassen, linalg.
### Side Effects          : 
### Other                 : 
				    
ConstructOutput := proc (	System	   :: list,
				variables  :: list, 
				parameters :: list	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :
			
    local Stck,dSysOndV,dSysOndP,
	  NbVariables, SLPArguments, VarSensitivityMatrix :

    #--------------------------------------------------------------------------
    # Some shortcuts
      NbVariables := nops(variables) :
    #--------------------------------------------------------------------------
      VarSensitivityMatrix := array(	1..NbVariables,
					1..NbVariables+nops(parameters)) :
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # The arguments of the slp
      SLPArguments 	:= variables	: 
    #--------------------------------------------------------------------------
    # Compute the partial derivative w.r.t. variables, parameters.
   Stck := linalg[matrix](BaurStrassen(System,[op(variables),op(parameters)])):
    #--------------------------------------------------------------------------
    # Two submatrix
      dSysOndV := linalg[submatrix](Stck,1..nops(System),1..NbVariables) :
    #--------------------------------------------------------------------------
    # Construction of the << sensitivity >> matrix w.r.t. parameters
    if 	parameters<>[] then
	dSysOndP := linalg[submatrix](	Stck,1..nops(System),
					NbVariables+1..linalg[coldim](Stck)) :
			
    	# We use the  << sensitivity >> matrix w.r.t. parameters
	Stck:= evalm(dSysOndV&*linalg[submatrix](   VarSensitivityMatrix,
						    1..NbVariables,
						    1..nops(parameters)
						) + dSysOndP) 		 :

    #--------------------------------------------------------------------------
    # We use now the << sensitivity >> matrix w.r.t. the initial condition
    # Construction of the << sensitivity >> matrix w.r.t. the initial condition
    Stck	 	:= linalg[augment](Stck,
    evalm(dSysOndV&*linalg[submatrix](   VarSensitivityMatrix,
					    1..NbVariables,
			    nops(parameters)+1..NbVariables+nops(parameters)
						))) :
    SLPArguments 	:= [op(SLPArguments),VarSensitivityMatrix]	:
    #--------------------------------------------------------------------------
    # This is not good. In fact, in the construction of linear variational
    # system, we use stck; this time we don't care about it.
    else Stck	 	:= evalm(dSysOndV&*linalg[submatrix](
					    VarSensitivityMatrix,
					    1..NbVariables,
			    nops(parameters)+1..NbVariables+nops(parameters)
						)) :
    SLPArguments 	:= [op(SLPArguments),VarSensitivityMatrix]	:
    fi:
    #--------------------------------------------------------------------------
    # Output: the system and the argument necessary to evaluate it.
    [Stck,SLPArguments]

end: # ConstructOutput
#______________________________________________________________________________


#------------------------------------------------------------------------------
### Name                  : System2SLP
### Input                 : System	    :: array,
###			    Specializations :: list of inputs and parameters
###			    MyPrime	    :: integer computation are done
###							modulo MyPrime.
### Output                : a Maple procedure coding the TangentialSystem
### Description           : 
### References            : 
### Date                  : Wed Aug 16 14:35:09 2000
### Last modified         : 
### Implementation        : 
### Used procedure 	  : use codegen[optimize], codegen[makeproc].
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 : The complexity of evaluation seems not to be 
###			    very important in the practical complexity.
###			    Thus, we do not try to get the smaller possible
###			    program. May be later this point have to be checked

System2SLP := proc (	System 		:: array,
			Specializations	:: list,
			ProcArguments	:: list,
			MyPrime		:: integer	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

	local Stck, output :

        #----------------------------------------------------------------------
	# We specialise inputs and parameters in the system.
          Stck := array(	op(2,op(1,System)),
				eval(op(3,op(1,System)),Specializations)  ):
        #----------------------------------------------------------------------
	# Construction of the SLP
	output := [codegen[optimize](Stck)] : unassign('Stck') 		      :
	output := map(proc(z) op(1,z) = modp(op(2,z),MyPrime) end,output)     :
	output := map(	proc(z) op(1,z) = 'modp'('series'(op(2,z),t),MyPrime) 
			end , output	):
	output := [op(output),Stck] 					      :
        #----------------------------------------------------------------------
	# Construction of a procedure and output
	codegen[makeproc](output,ProcArguments) 

end: # System2SLP
#______________________________________________________________________________



		#-----------------------------------------------#	
		# 	    Power series expansion		#
		#-----------------------------------------------#	


#------------------------------------------------------------------------------
### Name                  : InverseMatrixSeries
### Input                 : Mat :: matrix at precision LocalOrder/2
###			  : LocalOrder :: integer, MatSize :: Range,
###			    MyPrime :: integer.
### Output                : 1/Mat inverse of Mat at precision LocalOrder 
###			    modulo MyOrder
### Description           : Use a Newton Operator 
###				InvAA = Id + O(t) 
###				-(InvAA-Id)^2 = O(t^2)
###				2InvAA - InvAAInvAA = Id + O(t^2)
###				(2InvA - InvAAInvA) A = Id + O(t^2)
###
### References            : 
###				@Article{BrentKung:1978,
###				  author =	 {R. P. Brent and H. T. Kung},
###				  title =	 {Fast Algorithms for 
###						  Manipulating Formal Power
###				                  Series},
###				  journal =	 {Journal of the Association 
###						 for Computing Machinery},
###				  year =	 1978,
###				  volume =	 25,
###				  number =	 4,
###				  pages =	 {581-595},
###				  month =	 oct
###			       }
### Date                  : Thu Sep 21 14:46:51 2000
### Last modified         : 
### Implementation        : recursive
### Used procedure	  : modp, add series
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 :  Computation are done modulo myprime.
###
###               	     When we work with series in Maple, we convert 
###			     them in a polynomial in order to perform a
###			     truncation for series. (see
###		          series(series(sin(t),t,2)*series(cos(t),t,10),t,10));
###			     This causes a loss in efficiency.
###
###			     The option remember is not used because this 
###			     procedure is called with the same Mat but with
###			     different LocalOrder.

InverseMatrixSeries := proc(  	Mat        :: matrix, 
				LocalOrder :: integer, 
				MatSize	   :: range,
				MyPrime	   :: integer		)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

           local Stck,Stck2,i,j,k,NewOrder :

	   if  LocalOrder = 1 then linalg[inverse](map(coeff,Mat,t,0))
	       else  
		NewOrder := 2^ceil(evalf(log[2](LocalOrder)-1))          :
		Stck 	   := map(series,Mat,t,NewOrder)		 :
	  	Stck 	   := InverseMatrixSeries(   Stck, NewOrder,
						     MatSize, MyPrime  ) :
		Stck 	   := map(series,Stck,t,LocalOrder)		 :
		#--------------------------------------------------------------
	        Stck2 := array(MatSize,MatSize,[seq([seq(modp(
		   	series(-add(Stck[i,k]*Mat[k,j],k=MatSize),t,LocalOrder)
			 ,MyPrime),j=MatSize)],i=MatSize)]) :
		#--------------------------------------------------------------
		# 2*identity - Stck2
		    for i from 1 to op(2,MatSize) 
			do Stck2[i,i] := series(2+Stck2[i,i],t,LocalOrder) od:
		#--------------------------------------------------------------
		# (2*identity- Stck*mat)*Stck
		map2(subs,O(1)=0,array(MatSize,MatSize,[seq([seq(modp(
           	series(add(Stck2[i,k]*Stck[k,j],k=MatSize),t,LocalOrder)
		,MyPrime),j=MatSize)],i=MatSize)])) 
	   fi

end : #InverseMatrixSeries
#------------------------------------------------------------------------------

#-----------------------------------------------------------------------------#
### Name                  : ExpMatrixSeries
### Input                 : Mat :: matrix at precision LocalOrder/2
###			  : LocalOrder :: integer, MatSize :: Range,
###			    MyPrime :: integer.
### Output                : exp(Series) at precision myorder
### Description           : use a Newton operator  
###  		            ExpMat*(1-int(diff(ExpMat,t)*InvExpMat,t) - Mat)
### References            : see InverseMatrixSeries
### Date                  : Fri Apr 21 16:19:01 2000
### Last modified         : 
### Implementation        : recursive
### Complexity		  : 
### Error Conditions      : 
### Environment Variables : Order
### Side Effects          : 
### Other                 : Computation are done modulo myprime.
###
###			    This procedure is made for Mat of valuation > 0.
###
###			    We compute the exponantial and the inverse in the
###			    same time.

ExpMatrixSeries := proc ( 	Mat 	   :: matrix, 
				LocalOrder :: integer,
				MatSize    :: range,
				MyPrime	   :: integer	  )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

	local Stck, Out1, Out2, i,j,k, NewOrder :

	   if  LocalOrder = 1 then [  	array(identity,MatSize,MatSize),
					array(identity,MatSize,MatSize)   ]

		else	
		NewOrder := 2^ceil(evalf(log[2](LocalOrder)-1)) :
		Stck 	 := map(series, Mat, t, LocalOrder)		:
	  	Stck     := ExpMatrixSeries(	Stck,    NewOrder,
						MatSize, MyPrime     ) 	:

		Stck[1] := map(series, Stck[1], t, LocalOrder) :
		Stck[2] := map(series, Stck[2], t, LocalOrder) :	
		#--------------------------------------------------------------

		#--------------------------------------------------------------
		# Inversion of Stck[1] using Stck[2].
		## - Stck[1] Stck[2]
	  	Out2 := array(MatSize, MatSize, [seq( [seq( modp(
		series(-add(Stck[1][i,k]*Stck[2][k,j],k=MatSize),t,NewOrder)
	   		,MyPrime), j=MatSize)], i=MatSize)]) :

		## 2*identity - Stck[1] Stck[2]  
		    for i from 1 to op(2,MatSize) 
			do Out2[i,i] := series(2+Out2[i,i],t,NewOrder) od:

		## Stck[2](2*identity - Stck[1] Stck[2])
		    Out2 := map2(subs,O(1)=0,array(MatSize,MatSize,
			[seq([seq(modp(
		series(add(Stck[2][i,k]*Out2[k,j],k=MatSize),t,NewOrder)
			,MyPrime),j=MatSize)],i=MatSize)])) :
		#--------------------------------------------------------------

		#--------------------------------------------------------------
		# Computation of the exponantial 
		## Computation of 
		    Out1 := array(MatSize,MatSize,[seq([seq(modp(
	    series(add(diff(Mat[i,k],t)*Stck[1][k,j],k=MatSize),t,LocalOrder)
			    ,MyPrime),j=MatSize)],i=MatSize)]) :

		    Out1 := array(MatSize,MatSize,[seq([seq(modp(
			    series(Out1[i,j]-diff(Stck[1][i,j],t),t,LocalOrder)
			    ,MyPrime),j=MatSize)],i=MatSize)]) :

		    Out1 := array(MatSize,MatSize,[seq([seq(modp(
				    series(int(series(
					    add(Out2[i,k]*Out1[k,j],k=MatSize)
				    ,t,LocalOrder),t),t,LocalOrder)
			    ,MyPrime),j=MatSize)],i=MatSize)]) :
		##  
		   for i from 1 to op(2,MatSize) 
		      do Out1[i,i] := series(1+Out1[i,i],t,LocalOrder) od:

		## Stck[1](identity + Mat - log(Stck[1]) )    
		   Out1 := map2(subs,O(1)=0,array(MatSize,MatSize,
			[seq([seq(modp(
		    series(add(Stck[1][i,k]*Out1[k,j],k=MatSize),t,LocalOrder)
			    ,MyPrime),j=MatSize)],i=MatSize)])) :
		#--------------------------------------------------------------

		[ Out1, Out2 ]
	   fi
end: # ExpMatrixSeries
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Name                  : HomogeneousResolution
### Input                 : LogDer :: matrix, the logarithmic derivative 
###		            associated to the system of ODE.
###			    LocalOrder :: integer, MyPrime :: integer.	
### Output                : a matrix M solution of homogenous system 
###				dM/dt + LogDer&*M = 0
###			    at precision LocalOrder
### Description           : 
### References            : 
### Date                  : Fri Sep 15 10:39:56 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 : The initial condition of the homogeneous solution
###			    is zero.

HomogeneousResolution := proc (   LogDer     :: matrix, 
				  LocalOrder :: integer,
				  MyPrime    :: integer   )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local Stck,MatSize :

    Stck   	:= map(series,map(-int,LogDer,t),t):    

    MatSize	:= 1..linalg[coldim](Stck) :

    ExpMatrixSeries(	Stck,	 LocalOrder,
			MatSize, MyPrime  	)[1] :

end: # HomogeneousResolution
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Name                  : PreprocessingForVariationOfConstants
### Input                 : HomSol :: table, InvSeparant :: matrix, 
###			    LocalOrder :: integer, MyPrime :: integer
### Output                : A matrix used in variation of constants
### Description           : 
### References            : 
### Date                  : Fri Sep 15 11:27:15 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 : This procedure allows to avoid the repetition of
###			    a computation in VariationOfConstants.

PreprocessingForVariationOfConstants := proc (  HomSol      :: table,
					        InvSeparant :: matrix,
						LocalOrder  :: integer,
						MyPrime	    :: integer  )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local Stck,i,j,k,matsize :

    matsize := 1..linalg[rowdim](InvSeparant)				:
    Stck    := InverseMatrixSeries(HomSol,LocalOrder,matsize,MyPrime)   :

	array(matsize,matsize,[seq([seq(modp(
		series(add(Stck[i,k]*InvSeparant[k,j],k=matsize),t,LocalOrder)
	,MyPrime),j=matsize)],i=matsize)])

end: # PreprocessingForVariationOfConstants
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
### Name                  : VariationOfConstants
### Input                 : HomSol :: matrix, InvHomSolInvA :: matrix
###			    SecondMember :: matrix, MyPrime :: integer.	
### Output                : - HomSol * int(InvHomSolInvA &* SecondMember,t)
### Description           : 
### References            : 
### Date                  : Fri Sep 15 12:13:34 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 : 

VariationOfConstants := proc (	HomSol        :: matrix, 
				InvHomSolInvA :: matrix,
				SecondMember  :: matrix,
				MyPrime	      :: integer	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :
				
    local Stck, colsize, rowsize :

    colsize := 1..linalg[coldim](SecondMember) :
    rowsize := 1..linalg[rowdim](SecondMember) :

    Stck := array(rowsize,colsize,[seq([seq(modp(
	-int(series(add(InvHomSolInvA[i,k]*SecondMember[k,j],k=rowsize),t),t)
	,MyPrime),j=colsize)],i=rowsize)]) :

    map2(subs,O(1)=0,array(rowsize,colsize,[seq([seq(modp(
    series(add(HomSol[i,k]*Stck[k,j],k=rowsize),t)
    ,MyPrime),j=colsize)],i=rowsize)]))
    
end: # VariationOfConstants
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
### Name                  : SecondMembersEvaluation
### Input                 : Sol :: matrix, SolParameters :: matrix, 
###			    SolInitialCondition :: matrix, SLP :: procedure.
### Output                : 
### Description           : This procedure evaluate the power series Sol(...)
###			    on SLP and return a list such that:
###				- list[1] is the inverse of the separent;
###				- list[2] is the logarithmic derivative;
###				- list[3] is the system;
###				- list[4] is the second member of linear 
###				variational system.
### References            : 
### Date                  : Wed Sep 13 12:17:10 2000
### Last modified         : 
### Implementation        : 
### Error Conditions      : 
### Environment Variables : The SLP coding the linear variational system
### Side Effects          : 
### Other                 : 

SecondMembersEvaluation := proc ( Sol 		      :: matrix, 
				  SensitivityMatrix   :: matrix,	
				  SLP   	      :: procedure 	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

	local Out1, Out2, Out3, Out4, Stck, NbVariables,NbSensitivityVar,i,j : 

	#----------------------------------------------------------------------
	# Some shortcuts
	  NbVariables  := linalg[rowdim](Sol) 		:
	  NbSensitivityVar := linalg[coldim](SensitivityMatrix) :
	#----------------------------------------------------------------------
	# Preparation of the SLP's inputs
	  Stck	:= convert(convert(Sol,vector),list) 		    :  
	  Stck 	:= op(map(diff,Stck,t))		, 	    op(Stck),
		   map(diff,SensitivityMatrix,t), SensitivityMatrix :
	#----------------------------------------------------------------------
	# Evaluation of the SLP coding the linear variational system
	  Stck	:= SLP(Stck)	:
	#----------------------------------------------------------------------
	# Split the output of the SLP coding the linear variational system

	Out1	:= linalg[matrix](NbVariables,NbVariables,
	    [seq(seq(Stck[i,j],j=2..NbVariables+1),i=1..NbVariables)] ): 

	Out2	 := linalg[matrix](NbVariables,NbVariables,
      [seq(seq(Stck[i,j],j=2+NbVariables..2*NbVariables+1),i=1..NbVariables)]):

       Out3 := linalg[matrix](NbVariables,1,[seq(Stck[i,1],i=1..NbVariables)]):
			
	Out4 := linalg[matrix](NbVariables,NbSensitivityVar,
	[seq(seq(Stck[i,j],j=2+2*NbVariables..1+2*NbVariables+NbSensitivityVar)
		,i=1..NbVariables)]) :
	#----------------------------------------------------------------------
 	[Out1,Out2,Out3,Out4]

end : # SecondMembersEvaluation
#------------------------------------------------------------------------------



		#-----------------------------------------------#	
		# 		Postprocessing			#
		#-----------------------------------------------#	



#------------------------------------------------------------------------------
### Name                  : ObservabilityAnalysis
### Input                 : JacobianMatrix :: matrix, Variables :: list, 
###			    Parameters :: list, MyPrime :: integer 
### Output                : 
### Description           : this procedure analyses the << sensitivity >>
###			    matrix of the output w.r.t.\ the initial condition
###			    in order to determine the transcendance degree
###			    of the field extension k(Y,Theta) -> k(Y,Theta,X).
###
###			    If this degree is zero than all the variables are
###			    observable and the procedure return a list composed
###			    of zero, of an empty list and the list of 
###                         observable<< variables >>.
###
###			    In the other case, this procedure computes the
###			    transcendence degree and return a list composed of
###				- l[1] the transcendance degree;
###				- l[2] the non observable << variables >>;
###				- l[3] the observable << variables >>.
###				<< variables >>:= variables + parameters.
###
###			    The computations are done modulo MyPrime.
###
### 		The matrix JacobianMatrix is encoded as below:
###
###				         Parameters  Variables
###		order 0 relations    | 				|
###		order 1 relations    |				|
###		     ....	     |				|
###		order n+l relations  |				|
###
### References            : 
### Date                  : Wed Sep 13 12:17:10 2000
### Last modified         : 
### Implementation        : 
### Used Procedure	  : Gausselim, linalg[delcols]
### Error Conditions      : 
### Environment Variables : 
### Side Effects          : 
### Other                 : 


AlgebraicDependence := proc( JacobianMatrix :: matrix,
			     Element	    :: list(name),
			     MyPrime	    :: integer,
			     TransDegree    :: integer	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local i, NonObservable, Rank :
    
    NonObservable := Element : 
    for i from 1 to nops(Element) do
       	userinfo(3,observabilityTest,`Checking`,Element[i]):
	Gausselim(linalg[delcols](JacobianMatrix,i..i),'Rank') mod MyPrime :

	if  evalb(nops(Element)-1-TransDegree = Rank) then 
	    NonObservable := remove(has,NonObservable,Element[i])
	fi
    od:

    NonObservable :
			         	     
end : # AlgebraicDependence


AlgebraicIndependence := proc( 	JacobianMatrix 	:: matrix,
				Element	    	:: list(name),
				MyPrime	    	:: integer,
				TransDegree    	:: integer	)

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local i, NonObservable, Rank :
    
    NonObservable := [] : 
    for i from 1 to nops(Element) do
       	userinfo(3,observabilityTest,`Checking`,Element[i]):
	Gausselim(linalg[delcols](JacobianMatrix,i..i),'Rank') mod MyPrime :

	if  evalb(nops(Element)-TransDegree = Rank) then 
	    NonObservable := [op(NonObservable),Element[i]] 
	fi
    od:

    NonObservable :
			         	     
end : # AlgebraicIndependence

SymetriesTower:=proc (	 JacMat 	:: matrix,
			 Variables      :: list, 
			 TransDegree	:: nonnegint,
			 MyPrime        :: nonnegint  )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

local i, Stck, output :

    if TransDegree = 1
     then   AlgebraicIndependence(JacMat,Variables,MyPrime,1):
     else   Stck := nops(Variables)-TransDegree :
	    for i from 1 to nops(Variables) 
	    while Stck=nops(Variables)-1-TransDegree do
		Gausselim(linalg[delcols](JacMat,i..i),'Stck') mod MyPrime :
	    od :
	    
	    SymetriesTower(	linalg[delcols](JacMat,i..i)	,
				subsop(i=NULL,Variables)	,
				TransDegree - 1			,
				MyPrime				) :
    fi 
    
end: # SymetriesTower

ObservabilityAnalysis := proc (	 JacobianMatrix :: matrix,
				 Variables      :: list, 
				 Parameters     :: list,		
				 MyPrime        :: integer  )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. Used by  observabilityTest" :

    local     	i, j, Rank, Rank2, TransDegree, MatrixWithoutACol, Mat, 
		NonObservable, Stck, LOrder,JacMat :

    Stck := [op(Parameters),op(Variables)] :
    #--------------------------------------------------------------------------
    JacMat := Gausselim(JacobianMatrix,'Rank') mod MyPrime :
    
    TransDegree := nops(Variables) + nops(Parameters) - Rank :
    
    userinfo(	1,observabilityTest,
	     	`The transcendence degree of k(U,Y) --> k(U,Y,X,Theta) is`,
    		TransDegree,cat(currentdir(),"/release")					    ) :
    #--------------------------------------------------------------------------
    if evalb(TransDegree = 0) then [TransDegree,[],Stck,0] 
	else
    NonObservable := AlgebraicIndependence(JacMat,Stck,MyPrime,TransDegree):
						
    userinfo(1,observabilityTest,NonObservable,`are not observable.\n`):
    #-------------------------------------------------------------------------
    Mat := [] : j := linalg[rowdim](JacobianMatrix) :
    for i from 1 to nops(Stck) do
        if  has(NonObservable,Stck[i]) then  Mat := [op(Mat),i]	fi
    od;
    JacMat := linalg[submatrix](JacobianMatrix,1..j,Mat):
    #--------------------------------------------------------------------------
    i := true : LOrder := 1 : j := linalg[coldim](JacMat) :
    Gausselim(JacMat,'Rank') mod MyPrime  :
    while i do
     Gausselim(	linalg[submatrix](JacMat,1..LOrder,1..j),'Rank2') 
		mod MyPrime  :
     if evalb(Rank2 = Rank) 
	then i := false  
	else LOrder := LOrder + 1
     fi  
    od:
    userinfo(3,observabilityTest,
	`least number of equations needed to determine the group of symetries`,
		LOrder,cat(currentdir(),"/release")):
    #--------------------------------------------------------------------------
#    if evalb(TransDegree = 1) then  
    	 [	TransDegree				,
		NonObservable			       	,
		remove(has,Stck,NonObservable)     	,
		LOrder					]
#    else
    #--------------------------------------------------------------------------
       	
#    [	TransDegree				    		,
#	SymetriesTower(JacMat,NonObservable,TransDegree,MyPrime),
#	remove(has,Stck,NonObservable)		    		,
#	LOrder					     			]
#	    fi :
    fi ;
end: # ObservabilityAnalysis
#------------------------------------------------------------------------------


		#-----------------------------------------------#	
		# 		Main procedure			#
		#-----------------------------------------------#	


#-----------------------------------------------------------------------------#
### Name                  : observabilityTest
###
###			    Remember that we study a system
###			       	diff(Theta,t) =  0	
###				  diff(X,t)   =  F(X,Theta,U),
###				       Y      =  G(X,Theta,U).
###
### Input                 : VectorField  :: list  F(X,Theta,U)		   ;
###			    Variables    :: list  the state variables X    ;
### 		WARNING	    diff(Variables[i],t) = VectorField[i].
###			    OutputSystem :: list  the outputs G(X,Theta,U) ;
###			    Parameters   :: list  the parameters Theta	   ;
###			    Inputs 	 :: list  the intputs U            .
###				
### Output                : see ObservabilityAnalysis
### Description           : 
###
###		          For the moment, the prime number is nextprime(10^10).
###
### References            : 
### Date                  :
### Last modified         : Wed Jul 19 13:00:08 CEST 2006
### Implementation        : 
### Used Procedure	  : nextprime, rand, randomize, ranpoly.
###			    ConstructLinearVariationalSystem,
###			    ConstructOutput, System2SLP,
###			    HomogeneousResolution, VariationOfConstants,
###			    PreprocessingForVariationOfConstants.
### Error Conditions      : error if nops(VectorsField) <> nops(Variables) 
###			    error if assigned(t)	
### Environment Variables :
### Side Effects          : 
### Complexity		  :
### Other                 : infolevel is checked at the begining of the 
###			    procedure. Hence, it is used as a global variable.
###			    Order is changed (it is an environment variable).
###
### To Do		  : the user must have the possibility to choose his
###			    own prime number.
###			    remove the global variables t.
	
observabilityTest := proc (	VectorField  :: list,
				Variables    :: list,
				OutputSystem :: list,
				Parameters   :: list,
				Inputs	     :: list	    )

description "Copyleft 2000 by Alexandre.Sedoglavic@lifl.fr. See ? observabilityTest" :

	global infolevel,t :

	local NbOutputs, NbVariables, NbParameters,
	      MyPrime, Hasard,
	      SpecializedVariables, Specialization,
	      i,Stck,T,
	      MainSLP, OutputSLP,
	      Sol, SensitivityMatrix,
	      SLPResult,
	      Bound, OneMoreLoop,
	      ComputationTime, TotalTime,
	      HomSol,InvHomSolInvA,
	      OldOrder, NewOrder, 
	      JacobianMatrix, ObservabilityResult	  :

	#----------------------------------------------------------------------
	# Check if all the quantities are defined
	Stck := indets(VectorField) union indets(OutputSystem) :
	Stck := Stck minus {op(Parameters),op(Variables),op(Inputs)};
	if Stck<>{} 
	    then  ERROR(`The quantities `,Stck,`are not defined`)
	fi :
			
	Stck := {op(Parameters),op(Variables),op(Inputs)} :
	Stck := Stck minus indets(VectorField) minus indets(OutputSystem) :
	if Stck<>{} 
	    then userinfo(1,observabilityTest,`The quantities `,Stck,
			    `are defined but not used`)
	fi :
	#----------------------------------------------------------------------
	# avoid maple garbage collecting message
	kernelopts(printbytes=false):

	    readlib(BaurStrassen) : readlib(dagnormal) :
	    
	#----------------------------------------------------------------------
	#Check the used infolevel.
	# if not(member([observabilityTest],[indices(infolevel)]))
	#	then infolevel[observabilityTest] := 1 
	# fi :
	 
	userinfo(1,observabilityTest,`Modular computation with version 0.0\n`):
	
	#----------------------------------------------------------------------
	# Some little test
	if nops(VectorField) <> nops(Variables) 
		then ERROR(`Variables must be ordered as vector field`) fi :
	if assigned(t) 
		then ERROR(`Sorry, for the moment t must be unassigned`) fi :
	#----------------------------------------------------------------------

       	#----------------------------------------------------------------------
	userinfo(1,observabilityTest,`Some informations about the system`)    :
	userinfo(1,observabilityTest,`Nb Inputs       U  :`,nops(Inputs))     :
	NbOutputs 	:= nops(OutputSystem)			              :
	userinfo(1,observabilityTest,`Nb Outputs      Y  :`,NbOutputs)        :
	NbVariables	:= nops(Variables) 				      :
	userinfo(1,observabilityTest,`Nb Variables    X  :`,NbVariables)      :
	NbParameters 	:= nops(Parameters)			      	      :
	userinfo(1,observabilityTest,`Nb Parameters Theta:`,NbParameters,cat(currentdir(),"/release")):
       	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# The probabilistic aspect of the algorithm
	# Choice of a prime number
	  MyPrime := nextprime(10^13)	:
	# Choice of the range where the specialization are randomly choosen
	# Set the seed for the random number generator
	  readlib(randomize)()		:
	  Hasard  := rand(0..MyPrime)	:
	#----------------------------------------------------------------------
	  userinfo(1,observabilityTest,
	    `The computation are done modulo the following prime number`,
	    MyPrime,cat(currentdir(),"/release")):
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# We choosen randomly some specialization according to the above choice
	  userinfo(3,observabilityTest,`Specialization are choosen`)	:


	  SpecializedVariables 	:= [seq( Hasard(), i=Variables )]	:    

	  Specialization := [ seq(i = Hasard(), i=Parameters), 
	  		      seq(i = randpoly(t,coeffs=rand(0..MyPrime),dense,
				             degree=1+NbParameters+NbVariables)
				, i = Inputs )		] : 
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# We construct the SLP used in the integration.
	  userinfo(1,observabilityTest,`System treatment`) : T := time():
	  userinfo(3,observabilityTest,`Construction of the main SLP`) :
	  Stck := ConstructLinearVariationalSystem( VectorField, Variables,
     						    Parameters	             ):
	  MainSLP := System2SLP(Stck[1],[op(Specialization)],Stck[2],MyPrime):
	#----------------------------------------------------------------------
	 userinfo(1,observabilityTest,`End of system treatment`,time()-T,cat(currentdir(),"/release")):
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# Initial conditions
	  Sol		:= linalg[matrix](NbVariables,1,SpecializedVariables):
  	  SensitivityMatrix:=linalg[augment](
				linalg[matrix](NbVariables,NbParameters,0),
				array(identity,1..NbVariables,1..NbVariables)):
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# Evaluation of the second member
	  SLPResult := SecondMembersEvaluation(Sol,SensitivityMatrix,MainSLP) :
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	# Loops counter
	  Order 	:= 2 				:
	  Bound 	:= NbVariables + NbParameters 	:
	  OneMoreLoop	:= true 			:
	#----------------------------------------------------------------------
	# Some landmark 
	  ComputationTime := time() 	:
	  TotalTime 	  := time() 	:
	#----------------------------------------------------------------------


         userinfo(1,observabilityTest,"Begining of integration process"):
         userinfo(1,observabilityTest,
		"For time estimation of this process, set current infolevel to 2"):
	userinfo(2,observabilityTest,
	"Convergence is quadratic. Thus, each step takes two time more then preceding ones"
	):
	#----------------------------------------------------------------------
	# Begining of the main loop
         userinfo(2,observabilityTest,`Power series expansion at order`,Bound):
	 while OneMoreLoop do
		userinfo(2,observabilityTest,`->  Order`,Order) :
		if Order = Bound then OneMoreLoop := false fi 		:
		#--------------------------------------------------------------
		HomSol := HomogeneousResolution(SLPResult[2],Order,MyPrime) :
		InvHomSolInvA := PreprocessingForVariationOfConstants(
					HomSol,SLPResult[1], Order, MyPrime ):
		Stck  := VariationOfConstants(	HomSol,InvHomSolInvA,
						SLPResult[3], MyPrime	):
		Sol   := evalm(Sol + Stck) :
		#--------------------------------------------------------------
		OldOrder := Order 		:
		if 2*Order < Bound then Order:=2*Order else Order:=Bound fi : 
		NewOrder := Order 		:
		#--------------------------------------------------------------
		SLPResult := SecondMembersEvaluation(	Sol,SensitivityMatrix,
							MainSLP	 ):
		Order  := OldOrder :
		#--------------------------------------------------------------
		Stck := VariationOfConstants(	HomSol,InvHomSolInvA,
						SLPResult[4], MyPrime	):
		SensitivityMatrix := evalm(SensitivityMatrix + Stck) :
		Order  := NewOrder :
		#--------------------------------------------------------------
		userinfo(2,observabilityTest,
		`->  Computation time`,time() - ComputationTime)
	  od :
  	  userinfo(1,observabilityTest,`End of integration \n`) :
	# End of the main loop
	#----------------------------------------------------------------------

	#----------------------------------------------------------------------
	  userinfo(3,Observabilitytest,`Construction of the output SLP`) :
	  Stck      := ConstructOutput(OutputSystem,Variables,Parameters):
	  OutputSLP :=System2SLP(Stck[1],[op(Specialization)],Stck[2],MyPrime):
	#----------------------------------------------------------------------
	  userinfo(1,observabilityTest,`Evaluation of output system`) :
	  Stck := convert ( convert(Sol,vector),list ) : 
	  Stck := OutputSLP ( op(Stck),SensitivityMatrix ) :
	JacobianMatrix := linalg[matrix](NbOutputs, NbParameters+NbVariables,
			 [seq( seq(Stck[i,j],j=1..NbParameters+NbVariables)
			       ,i=1..NbOutputs )           		    ]):
	JacobianMatrix := linalg[stackmatrix](
			seq(map(coeff,JacobianMatrix,t,i),i=0..OldOrder-1)) :

	userinfo(1,observabilityTest,`End of evaluation of output system \n`):
	#----------------------------------------------------------------------
	ObservabilityResult := ObservabilityAnalysis( JacobianMatrix,Variables,
						      Parameters,MyPrime    ) :
	#----------------------------------------------------------------------
	# Take in account the number of output in the necessary 
	# order of derivation
	ObservabilityResult := subsop(	
	    4=ceil(op(4,ObservabilityResult)/nops(OutputSystem)),
		ObservabilityResult					) :
	#----------------------------------------------------------------------
	TotalTime := time() - TotalTime :
	userinfo(1,observabilityTest,`Total used time`,TotalTime,cat(currentdir(),"/release")):
	#----------------------------------------------------------------------

	# set maple garbage collecting message
	  kernelopts(printbytes=true):
	#----------------------------------------------------------------------

	ObservabilityResult

end : #ObservabilityTest
#-----------------------------------------------------------------------------#

savelib(
	`ConstructLinearVariationalSystem`,
	`ConstructOutput`,
	`System2SLP`,
	`InverseMatrixSeries`,
	`ExpMatrixSeries`,
	`HomogeneousResolution`,
	`PreprocessingForVariationOfConstants`,
	`VariationOfConstants`,
	`SecondMembersEvaluation`,
	`AlgebraicDependence`,
	`AlgebraicIndependence`,
	`SymetriesTower`,
	`ObservabilityAnalysis`,
	`observabilityTest`
	) ;
quit :

