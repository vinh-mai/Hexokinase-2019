#======================#
# SENSITIVITY ANALYSIS #
#======================#

from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
from scipy.integrate import odeint
from scipy import integrate

# SOLUTION TO THE MODEL
#
def SOL(p, initial_cond, t0, t_end, stpz):
	#
    k0, k1, k_1, k2, k_2, k3, k_3, k4, k_4, k5, k_5, \
	k6, k_6, k7, k_7, k8, k_8, k9, k_9 = p
 	#
    t = np.arange(t0, t_end, stpz)
    #
    def MODEL(y,t):
    	#
		#
    	# Define y
    
        _E_, _0_, _1_, _2_, _3_, _4_,\
			\
    	_0E_,_1E_, _2E_, _3E_,\
			\
    	_E0_, _E1_, _E2_, _E3_,\
			\
    	_01E_, _02E_, _03E_,\
			\
    	_0E0_, _0E1_, _0E2_, _0E3_,\
			\
    	_1E0_, _1E1_, _1E2_, _1E3_,\
			\
    	_2E0_, _2E1_, _2E3_,\
			\
    	_3E0_, _3E1_, _3E2_, _3E3_,\
			\
    	_E01_, _E02_, _E03_,\
			\
    	_01E0_, _01E1_, _01E2_, _01E3_,\
			\
    	_02E0_, _02E1_, _02E3_,\
			\
    	_03E0_, _03E1_, _03E2_, _03E3_,\
			\
    	_0E01_, _1E01_, _2E01_, _3E01_,\
			\
    	_0E02_, _1E02_, _3E02_,\
			\
    	_0E03_, _1E03_, _2E03_, _3E03_,\
			\
    	_01E01_, _01E02_, _01E03_, _02E01_,\
			\
    	_02E03_, _03E01_, _03E02_, _03E03_ = y 
		#
		#
    	# _E_: Hexokinase 1; _0_: Glucose
		# _1_: ATP; _2_: G6P; _3_: Pi; _4_: ADP
		# _xEy_: x, y substances bound at N, C domains, respectively.
		# Define dydt
        dydt=[]
		#
    	#1 Eq. for enyme _E_
        dydt.append(k0*_E01_ + k_1*_0E_ + k_2*_E0_ + k_3*_1E_ \
        
		+ k_4*_E1_ + k_5*_2E_ + k_6*_E2_ + k_7*_3E_	+ k_8*_E3_ \
            
		- _E_*((k1 + k2)*_0_ + (k3 + k4)*_1_ + (k5 + k6)*_2_ \
		
			+ (k7 + k8)*_3_))

	    #2 Eq. for _0_
        dydt.append(k_1*(_0E_ + _01E_ + _02E_ + _03E_ + _0E0_ \
    
    		+ _0E1_ + _0E2_ + _0E3_ + _01E0_ + _01E1_ + _01E2_ \
    		
    		+ _01E3_ + _02E0_ + _02E1_ + _02E3_	+ _03E0_ \
    		
    		+ _03E1_ + _03E2_ + _03E3_ + _0E01_ + _0E02_ \
    		
    		+ _0E03_ + _01E01_ + _01E02_ + _01E03_ + _02E01_ \
    		
    		+ _02E03_ + _03E01_	+ _03E02_ + _03E03_ )\
			
		+ k_2*(_E0_ + _0E0_ + _1E0_ + _2E0_ + _3E0_ + _E01_ \
		
			+ _E02_ + _E03_ + _01E0_ + _02E0_ + _03E0_ + _0E01_ \
			
			+ _1E01_ + _2E01_ + _3E01_ +_0E02_ + _1E02_ + _3E02_ \
			
			+ _0E03_ + _1E03_ + _2E03_ + _3E03_	+ _01E01_ \
			
			+ _01E02_ + _01E03_ + _02E01_ + _02E03_ + _03E01_ \

			+ _03E02_ + _03E03_ )\
		
		- _0_*((k1 + k2)*(_E_ + _1E_ + _E1_ + _2E_ + _E2_ + _3E_ \
		
				+ _E3_ + _1E1_ + _1E2_ + _1E3_ + _2E1_ + _2E3_ \
				
				+ _3E1_ + _3E2_ + _3E3_)\

			+ k1*(_E0_ + _1E0_ + _2E0_ + _3E0_ + _E01_ + _E02_ \
			
				+ _E03_ + _1E01_ + _2E01_ + _3E01_ + _1E02_ \
				
				+ _3E02_ + _1E03_ + _2E03_ + _3E03_) \

			+ k2*(_0E_ + _01E_ + _02E_ + _03E_ + _0E1_ + _0E2_ \
			
				+ _0E3_ + _01E1_ + _01E2_ + _01E3_ + _02E1_ \
				
				+ _02E3_ + _03E1_ + _03E2_ + _03E3_)))

	    #3 Eq. for _1_
        dydt.append(k_3*(_1E_ + _01E_ + _1E0_ + _1E1_ + _1E2_ + _1E3_ \
    
    		+ _01E0_ + _01E1_ + _01E2_ + _01E3_ + _1E01_ + _1E02_ \
    	
    		+ _1E03_ + _01E01_ + _01E02_ + _01E03_) \

		+ k_4*(_E1_ + _0E1_ + _1E1_ + _2E1_ + _3E1_ + _E01_ \
		
			+ _01E1_ + _02E1_ + _03E1_ + _0E01_ + _1E01_ \
		
			+ _2E01_ + _3E01_ + _01E01_ + _02E01_ + _03E01_) \

		- _1_*(k3*(_E_ + _0E_ + _E0_ + _E1_ + _E2_ + _E3_ + _0E0_ \
		
				+ _0E1_ + _0E2_ + _0E3_ + _E01_ + _E02_ + _E03_ \
			
				+ _0E01_ + _0E02_ + _0E03_) \

			+ k4*(_E_ + _0E_ + _1E_ + _3E_ + _E0_ + _01E_ + _03E_ \
			
				+ _0E0_ + _1E0_ + _3E0_ + _01E0_ + _03E0_)))		

    	#4 Eq. for _2_
        dydt.append(k0*(_E01_ + _0E01_ + _1E01_ + _2E01_ + _3E01_ \

    			+ _01E01_ + _02E01_ + _03E01_) \
    
			+ k_5*(_2E_ + _02E_ + _2E0_ + _2E1_ + _2E3_ + _02E0_ \

				+ _02E1_ + _02E3_ + _2E01_ + _2E03_ + _02E01_ \
				
				+ _02E03_) \

			+ k_6*(_E2_ + _0E2_ + _1E2_ + _E02_ + _01E2_ \
				
				+ _0E02_ + _1E02_ + _01E02_) \

			+ k_9*(_3E2_ + _03E2_ + _3E02_ +  _03E02_ ) \

		- _2_*(k5*(_E_ + _0E_ + _E0_ + _E1_ + _E3_ + _0E0_ \
		
				+ _0E1_ + _0E3_ + _E01_ + _E03_ + _0E01_ + _0E03_) \

			+ k6*(_E_ + _0E_ + _1E_ + _E0_ + _01E_ + _0E0_ \
			
				+ _1E0_ + _01E0_) \

			+ k9*(_3E_ + _03E_ + _3E0_ + _03E0_)))

    	#5 Eq. for _3_
        dydt.append(k_7*(_3E_ + _03E_ + _3E0_ + _3E1_ + _3E2_ + _3E3_ \

			+ _03E0_ + _03E1_ + _03E2_ + _03E3_ + _3E01_ \
			
			+ _3E02_ + _3E03_ + _03E01_ + _03E02_ + _03E03_) \

		+ k_8*(_E3_ + _0E3_ + _1E3_ + _2E3_ + _3E3_ + _E03_ \

			+ _01E3_ + _02E3_ + _03E3_ + _0E03_ + _1E03_ \
			
			+ _2E03_ + _3E03_ + _01E03_ + _02E03_ + _03E03_) \

		- _3_*(k7*(_E_ + _0E_ + _E0_ + _E1_ + _E2_ + _E3_ + _0E0_ \

				+ _0E1_ + _0E2_ + _0E3_ + _E01_ + _E02_ + _E03_ \

				+ _0E01_ + _0E02_ + _0E03_) \

			+ k8*(_E_ + _0E_ + _1E_ + _3E_ + _E0_ + _01E_ + _03E_ \
			
				+ _0E0_ + _1E0_ + _3E0_ + _01E0_ + _03E0_)))

    	#6 Eq. for _4_
        dydt.append(k0*(_E01_ + _0E01_ + _1E01_ + _2E01_\
                + _3E01_ + _01E01_ + _02E01_ + _03E01_))
    	#
		#
		#
    	#7 Eq. for _0E_
        dydt.append(k0*_0E01_ + k1*_0_*_E_ + k_2*_0E0_ + k_3*_01E_ \
    
    		+ k_4*_0E1_ + k_5*_02E_ + k_6*_0E2_ + k_7*_03E_ \
    		
    		+ k_8*_0E3_ \

		- _0E_*(k_1 + k2*_0_ + (k3 + k4)*_1_ + (k5 + k6)*_2_ \
		
			+ (k7 + k8)*_3_))

	    #8 Eq. for _1E_
        dydt.append(k0*_1E01_ + k_1*_01E_ + k_2*_1E0_ + k_4*_1E1_ \
    
    		+ k_6*_1E2_ + k_8*_1E3_ + _1_*k3*_E_ \

		- _1E_*(k_3 + (k1 + k2)*_0_ + k4*_1_ + k6*_2_ + k8*_3_))

		#9 Eq. for _2E_
        dydt.append(k0*_2E01_ + k_1*_02E_ + k_2*_2E0_ + k_4*_2E1_ \
    
    		+ k_8*_2E3_ + _2_*k5*_E_ \

		- _2E_*(k_5 + (k1 + k2)*_0_))

    	#10 Eq. for _3E_
        dydt.append(k0*_3E01_ + k_1*_03E_ + k_2*_3E0_ + k_4*_3E1_ \
    
    		+ k_9*_3E2_ + k_8*_3E3_ + _3_*k7*_E_ \

		- _3E_*(k_7 + (k1 + k2)*_0_ + k4*_1_ + k9*_2_ + k8*_3_))

	    #11 Eq. for _E0_
        dydt.append(k2*_E_*_0_ + k_1*_0E0_ + k_3*_1E0_ + k_5*_2E0_ \

			+ k_7*_3E0_	+ k_4*_E01_	+ k_6*_E02_  + k_8*_E03_ \

		- _E0_*(k_2 + k1*_0_ + (k3 + k4)*_1_ + _2_*(k5 + k6) \
		
			+ (k7 + k8)*_3_)) 

		#12 Eq. for _E1_
        dydt.append(k_1*_0E1_ + k_2*_E01_ + k_3*_1E1_ + k_5*_2E1_ \
    
    		+ k_7*_3E1_ + _1_*_E_*k4 \

		- _E1_*(k_4 + (k1 + k2)*_0_ + k3*_1_ + k5*_2_ + k7*_3_))

		#13 Eq. for _E2_
        dydt.append(k_1*_0E2_ + k_2*_E02_ + k_3*_1E2_ + k_7*_3E2_ \

			+ _2_*k6*_E_ \

		-_E2_*(k_6 + (k1 + k2)*_0_ + k3*_1_ + k7*_3_))

		#14 Eq. for _E3_
        dydt.append(k_1*_0E3_ + k_2*_E03_ + k_3*_1E3_ + k_5*_2E3_ \
    
    		+ k_7*_3E3_ + _3_*k8*_E_ \

		- _E3_*(k_8 + (k1 + k2)*_0_ + k3 *_1_ + k5*_2_ + k7*_3_))

		#
		#
		#15 Eq. for _01E_
        dydt.append(k0*_01E01_ + k_2*_01E0_ + k_4*_01E1_ \
        
        	+ k_6*_01E2_ + k_8*_01E3_ + _0_*_1E_*k1 + _1_*k3*_0E_ \

		- _01E_*(k_1 + k_3 + k2*_0_ + k4*_1_ + k6*_2_ + k8*_3_))

		#16 Eq. for _02E_
        dydt.append(k0*_02E01_ + k_2*_02E0_ + k_4*_02E1_ \
        
        	+ k_8*_02E3_ + k1*_0_*_2E_ + _2_*k5*_0E_ \

		- _02E_*(k_1 + k_5 + k2*_0_))

		#17 Eq. for _03E_
        dydt.append(k0*_03E01_ + k_2*_03E0_ + k_4*_03E1_ + k_9*_03E2_ \
    		
    		+ k_8*_03E3_ + _0_*_3E_*k1 + _3_*k7*_0E_ \

		- _03E_*(k_1 + k_7 + k2*_0_ + k4*_1_ + k9*_2_ + k8*_3_))
		
		#
    	#18 Eq. for _0E0_
        dydt.append(k_3*_01E0_ + k_4*_0E01_ + k_5*_02E0_ + k_6*_0E02_ \
    
    		+ k_7*_03E0_ + k_8*_0E03_ + _0_*(k1*_E0_ + k2*_0E_) \

		- _0E0_*(k_1 + k_2 + (k3 + k4)*_1_ + (k5 + k6)*_2_ \
		
			+ (k7 + k8)*_3_))

		#19 Eq. for _0E1_
        dydt.append(k_2*_0E01_ + k_3*_01E1_ + k_5*_02E1_ + k_7*_03E1_ \
    
	    	+ k1*_0_*_E1_ + _1_*k4*_0E_ \

		- _0E1_*(k_1 + k_4 + k2*_0_ + k3*_1_ + k5*_2_ + k7*_3_))

		#20 Eq. for _0E2_
        dydt.append(k_2*_0E02_ + k_3*_01E2_ + k_7*_03E2_ \
        
        	+ k1*_0_*_E2_ + _2_*k6*_0E_ \

		- _0E2_*(k_1 + k_6 + k2*_0_ + k3*_1_ + k7*_3_))

    	#21 Eq. for _0E3_
        dydt.append(k_2*_0E03_ + k_3*_01E3_ + k_5*_02E3_ \
        
        	+ k_7*_03E3_ + k1*_0_*_E3_ + _3_*k8*_0E_ \

		- _0E3_*(k_1 + k_8 + k2*_0_ + k3*_1_ + k5*_2_ + k7*_3_))

		#
		#22 Eq. for _1E0_
        dydt.append (k_1*_01E0_ + k_4*_1E01_ + k_6*_1E02_ \
    
    		+ k_8*_1E03_ + k2*_0_*_1E_ + _1_*k3*_E0_ \

		- _1E0_*(k_2 + k_3 + k1*_0_ + k4*_1_ + k6*_2_ + k8*_3_))

		#23 Eq. for _1E1_
        dydt.append(k_1*_01E1_ + k_2*_1E01_ + _1_*(k3*_E1_ + k4*_1E_) \

			- _1E1_*(k_3 + k_4 + (k1 + k2)*_0_))

		#24 Eq. for _1E2_
        dydt.append(k_1*_01E2_ + k_2*_1E02_ + _1_*k3*_E2_ \
        	
        	+ _2_*k6*_1E_ \

			- _1E2_*(k_3 + k_6 + (k1 + k2)*_0_))

		#25 Eq. for _1E3_
        dydt.append(k_1*_01E3_ + k_2*_1E03_ \
        
        	+ _1_*k3*_E3_ + _3_*k8*_1E_ \

			- _1E3_*(k_3 + k_8 + (k1 + k2)*_0_))

		#
		#26 Eq. for _2E0_
        dydt.append(k_1*_02E0_ + k_4*_2E01_ + k_8*_2E03_ \
    
    		+ k2*_0_*_2E_ + _2_*k5*_E0_ \
    		
	    	- _2E0_*(k_2 + k_5 + k1*_0_))

		#27 Eq. for _2E1_
        dydt.append(k_1*_02E1_ + k_2*_2E01_ + _2_*k5*_E1_ \

			- _2E1_*(k_4 + k_5 + (k1 + k2)*_0_))

    	#28 Eq. for _2E3_
        dydt.append(k_1*_02E3_ + k_2*_2E03_ + _2_*k5*_E3_ \

			- _2E3_*(k_5 + k_8 + (k1 + k2)*_0_))

		#
		#29 Eq. for _3E0_
        dydt.append(k_1*_03E0_ + k_4*_3E01_ + k_9*_3E02_ \
    
    		+ k_8*_3E03_ + _0_*k2*_3E_ + _3_*k7*_E0_ \

		- _3E0_*(k_2 + k_7 + k1*_0_ + k4*_1_ + k9*_2_ + k8*_3_))

		#30 Eq. for _3E1_
        dydt.append(k_1*_03E1_ + k_2*_3E01_ + _1_*k4*_3E_ \
        
        	+ _3_*k7*_E1_ \

			- _3E1_*(k_4 + k_7 + (k1 + k2)*_0_))

		#31 Eq. for _3E2_
        dydt.append(k_1*_03E2_ + k_2*_3E02_ + _2_*k9*_3E_ \
        
        	+ _3_*k7*_E2_ \

			- _3E2_*(k_7 + k_9 + (k1 + k2)*_0_))

		#32 Eq. for _3E3_
        dydt.append(k_1*_03E3_ + k_2*_3E03_ \
        
        	+ _3_*(k7*_E3_ + k8*_3E_) \

			- _3E3_*(k_7 + k_8 + (k1 + k2)*_0_))

		#
		#33 Eq. for _E01_
        dydt.append(k_1*_0E01_ + k_3*_1E01_ + k_5*_2E01_ \
    
    		+ k_7*_3E01_ + _0_*k2*_E1_ + _1_*k4*_E0_ \

		- _E01_*(k0 + k_2 + k_4 + k1*_0_ + k3*_1_ + k5*_2_ + k7*_3_))

		#34 Eq. for _E02_
        dydt.append(k_1*_0E02_ + k_3*_1E02_ + k_7*_3E02_ \
    
    		+ _0_*k2*_E2_ + _2_*k6*_E0_ \

			- _E02_*(k_2 + k_6 + k1*_0_ + k3*_1_ + k7*_3_))

		#35 Eq. for _E03_
        dydt.append(k_1*_0E03_ + k_3*_1E03_ + k_5*_2E03_ \
    
    		+ k_7*_3E03_ + _0_*k2*_E3_ + _3_*k8*_E0_ \

		- _E03_*(k_2 + k_8 + _0_*k1 + _1_*k3 + _2_*k5 + _3_*k7))

		#
		#
		#36 Eq. for _01E0_
        dydt.append(k_4*_01E01_ + k_6*_01E02_ + k_8*_01E03_ \
    
    		+ _0_*(k1*_1E0_ + k2*_01E_) + _1_*k3*_0E0_ \

		- _01E0_*(k_1 + k_2 + k_3 + _1_*k4 + _2_*k6 + _3_*k8))

		#37 Eq. for _01E1_
        dydt.append(k_2*_01E01_ + _0_*k1*_1E1_ \
    
    		+ _1_*(k3*_0E1_ + k4*_01E_) \

			- _01E1_*(k_1 + k_3 + k_4 + _0_*k2))

		#38 Eq. for _01E2_
        dydt.append(k_2*_01E02_ + _0_*k1*_1E2_ + _1_*k3*_0E2_ \

			+ _2_*k6*_01E_ \

			- _01E2_*(k_1 + k_3 + k_6 + _0_*k2))

		#39 Eq. for _01E3_
        dydt.append(k_2*_01E03_ + _0_*k1*_1E3_ \
    
    		+ _1_*k3*_0E3_ + _3_*k8*_01E_ \

			- _01E3_*(k_1 + k_3 + k_8 + _0_*k2 ))

		#
		#40 Eq. for _02E0_
        dydt.append(k_4*_02E01_ + k_8*_02E03_ \
    
    		+ _0_*(k1*_2E0_ + k2*_02E_) + _2_*k5*_0E0_ \

		- _02E0_*(k_1 + k_2 + k_5))

		#41 Eq. for _02E1_
        dydt.append(k_2*_02E01_ + _0_*k1*_2E1_ + _2_*k5*_0E1_ \

			- _02E1_*(k_1 + k_4 + k_5 + _0_*k2))

		#42 Eq. for _02E3_
        dydt.append(k_2*_02E03_ + _0_*k1*_2E3_ + _2_*k5*_0E3_ \

			- _02E3_*(k_1 + k_5 + k_8 + _0_*k2))

		#
		#43 Eq. for _03E0_
        dydt.append(k_4*_03E01_ + k_9*_03E02_ + k_8*_03E03_ \

    		+ _0_*(k1*_3E0_ + k2*_03E_) + _3_*k7*_0E0_ \

		- _03E0_*(k_1 + k_2 + k_7 + _1_*k4 + _2_*k9 + _3_*k8))

		#44 Eq. for _03E1_
        dydt.append(k_2*_03E01_ + _0_*k1*_3E1_ + _1_*k4*_03E_ \
    
    		+ _3_*k7*_0E1_ \

			- _03E1_*(k_1 + k_4 + k_7 + _0_*k2))

		#45 Eq. for _03E2_
        dydt.append(k_2*_03E02_ + _0_*k1*_3E2_ \
    
    		+ _2_*k9*_03E_ + _3_*k7*_0E2_ \

			- _03E2_*(k_1 + k_9 + k_7 + _0_*k2))

		#46 Eq. for _03E3_
        dydt.append(k_2*_03E03_ + _0_*k1*_3E3_ \
    
    		+ _3_*(k7*_0E3_ + k8*_03E_) \

			- _03E3_*(k_1 + k_7 + k_8 + _0_*k2))

		#
		#47 Eq. for _0E01_
        dydt.append(k_3*_01E01_ + k_5*_02E01_ + k_7*_03E01_ \
    
    		+ _0_*(k1*_E01_ + k2*_0E1_) + _1_*k4*_0E0_ \

		- _0E01_*(k0 + k_1 + k_2 + k_4 + _1_*k3 + _2_*k5 + _3_*k7))

		#48 Eq. for _1E01_
        dydt.append(k_1*_01E01_ + _0_*k2*_1E1_ \
    
    		+ _1_*(k3*_E01_ + k4*_1E0_) \

			- _1E01_*(k0 + k_2 + k_3 + k_4 + _0_*k1))

		#49 Eq. for _2E01_
        dydt.append(k_1*_02E01_ + _0_*k2*_2E1_ + _2_*k5*_E01_ \

			- _2E01_*(k0 + k_2 + k_4 + k_5 + _0_*k1))

		#50 Eq. for _3E01_
        dydt.append(k_1*_03E01_ + _0_*k2*_3E1_ + _1_*k4*_3E0_ \
    
    		+ _3_*k7*_E01_ \

			- _3E01_*(k0 + k_2 + k_4 + k_7 + _0_*k1))

		#
		#51 Eq. for _0E02_
        dydt.append(k_3*_01E02_ + k_7*_03E02_ 
    
    		+ _0_*(k1*_E02_ + k2*_0E2_) + _2_*k6*_0E0_ \

			- _0E02_*(k_1 + k_2 + k_6 + _1_*k3 + _3_*k7))

		#52 Eq. for _1E02_
        dydt.append(k_1*_01E02_ + _0_*k2*_1E2_ + _1_*k3*_E02_ \
    
    		+ _2_*k6*_1E0_ \

			- _1E02_*(k_2 + k_3 + k_6 + _0_*k1))

		#53 Eq. for _3E02_
        dydt.append(k_1*_03E02_ + _0_*k2*_3E2_ + _2_*k9*_3E0_ \
    
    		+ _3_*k7*_E02_ \

			- _3E02_*(k_2 + k_9 + k_7 + _0_*k1))

		#
		#54 Eq. for _0E03_
        dydt.append(k_3*_01E03_ + k_5*_02E03_ + k_7*_03E03_ \
    
    		+ _0_*(k1*_E03_ + k2*_0E3_) + _3_*k8*_0E0_ \

			- _0E03_*(k_1 + k_2 + k_8 + _1_*k3 + _2_*k5 + _3_*k7))

		#55 Eq. for _1E03_
        dydt.append(k_1*_01E03_ + _0_*k2*_1E3_ + _1_*k3*_E03_ \
    
    		+ _3_*k8*_1E0_ \

			- _1E03_*(k_2 + k_3 + k_8 + _0_*k1))

		#56 Eq. for _2E03_
        dydt.append(k_1*_02E03_ + _0_*k2*_2E3_ + _2_*k5*_E03_ \

			- _2E03_*(k_2 + k_5 + k_8 + _0_*k1))

		#57 Eq. for _3E03_
        dydt.append(k_1*_03E03_ + _0_*k2*_3E3_ \
    
    		+ _3_*(k7*_E03_ + k8*_3E0_) \

			- _3E03_*(k_2 + k_7 + k_8 + _0_*k1))

		#
		#
		#58 Eq. for _01E01_
        dydt.append(_0_*(k1*_1E01_ + k2*_01E1_) \
    
    		+ _1_*(k3*_0E01_ + k4*_01E0_) \

			- _01E01_*(k0 + k_1 + k_2 + k_3 + k_4))

		#59 Eq. for _01E02_
        dydt.append(_0_*(k1*_1E02_ + k2*_01E2_) \
    
    		+ _1_*k3*_0E02_ + _2_*k6*_01E0_ \

			- _01E02_*(k_1 + k_2 + k_3 + k_6))

		#60 Eq. for _01E03_
        dydt.append(_0_*(k1*_1E03_ + k2*_01E3_) \
    
    		+ _1_*k3*_0E03_  + _3_*k8*_01E0_ \

			- _01E03_*(k_1 + k_2 + k_3 + k_8))
		
		#
		#61 Eq. for _02E01_
        dydt.append(_0_*(k1*_2E01_ + k2*_02E1_) + _2_*k5*_0E01_ \

			- _02E01_*(k0 + k_1 + k_2 + k_4 + k_5))

		#62 Eq. for _02E03_
        dydt.append(_0_*(k1*_2E03_ + k2*_02E3_) + _2_*k5*_0E03_ \

			- _02E03_*(k_1 + k_2 + k_5 + k_8))

		#
		#63 Eq. for _03E01_
        dydt.append(_0_*(k1*_3E01_ + k2*_03E1_) \
    
    		+ _1_*k4*_03E0_ + _3_*k7*_0E01_ \

			- _03E01_*(k0 + k_1 + k_2 + k_4 + k_7))

		#64 Eq. for _03E02_
        dydt.append(_0_*(k1*_3E02_ + k2*_03E2_) \
    
    		+ _2_*k9*_03E0_ + _3_*k7*_0E02_ \

			- _03E02_*(k_1 + k_2 + k_9 + k_7))

		#65 Eq. for _03E03_
        dydt.append(_0_*(k1*_3E03_ + k2*_03E3_) \
    
    		+ _3_*(k7*_0E03_ + k8*_03E0_) \

			- _03E03_*(k_1 + k_2 + k_7 + k_8))

		#--------------#
		# Return dydt  #
		#--------------#
        return(dydt)
    #----------------------#
	# Integrate the model  #
	#----------------------#
    ds = integrate.odeint(MODEL, initial_cond, t)
    return(ds)
#


# -------------------#
# INITIAL CONDITIONS #
#     Unit: mM       #
#--------------------#

# E, 0, 1, 2, 3, 4
E, G, ATP, ADP = 6.65e-2, 2.5, 3.0, 0.0 

G6P = 2.0
Pi = [2.0, 10.0]

y0 = [[E, G, ATP, G6P, Pi[0], ADP,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0],\
 	 [E, G, ATP, G6P, Pi[1], ADP,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,\
     0.0, 0.0, 0.0, 0.0, 0.0]]
#
#===================================#
# PARAMETER VALUES FOR SIMULATIONS  #
#===================================#

k0 = 63.0

KmG = 0.053

a = 99.0
# For Glucose
k1, k_1 = (a + 1.0)*k0/KmG, a*k0

k2, k_2 = k1, k_1

# For ATP
# Km
KmA = 0.7 
k3, k_3 = (a + 1.0)*k0/KmA, a*k0

# Ki = 
k4, k_4 = k3, k_3

# For G6P

# Ki = 0.71 mM (N)
k5, k_5 = k3, 0.71*k3

# Ki = 54 microM = 0.054 mM(C)
k6, k_6 = k3, 0.054*k3

#d = 1.0e-0
# For Pi
# K_i = 0.022mM
k7, k_7 = k3, 0.022*k3

# Ki = 0.22mM
k8, k_8 = k3, 0.22*k3

# For G6P binding to _3E_, _03E_, and _03E0_
k9, k_9 = 0.1*k3, 0.01*k3

#===========================================#
#  MODEL FOR SENSITIVITY ANALYSIS			#
#===========================================#
## Model for sensitivity analysis
def model(p):
    ys = SOL(p, y0[0], 0.0, 10.05, 0.01)
    #Fs = ys[:,3]
    return(ys[:,3])
    #
#====================================#
a = 0.9
b = 1.1
## Define the problem of SA
problem = {
	'num_vars': 19,
	'names':['k0', 'k1', 'k_1', 'k2', 'k_2', 'k3', 'k_3', \
		'k4', 'k_4', 'k5', 'k_5', 'k6', 'k_6', 'k7', 'k_7', \
		'k8', 'k_8', 'k9', 'k_9'],
	'bounds': [[a*k0, b*k0], \
		[a*k1, b*k1], [a*k_1, b*k_1], \
		[a*k2, b*k2], [a*k_2, b*k_2], \
		[a*k3, b*k3], [a*k_3, b*k_3], \
		[a*k4, b*k4], [a*k_4, b*k_4], \
		[a*k5, b*k5], [a*k_5, b*k_5], \
		[a*k6, b*k6], [a*k_6, b*k_6], \
		[a*k7, b*k7], [a*k_7, b*k_7], \
		[a*k8, b*k8], [a*k_8, b*k_8], \
		[a*k9, b*k9], [a*k_9, b*k_9]]
}
#
#print(problem['bounds'])
#=====================================#
## Generate samples
### Number of samples
n = 1000
##
param_values = saltelli.sample(problem, n,\
			                   calc_second_order=False)

#=====================================#
#
# Data poits for the GSA
# 
#L = [101, 201, 301, 401, 501, 601, 701, 801, 901, 1001]

Y1 = np.zeros([param_values.shape[0]])

Y2 = np.zeros([param_values.shape[0]])

Y3 = np.zeros([param_values.shape[0]])

Y4 = np.zeros([param_values.shape[0]])

Y5 = np.zeros([param_values.shape[0]])

Y6 = np.zeros([param_values.shape[0]])

Y7 = np.zeros([param_values.shape[0]])

Y8 = np.zeros([param_values.shape[0]])

Y9 = np.zeros([param_values.shape[0]])

Y10 = np.zeros([param_values.shape[0]])

#==========================================#
for j, X in enumerate(param_values):
    K      = model(X)
    Y1[j]  = K[101]
    Y2[j]  = K[201]
    Y3[j]  = K[301]
    Y4[j]  = K[401]
    Y5[j]  = K[501]
    Y6[j]  = K[601]
    Y7[j]  = K[701]
    Y8[j]  = K[801]
    Y9[j]  = K[901]
    Y10[j] = K[1001]
    
#==========================================#
Y1
Y2
Y3
Y4
Y5
Y6
Y7
Y8
Y9
Y10

#==========================================#
np.savetxt(str(Pi[0]) + "Pi_Y1_outputs.txt", Y1)
np.savetxt(str(Pi[0]) + "Pi_Y2_outputs.txt", Y2)
np.savetxt(str(Pi[0]) + "Pi_Y3_outputs.txt", Y3)
np.savetxt(str(Pi[0]) + "Pi_Y4_outputs.txt", Y4)
np.savetxt(str(Pi[0]) + "Pi_Y5_outputs.txt", Y5)
np.savetxt(str(Pi[0]) + "Pi_Y6_outputs.txt", Y6)
np.savetxt(str(Pi[0]) + "Pi_Y7_outputs.txt", Y7)
np.savetxt(str(Pi[0]) + "Pi_Y8_outputs.txt", Y8)
np.savetxt(str(Pi[0]) + "Pi_Y9_outputs.txt", Y9)
np.savetxt(str(Pi[0]) + "Pi_Y10_outputs.txt", Y10)

#===========================================#
# Open file to write results
f = open('SA_g6p_' + str(G6P) + '_n_' + str(n) \
		 + '_Pi_' + str(Pi[0]) + '.txt', 'a+')

# Perform analysis Y2
print('Y2')
Si = sobol.analyze(problem, Y2, calc_second_order=False,\
				   print_to_console=True)

### Record the results to the file
f.write('For i = Y2:\n')
for key, value in Si.items():
    f.write('%s: %s\n' % (key, value))

f.write('\n')
f.write('================================\n')

# Perform analysis Y4
print('Y4')
Si = sobol.analyze(problem, Y4, calc_second_order=False,\
				   print_to_console=True)

### Record the results to the file
f.write('For i = Y4:\n')
for key, value in Si.items():
    f.write('%s: %s\n' % (key, value))

f.write('\n')
f.write('================================\n')


# Perform analysis Y6
print('Y6')
Si = sobol.analyze(problem, Y6, calc_second_order=False,\
				   print_to_console=True)

### Record the results to the file
f.write('For i = Y6:\n')
for key, value in Si.items():
    f.write('%s: %s\n' % (key, value))

f.write('\n')
f.write('================================\n')

# Perform analysis Y8
print('Y8')
Si = sobol.analyze(problem, Y8, calc_second_order=False,\
				   print_to_console=True)

### Record the results to the file
f.write('For i = Y8:\n')
for key, value in Si.items():
    f.write('%s: %s\n' % (key, value))

f.write('\n')
f.write('================================\n')

# Perform analysis Y10
print('Y10')
Si = sobol.analyze(problem, Y10, calc_second_order=False,\
				   print_to_console=True)

### Record the results to the file
f.write('For i = Y10:\n')
for key, value in Si.items():
    f.write('%s: %s\n' % (key, value))

f.write('\n')
f.write('================================\n')

# Close result file #
f.close()
  
#=========================================================#
#==			       THE END  							==#
#=========================================================#
