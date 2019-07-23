#
# Import some needed packages
#
from scipy.integrate import odeint
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
from _model_ import *
#from _simplest_model_ import *
#from test import *
#
#
#
# Define numerical solutions
# Step size
h = 0.01
#
# Solution
def solution(model, initial_cond, t0, t1, p):
    # p: parameters
    t = np.arange(t0, t1, h)
    sl = integrate.odeint(model, initial_cond, t, args=p)
    return(sl)
#
#
# Time grid
t0 = 0.0
t1 = 8.01
t = np.arange(t0, t1, h)
#t = np.linspace(t0, t11, 1501)
#
# Solutons to the full model
#
ys = []
#
ys00 = []
for i in range(len(y0)):
    ys.append(solution(MODEL, y0[i], t0, t1, p0))
    ys00.append(solution(MODEL, y0[i], t0, t1, p00))
#
# Test conservation laws of the full model
#--# Glucose: _0_
#glc = ys[:,1] + ys[:,3] + ys[:,6] + ys[:,8] + ys[:,10] \
#	\
#	+ ys[:,12] + ys[:,14] + 2*ys[:,15] + ys[:,16] \
#	\
#	+ 2*ys[:,17] + ys[:,18] + 2*ys[:,19] + ys[:,20] \
#	\
#	+ ys[:,21] + ys[:,23] + 2*ys[:,25] + ys[:,26] + ys[:,27] \
#	\
#	+ ys[:,28] + ys[:,30] + ys[:,32] + 2*ys[:,33] \
#	\
#	+ ys[:,34]  + 2*ys[:,35] + ys[:,36] + 2*ys[:,37] \
#	\
#	+ ys[:,38] + 3*ys[:,39] + 2*ys[:,40] + 2*ys[:,41] \
#	\
#	+ 2*ys[:,42] + ys[:,43] + 2*ys[:,44] + ys[:,45] \
#	\
#	+ 2*ys[:,46] + ys[:,47] + 2*ys[:,48] + ys[:,49] \
#	\
#	+ 3*ys[:,50] + 2*ys[:,51] + 2*ys[:,52] + 2*ys[:,53] \
#	\
#	+ ys[:,54] + 2*ys[:,54] + 2*ys[:,55] + ys[:,56] \
#	\
#	+ 2*ys[:,57] + 3*ys[:,58] + 2*ys[:,59] + 3*ys[:,60] \
#	\
#	+ 3*ys[:,61] + 2*ys[:,62] + 3*ys[:,63] + 2*ys[:,64]

# ATP: _1_
#atp = ys[:,2] + ys[:,5]  + ys[:,7] + ys[:,11] + ys[:,14] \
#	\
#	+ ys[:,18] + ys[:,21] + ys[:,22]*2 + ys[:,23] + ys[:,24] \
#	\
#	+ ys[:,26] + ys[:,29] + ys[:,32] + ys[:,35] + ys[:,36]*2 \
#	\
#	+ ys[:,37] + ys[:,38] + ys[:,40] + ys[:,43] + ys[:,46] \
#	\
#	+ ys[:,47]*2 + ys[:,48] + ys[:,49] + ys[:,51] + ys[:,54] \
#	\
#	+ ys[:,57]*2 + ys[:,58] + ys[:,59] + ys[:,60] + ys[:,62]
#
# Pi: _3_
#Pi = ys[:,4] + ys[:,9] + ys[:,13] + ys[:,16] + ys[:,20] \
#	\
#	+ ys[:,24] + ys[:,27] + ys[:,28] + ys[:,29] + ys[:,30] \
#	\
#	+ ys[:,31]*2 + ys[:,34] + ys[:,38] + ys[:,41] + ys[:,42] \
#	\
#	+ ys[:,43] + ys[:,44] + ys[:,45]*2 + ys[:,49] + ys[:,52] \
#	\
#	+ ys[:,53] + ys[:,54] + ys[:,55] + ys[:,56]*2 + ys[:,59] \
#	\
#	+ ys[:,61] + ys[:,62] + ys[:,63] + ys[:,64]*2

#
#---SA illustration---

#plt.figure(figsize=(7, 5))
#for i in range(len(xl)):
#    plt.plot(t, ys[i][:,3], label='%s$k_{-2}$' % xl[i])
#plt.xlabel('Time ($s$)')
#plt.ylabel('$G6P$ concentration $mM$')
#
list_color = ['red', 'black', 'green', 'blue', 'darkorange', 'olive']
#
plt.figure(figsize=(7,5))
plt.plot(t, ys[0][:,3], color=list_color[0], label='$[P_i]$ = ' + str(Pi[0]) + ' $mM$')
plt.plot(t, ys00[0][:,3], color=list_color[0], label='$[P_i]$ = ' + str(Pi[0]) + ' $mM$', linestyle='-.')
for i in range(1, len(y0)):
    plt.plot(t, ys[i][:,3], color=list_color[i], label='$[P_i]$ = ' + str(Pi[i]) + ' $mM$')
    plt.plot(t, ys00[i][:,3], color=list_color[i], label='$[P_i]$ = ' + str(Pi[i]) + ' $mM$', linestyle='-.')
plt.xlabel('Time ($s$)')
plt.ylabel('$G6P$ concentration $mM$')
plt.xticks(np.arange(0.0, t1, 2.0))
#
plt.legend()
plt.show()
#=======================================================#
