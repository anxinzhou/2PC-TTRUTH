from math import log
import numpy as np 
import math
def log_fit(s,e):
	a=np.linspace(s,e,1000)
	b=np.log(a)/log(2)
	coef= np.polyfit(a,b,1)
	# print(b)
	max_err = max(coef[0]*a+coef[1] - b)
	print("max error:",max_err)
	return coef

def sigmoid_fit(s,e):
	x = np.linspace(s,e,100)
	y = 1/(1+np.power([math.e],-x))
	coef = np.polyfit(x,y,1)
	max_err = max(coef[0]*x+coef[1] - y)
	print("max error:",max_err)
	return coef

def rep_square_root(s,e):
	x = np.linspace(s,e,100)
	y = 1/(x**0.5)
	coef = np.polyfit(x,y,1)
	max_err = max( abs((coef[0]*x+coef[1] - y))/y)
	print("max error:",max_err)
	return coef

# print(log_fit(0.7,0.85))
# print(log_fit(0.85,1))
# print(sigmoid_fit(-5,-3))
# print(sigmoid_fit(-3,-2))
# print(sigmoid_fit(-2,-1))
# print(sigmoid_fit(-1,1))
# print(sigmoid_fit(1,2))
# print(sigmoid_fit(2,3))
# print(sigmoid_fit(3,5))

# print(sigmoid_fit(5,))
print(rep_square_root(0.707,0.85))
print(rep_square_root(0.85,1))