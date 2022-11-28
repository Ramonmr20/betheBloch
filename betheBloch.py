import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.constants import c
import fit_config as fit

#Constantes
pi = 3.1415926535
Na = 6.022e23 	# mol^-1
re = 2.8179e-15 # m
me = 0.511 		#MeV
C0 = 0.1535 	#MeV g^-1 cm^2 <-- es el resultado de mult las ctes al principio de la ec

beta2 = lambda bg: bg**2/(bg**2 + 1) # bg es BetaGamma
gamma = lambda bg: np.sqrt(1 + bg**2)

def C(Z,hvp):
	return -2*np.log(I(Z)/hvp) - 1

def C_shell(bg,Z):
	aux1 = (0.422377*bg**(-2) + 0.0304043*bg**(-4) - 0.00038106*bg**(-6))*1e-6*I(Z)**2
	aux2 = (3.850190*bg**(-2) - 0.1667989*bg**(-4) + 0.00157955*bg**(-6))*1e-9*I(Z)**3

	return aux1 + aux2

def I(Z):
	if Z<13:
		aux =  12*Z + 7
	else:
		aux = 9.76*Z + 58.8*Z**(-0.19)

	return aux*1e-6 #Paso de eV a MeV
		
def Wmax(bg,M):
	return (2*me*bg**2)/(1 + 2*gamma(bg)*me/M + (me/M)**2)

def delta(bg,Z,M,delta0,x0,x1,a,m,hvp): #delta0,x0,x1,a,k estan tabulados
	x = np.log10(bg)

	if x >= x1:
		aux = 4.6052*x + C(Z,hvp)
	elif x >= x0:
		aux = 4.6052*x + C(Z,hvp) + a*(x1 - x)**m
	else:
		aux = delta0*10**(2*(x-x0))
	
	#print(x0,'<',x,'<',x1,'||',aux,'C',C(Z,hvp))

	return aux
	
	

def dEdx(bg,z,M,a,m,x0,x1,delta0,Z,A,hvp):
	aux = C0*(Z/A)*(z**2/beta2(bg))*(np.log(2*me*bg*bg*Wmax(bg,M)/I(Z)**2) - 2*beta2(bg) - delta(bg,Z,M,delta0,x0,x1,a,m,hvp)
  - C_shell(bg,Z)/Z)
	
	#print('log:  ',np.log(2*me*bg*bg*Wmax(bg,M)/I(Z)))
	#print('beta2: ', - 2*beta2(bg))
	#print('delta',- delta(bg,Z,M,delta0,x0,x1,a,m,hvp) )
	#print('cshell',- C_shell(bg,Z)/Z)
	#print(aux)

	return aux

 
# Graficas para distintos materiales

'''
print('#'*30)
print('Tipo de material a considerar')
print('1 -- Metales')
ind = input('Introduce el número: ')
print('#'*30)
'''
ind = '1'####

data = []
if ind == '1':
	data = np.genfromtxt('dataMet.txt',dtype=None)
else:
	print('Error')
	exit()


ind = 15
a = float(data[ind][1])
m = float(data[ind][2])
x0 = float(data[ind][3])
x1 = float(data[ind][4])
delta0 = float(data[ind][5])
Z = float(data[ind][6])
A = float(data[ind][7])
hvp = float(data[ind][8])*1e-6
rho = float(data[ind][9])




'''
print('#'*30)
print('Material a considerar')
for ii in range(len(data)):
	print(ii,' -- ',data[ii][0])
ind = int(input('Introduce el número: '))
print('#'*30)
'''

def plotBethe():
	a = float(data[ind][1])
	m = float(data[ind][2])
	x0 = float(data[ind][3])
	x1 = float(data[ind][4])
	delta0 = float(data[ind][5])
	Z = float(data[ind][6])
	A = float(data[ind][7])
	hvp = float(data[ind][8])*1e-6
	rho = float(data[ind][9])
	
	
	'''
	print('#'*30)
	print('Introduce los siguientes valores:')
	z = float(input('Carga de la partícula incidente: '))
	M = float(input('Masa de la partícula indicente: '))
	bg0 = float(input('BetaGamma minimo: '))
	bgm = float(input('BetaGamma max: '))
	'''
	bg0 = 0.05
	bgm = 1e5
	
	# Grafica Bethe-Bloch
	min_exponent = np.log10(bg0)
	xx = np.logspace(min_exponent,5,1000)
	yy = []
	for ii in xx:
		yy.append(dEdx(ii,z,M,a,m,x0,x1,delta0,Z,A,hvp)*rho)
	
	return xx,yy
	
#print(xx)
#print(yy)

plt.figure(0)



fit.figure_features()
ax = plt.axes()
'''
ind=0 # Li
xx,yy = plotBethe()
plt.plot(xx,yy,'g',lw=3,label='Li')

ind=15 # Ni
xx,yy = plotBethe()
plt.plot(xx,yy,'r',lw=3,label='Ni')

ind=16 # Ni
xx,yy = plotBethe()
plt.plot(xx,yy,'b',lw=3,label='Pyrex Glass')
'''

ind=15 # Ni
#Muon
z = 1
M = 105
xx,yy = plotBethe()
plt.plot(xx,yy,'r',lw=3,label='Muón')
#B+c
z = 1
M = 6274
xx,yy = plotBethe()
plt.plot(xx,yy,'g',lw=3,label=r'Mesón $B^+_c$')
#Alpha
z = 2
M = 3723
xx,yy = plotBethe()
plt.plot(xx,yy,'b',lw=3,label=r'$\alpha$')


fit.add_grid(ax)
plt.yscale('log')
plt.xscale('log')

plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel(r'$\beta\gamma$',fontsize=40)
plt.ylabel(r'$-\left<\frac{dE}{dx}\right>$ (MeV cm$^-1$)',fontsize=40)
plt.legend(fontsize=30)

#Grafica Bragg
def bragg():
	ind = 15
	a = float(data[ind][1])
	m = float(data[ind][2])
	x0 = float(data[ind][3])
	x1 = float(data[ind][4])
	delta0 = float(data[ind][5])
	Z = float(data[ind][6])
	A = float(data[ind][7])
	hvp = float(data[ind][8])*1e-6
	rho = float(data[ind][9])

	E = 5000
	M = 100
	paso = 1e-12 #s
	bg = np.sqrt(E**2/M**2 - 1)
	d_old = 0
	d_a = [0]
	E_a = [E]
	dE_a = [0]
	while E>M:
		d = paso*c*100*bg/(np.sqrt(bg**2 + 1))
		#print(d)
		
		perdida = dEdx(bg,z,M,a,m,x0,x1,delta0,Z,A,hvp)*rho
		if (random.uniform(0,1)>d):
			perdida = 0
		else:
			if (perdida>(E-M)):
				perdida = E
				E = M
				bg = 0
			else:
				E -= perdida
				bg = np.sqrt(E**2/M**2 - 1)

		d_a.append(d_old+d)
		E_a.append(E)
		dE_a.append(perdida)
		d_old += d
		
	return d_a,E_a,dE_a
	
plt.figure(1)	
'''
d_array = []
E_array = []
dE_array = []
for ii in range(100):
	d_a,E_a,dE_a = bragg()
	d_array.append(d_a)
	E_array.append(E_a)
	dE_array.append(dE_a)


d_plot = np.linspace(0,110,40)
E_plot = [0]*len(d_plot)
dE_plot = [0]*len(d_plot)
counter = [0]*len(d_plot)

for p in range(len(d_array)):
	print(p)
	for ii in range(len(d_array[p])):
		for jj in range(len(d_plot)):
			if d_array[p][ii] < d_plot[jj]:
				E_plot[jj] += E_array[p][ii]				
				dE_plot[jj] += dE_array[p][ii]				
				counter[jj] += 1
				break
'''

def perdidaE(ind,M,z):
	a = float(data[ind][1])
	m = float(data[ind][2])
	x0 = float(data[ind][3])
	x1 = float(data[ind][4])
	delta0 = float(data[ind][5])
	Z = float(data[ind][6])
	A = float(data[ind][7])
	hvp = float(data[ind][8])*1e-6
	rho = float(data[ind][9])

	
	E = 10000
	paso = 1e-11 #s
	bg = np.sqrt(E**2/M**2 - 1)
	d_old = 0
	d_a = [0]
	E_a = [E]
	dE_a = [0]
	while E>M:
		d = paso*c*100*bg/(np.sqrt(bg**2 + 1))
		print(d)
		
		perdida = dEdx(bg,z,M,a,m,x0,x1,delta0,Z,A,hvp)*rho*d
		if (perdida>(E-M)):
			perdida = E
			E = M
			bg = 0
		else:
			E -= perdida
			bg = np.sqrt(E**2/M**2 - 1)

		d_a.append(d_old+d)
		E_a.append(E)
		dE_a.append(perdida)
		d_old += d
		
	return d_a,E_a,dE_a

'''	
d_a,E_a,dE_a = perdidaE(15,105,1)
fit.figure_features()
ax = plt.axes()


plt.plot(d_a,E_a,label='Muón',lw=3)
fit.add_grid(ax)

d_a,E_a,dE_a = perdidaE(15,1050,2)
plt.plot(d_a,E_a,label=r'$\alpha$',lw=3)

d_a,E_a,dE_a = perdidaE(15,6000,1)
plt.plot(d_a,E_a,label=r'$B^+_c$',lw=3)
'''

#plt.plot(d_a,E_a,label='Niquel',lw=3)

ax = plt.axes()
d_a,E_a,dE_a = perdidaE(1,105,1)
#plt.plot(d_a,E_a,label=r'Litio',lw=3)

d_a,E_a,dE_a = perdidaE(13,105,1)
plt.plot(d_a,E_a,label='Muón en Hierro',lw=3)
print(d_a[-1])

fit.add_grid(ax)
plt.xlabel(r'$x$ (cm)',fontsize=40)
plt.ylabel(r'$E$ (MeV)',fontsize=40)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

plt.yscale('log')
plt.legend(fontsize=30)
plt.show()



