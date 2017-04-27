from __future__ import division
from scipy.integrate import quad
from math import sqrt
from math import log10
from math import log
import astropy.io.fits as fitsio
import glob
import pylab as P
import numpy as N 
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
from scipy.optimize import fmin
from operator import add
from numpy.linalg import inv
import sys

#Ouverture du fichier de donnees des 740 supernovae:
datafile=open('data/jla_lcparams.txt','r')

#Lecture de la 1ere ligne (noms des colonnes) puis passage a la 2eme:
dataline=datafile.readline()
dataline=datafile.readline()
#   Formation d'une liste composee des elements de la 2eme ligne (1ere supernova) ; "split" separe les donnees des differentes colonnes de la ligne :
liste=dataline.split()
#   Listes a representer graphiquement & initialisation des Variables
zcmb=[]
mB=[]
dmB=[]
X1=[]
dX1=[]
C=[]
dC=[]
M_stell=[]
exp=[]

c=299792.458
H=0.000070
omgM=0.295
alpha=0.141
beta=3.101
M_b=-19.05
delta_M=-0.070

#Boucle -> remplissage des tableaux de donnees extraites des colonnes du tableau des 740 supernovae
while dataline != "" :
#"append" permet d'ajouter un element a la liste a laquelle il est applique (ex: zmcb, mB...), "float" permet de convertir une chaine de nombre (telle que lue dans le tableau des 740 supernovae) et de la convertir en float. 
	zcmb.append(float(liste[1]))
	mB.append((float(liste[4])))
	dmB.append(float(liste[5]))
	X1.append(float(liste[6]))
	dX1.append(float(liste[7]))
	C.append(float(liste[8]))
	dC.append(float(liste[9]))
	M_stell.append(float(liste[10]))
	exp.append(float(liste[17]))
#passage a la ligne suivante du tableau :       
	dataline=datafile.readline()
	liste=dataline.split()


#creation de sigma coh pour etre en accord avec le papier de betoule
sigma_coh=N.zeros(shape=(740))
for i in range(len(exp)):
	if exp[i]==1.0:
		sigma_coh[i]=0.08
	elif exp[i]==2.0:
		sigma_coh[i]=0.11
	elif exp[i]==3.0:
		sigma_coh[i]=0.12
	elif exp[i]==4.0:
		sigma_coh[i]=0.11

#Fonctions

#Fonction donnant les mu theoriques sous forme d'un tableau (fonction de omega_matter et la distance-luminosite)
def intfun(y,x):
	return 1/sqrt(((1+y*x)*((1+x)**2))-(1-y)*x*(2+x))
def fitfundL(z,omgM):
	MU=[]
	for i in range (len(z)): 
	        zz=z[i]
	        MU.append(dL_z(zz,omgM))
	return MU  

def dL_z(zz,omgM): 
	mu_zz = 5*log10((1+zz)*c*(quad(intfun,0,zz,args=(omgM))[0]/(10*H)))
	return mu_zz  

#Fonction de construction des muexp sous forme d'un tableau (formule du papier de 2014)
def muexp(mB,X1,C,alpha,beta,M_b,delta_M,M_stell):
	mu=[]
	for i in range(len(mB)):
	        if M_stell[i]<10:
		       	mu.append(mB[i]-M_b+alpha*X1[i]-beta*C[i])
		else : 
		       	mu.append(mB[i]-M_b-delta_M+alpha*X1[i]-beta*C[i])
	return mu

#Fonction de construction des dmuexp sous forme d'un tableau (formule de propagation des incertitudes)
def dmuexp(dmB,dX1,dC,alpha,beta):
	dmu=[]
	for i in range(len(dmB)):
	        dmu.append(sqrt(dmB[i]**2+(alpha*dX1[i])**2+(beta*dC[i])**2))
	return dmu

#Construction des mu et dmu a l'aide des fonctions muexp,dmuexp
mu=muexp(mB,X1,C,alpha,beta,M_b,delta_M,M_stell)
dmu=dmuexp(dmB,dX1,dC,alpha,beta)
 
#Fonction de construction de la matrice de covariance selon Betoule et al (2014) entre les supernovae pour obtenir le chi2
def mu_cov(alpha, beta): 
	 
	#Assemble the full covariance matrix of distance modulus, See Betoule et al. (2014), Eq. 11-13 for reference
	Ceta = sum([fitsio.open(mat)[0].data for mat in glob.glob('data/covmat/C*.fits')])

	Cmu = N.zeros_like(Ceta[::3,::3])
	for i, coef1 in enumerate([1., alpha, -beta]):
		for j, coef2 in enumerate([1., alpha, -beta]):
			Cmu += (coef1 * coef2) * Ceta[i::3,j::3]

	# Add diagonal term from Eq. 13
	sigma = N.loadtxt('data/covmat/sigma_mu.txt')
	sigma_pecvel = (5 * 150 / 3e5) / (N.log(10.) * sigma[:, 2])
	Cmu[N.diag_indices_from(Cmu)] += sigma[:, 0] ** 2 + sigma[:, 1] ** 2 + sigma_pecvel ** 2
	#Cmu[N.diag_indices_from(Cmu)] += sigma_coh ** 2 + sigma[:, 1] ** 2 + sigma_pecvel ** 2

	return Cmu

#je ne comprend pas le sens de : 
'''
if __name__ == "__main__":
    Cmu = get_mu_cov(0.13, 3.1)
'''       

##Fonction de construction du Chi2 with the covariance matrix implemented as Betoule et al (2014)
def chi2(p,z,mB,dmB,X1,dX1,C,dC,M_stell):
	omgM = p[0]
	alpha = p[1]
	beta = p[2]
	M_b = p[3]
	delta_M = p[4]
	chi2tot=0
	Mcovar = mu_cov(alpha, beta)
	Mat = inv(Mcovar)
	dL=N.zeros(shape=(740))	
	mu_z=muexp(mB,X1,C,alpha,beta,M_b,delta_M,M_stell)

	#boucles de construction des matrices
	for i in range(len(z)):
		zz = z[i]
		dL[i] = dL_z(zz,omgM)


	#instruction de produit matriciel
	chi2tot =  P.dot((P.dot((mu_z-dL),Mat)),((mu_z-dL))) 
	return chi2tot


#Fonction de construction de la matrice de covariance des 5 parametres (en vue d'obtenir les incertitudes sur les 5 parametres)

def covar(param,z,mB,dmB,X1,dX1,C,dC,M_stell):
	# "covar" est initialisee comme matrice de taille 5x5 (nb de parametres x nb de parametres) remplie des
	covar=N.zeros(shape=(len(param),len(param)))
	P1=N.zeros(shape=(len(param)))
	P2=N.zeros(shape=(len(param)))
	P3=N.zeros(shape=(len(param)))
	E=N.zeros(shape=(len(param)))
	#l'erreur est fixee a 1% de la valeur des parametres ; cette valeur peut etre controlee et modifiee
	err = 0.0001

	for i in range(len(param)):
		P1[i]=param[i]
		P2[i]=param[i]
		P3[i]=param[i]

	for i in range(len(param)):
		E[i]=param[i]*err
		P1[i]=P1[i]*(1+err)
	

		for j in range(len(param)):

			if i==j:
				P2[j]=P2[j]*(1+2*err)
				#calcul de la derivee partielle double du chi2 par rapport au parametre i=j
				covar[i][j]=((chi2(P2,z,mB,dmB,X1,dX1,C,dC,M_stell)-2*chi2(P1,z,mB,dmB,X1,dX1,C,dC,M_stell)+chi2(param,z,mB,dmB,X1,dX1,C,dC,M_stell))/(2*E[i]*E[j]))
				
				P2[j]=P2[j]/(1+2*err)
				#print 'ici i =',i,'j=',j
			else:      
				P1[j]=P1[j]*(1+err)
				P2[i]=P2[i]*(1+err)
				P3[j]=P3[j]*(1+err)
				E[j]=param[j]*err
				#calcul de la derivee partielle croisee du chi2 par rapport aux parametres i et j
				covar[i][j]=((chi2(P1,z,mB,dmB,X1,dX1,C,dC,M_stell)-chi2(P2,z,mB,dmB,X1,dX1,C,dC,M_stell)-chi2(P3,z,mB,dmB,X1,dX1,C,dC,M_stell)+chi2(param,z,mB,dmB,X1,dX1,C,dC,M_stell))/(2*E[i]*E[j]))

				P1[j]=P1[j]/(1+err)
				P2[i]=P2[i]/(1+err)
				P3[j]=P3[j]/(1+err)
				#print 'ici i =',i,'j=',j

		P1[i]=P1[i]/(1+err)
		

	#la matrice construite est l'inverse de la matrice de covariance, il convient donc de l'inverser : 	
	covar=inv(covar)
	return covar


#Fonction de calcul des ecarts entre mu_exp et mu_theorique sous forme de tableaux
def ecarts(z,mu,omgM):
	ecarts=[]
	for i in range(len(z)):
		zz = z[i]
		ecarts.append(mu[i]-dL_z(zz,omgM))
	return ecarts





#Fits

#Fit with 1 free parameter
print '*** Fit with 1 free parameter \n'
omG = [0.295]
# "Curve_fit" a un parametre Omega_matter ; elle retourne la valeur du parametre optimale (popt) ainsi que la matrice de covariance (pcov)
popt, pcov = curve_fit(fitfundL,zcmb,mu,omG,dmu,True)
perr = N.sqrt(N.diag(pcov))
print 'Omega_M=', popt[0], ' +-',perr[0]

xfunc = N.linspace(0.001,2,1000)
yfunc = N.zeros(len(xfunc))
omgM = popt[0]

#Representation de la courbe theorique avec le parametre omgm obtenu :
yfunc = fitfundL(xfunc,omgM)

#calcul du chi2 initial
Mcovar=mu_cov(alpha,beta)
print 'ok 1 \n'

chi2curfit=chi2([omgM,alpha,beta,M_b,delta_M],zcmb,mB,dmB,X1,dX1,C,dC,M_stell)
print 'chi2 =',chi2curfit

print 'differences between expected and calculated'
print abs(chi2curfit - 682.9)

#Fit with 5 free parameters
print '\n *** Fit with 5 free parameters \n'
p0=[omgM,alpha,beta,M_b,delta_M]
#Minimisation du chi2 grace a "fmin", en fonction des 5 parametres de fit (retournant les parametres ideaux, mais sans les incertitudes associees)
p = fmin(chi2,p0,args=(zcmb,mB,dmB,X1,dX1,C,dC,M_stell))
print 'p = ',p
sys.exit()
print 'ok 2 \n'

chi25 = chi2(p,zcmb,mB,dmB,X1,dX1,C,dC,M_stell)
print 'chi2 final =', chi25
print 'differences between expected and calculated \n'
print (chi25-682.9)

#Calcul de la matrice de covariance entre les 5 parametres cosmologiques 
#Pcovar=mu_cov(alpha,beta)
Pcovar= covar(p,zcmb,mB,dmB,X1,dX1,C,dC,M_stell)
print 'ok 3 \n'

#Calcul des incertitudes sur les 5 parametres grace a la racine carree des valeurs de la diagonale de la matrice de covariance 
Errors=N.zeros(5)
Errors=N.sqrt(N.diag(Pcovar))

#affichage de la fonction de covariance et des erreurs si besoin
#print 'Pcovar =  \n' , Pcovar, '\n'
print  ' \n Errors = ' , Errors, '\n'

#Affichage des valeurs de fit des 5 parametres et de leurs incertitudes
print '\nOmega_M=', p[0], '+-', Errors[0]
print 'alpha=', p[1], '+-', Errors[1]
print 'beta=', p[2], '+-', Errors[2]
print 'M_b=', p[3], '+-', Errors[3]
print 'delta_M=', p[4], '+-', Errors[4], '\n'

O = 0.295 
dO = 0.034
al = 0.141 
dal = 0.006
bet = 3.101
dbet = 0.075
Mb = -19.05 
dMb = 0.02
delt = -0.070 
ddelt = 0.023


print 'differences between expected and calculated \n'
print abs(p[0] - O)
print abs(p[1] - al )
print abs(p[2] - bet )
print abs(p[3] - Mb)
print abs(p[4] - delt)
print 'differences between expect and calculated for uncertainties \n'
print abs(Errors[0] - dO)
print abs(Errors[1] - dal)
print abs(Errors[2] - dbet)
print abs(Errors[3] - dMb)
print abs(Errors[4] - ddelt)

#Nouveau calcul des valeurs des mu et des dmu avec les parametres ajustes
mufinal=muexp(mB,X1,C,p[1],p[2],p[3],p[4],M_stell)
dmufinal=dmuexp(dmB,dX1,dC,p[1],p[2])

#Calcul de ces ecarts entre experience et theorie
ecarts1=ecarts(zcmb,mu,omgM)
ecarts5=ecarts(zcmb,mufinal,p[0])

    
#Plots

#Premier scatter de points effectue sans corrections des effets (C, M_b, X1) et sans incertitudes, simplement mb en fonction de zcmb
P.figure()
P.scatter(zcmb, mB,s=3)
P.xlabel('zcmb')
P.ylabel('mb')
P.xscale('log')
P.title('Premier scatter')
#Premier fit a 1 parametre suite a l'utilisation de "Curve_fit"
P.figure()
P.scatter(zcmb,mu,s=1)
P.errorbar(zcmb,mu, linestyle='', xerr=0,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P.xlabel('zcmb')
P.ylabel('mu')
P.xscale('log')
P.plot(xfunc,yfunc,'r-')
P.title('Premier fit a 1 parametre')
#Fit a 5 parametres apres minimisation du chi2
P.figure()
P.plot(xfunc,fitfundL(xfunc,p[0]),'r-')
P.errorbar(zcmb,mufinal,linestyle='', xerr=0,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
P.xlabel('zcmb')
P.ylabel('mu')
P.xscale('log')
P.title('Fit a 5 parametres')
#Ecarts entre mu_exp et mu_theorique dans le fit a 1 parametre
P.figure()
P.scatter(zcmb,ecarts1,s=1)
P.errorbar(zcmb,ecarts1, linestyle='', xerr=0,yerr=dmu,ecolor='blue',alpha=1.0,zorder=0)
P.xscale('log')
P.xlabel('zcmb')
P.ylabel('muexp-muth')
P.title('Ecarts - fit a 1 parametre')
#Ecarts entre mu_exp et mu_theorique dans le fit a 5 parametres
P.figure()
P.scatter(zcmb,ecarts5,s=1)
P.errorbar(zcmb,ecarts5, linestyle='', xerr=0,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
#P.plot(zcmb,)
P.xscale('log')
P.xlabel('zcmb')
P.ylabel('muexp-muth')
P.title('Ecarts - fit a 5 parametres')

#Plot du fit a 5 parametres avec les ecarts (type publication 2014)
f, (P1, P2) = P.subplots(2, sharex=True, sharey=False, gridspec_kw=dict(height_ratios=[3,1]))
P1.set_title('Hubble diagram')
P1.scatter(zcmb,mu,label='SN Ia (740)')
dz = 0 #surement besoin de modifier 
P1.errorbar(zcmb,mufinal,linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
P1.set_xlabel('$z$')
P1.set_ylabel('$\\mu = m_b - (M - \\alpha X_1 - \\beta C)$')
P1.set_xscale('log')
P1.set_xlim(0.008,1.2)
P1.set_ylim(30,48)
P1.plot(xfunc,yfunc)
P1.legend(bbox_to_anchor=(0.5, 0.9))
P1.text(0.1, 34, 'string_chi2Ndof', ha='left', rotation=0)
P1.text(0.1, 32, 'string_pM', ha='left', rotation=0)
P2.scatter(zcmb,ecarts5,c='black',s=2)
P2.errorbar(zcmb,ecarts5, linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
P2.set_xscale('log')
P2.set_xlabel('$z$')
P2.set_ylabel('$\\mu - \\mu_{\\Lambda {\\rm CDM}}$')
P2.set_xlim(0.008,1.2)
P2.set_ylim(-1,1)
#P2.plot(zcmb,)
#P.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
f.subplots_adjust(hspace=0.15)

P.show()

#si il a besoin d'enregistrer plus tard la syntaxe sera ici
# savefig('../figures/scatter_ex.png',dpi=48)

#if __name__=='__main__':

	   

