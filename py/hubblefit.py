from __future__ import print_function
from __future__ import division
from scipy.integrate import quad
from math import *
try:                import pyfits
except ImportError: import astropy.io.fits as pyfits
import glob
import pylab as P
import numpy as N 
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
from scipy.optimize import fmin
from operator import add
from numpy.linalg import inv
from iminuit import Minuit, describe, Struct
#from Plot import *
try:    import cPickle as pickle
except: import pickle
import sys

###########
#constants
###########

clight=299792.458
H=0.000070
omgM=0.295
alpha=0.141
beta=3.101
Mb=-19.05
delta_M=-0.070

############
#Fonctions
############

def savepkl(dic,name='nomame'):
	'''
	Function that create a pkl file to export some data
	inputs:
		-dic : a dictonary containing the data
		-name : the name of the pkl file
	'''
	File = open('Results/pkl/' + name +'.pkl','w')
	pickle.dump(dic,File)
	File.close()

def Sx(X):
	'''
	Function that compute Sx of the distribution given in X
	For the matematical description of Sx, check the pdf file : 'Data/Equations.pdf'
	input : the distribution X
	output : Sx(X)
	'''
	return N.sqrt(abs((1/(len(X)-1))*(sum((X-N.mean(X))**2))))

def RMS(X):
	'''
	Funtction that compute the RMS on a distribution
	inputs :
		-X the distribution
	output :
		-The RMS
	'''
	rms = N.sqrt(abs((1/len(X))*(sum((X-N.mean(X))**2))))
	return rms

def RMSerr(X):
	'''
	Funtction that compute the error on the RMS on a distribution
	inputs :
		-X the distribution
	output :
		-The error of the RMS
	'''
	rmserr = Sx(X)/(N.sqrt(2*len(X)))
	return rmserr

def MEANerr(X):
	'''
	Function that compute the error ont he mean of the distribution given in X
	For the matematical description of Sx, check the pdf file : 'Data/Equations.pdf'
	inputs :
		-X the distibution
	output:
		-the error on the mean
	'''
	meanerr = Sx(X)* (1./N.sqrt(len(X)))
	return meanerr

def Histo(x,mean,sigma_mean,std,stderr,P,orient='horizontal',xlabel='',ylabel=''):
	'''
	Function that plot the histogramm of the distribution given in x
	imputs:
	-x is the distrib itself (array)
	-mean is the mean of the distrib (float)
	-sigma_mean : the error on the average (float)
	-std : the standard deviation (RMS) of the distribution (float)
	-stderr : the errot on the RMS (float)
	-P is the figure where the histogram will be plotted.
	-xylabel and y label are the name ofthe axis.
	'''
	numBins = 20
	P.hist(x,numBins,color='blue',alpha=0.8,orientation=orient,label='average = ' + str("%.5f" % mean) + '$\pm$' + str("%.5f" % sigma_mean) + '\n' +  'rms =' + str("%.5f" % std) + '$\pm$' +str("%.5f" % stderr))
	if xlabel == '':
		P.set_xlabel('number of SNe')
	else:
		P.set_xlabel(xlabel)
	if ylabel == '':
		P.set_ylabel('number of SNe')
	else:
		P.set_ylabel(ylabel)
	P.set_title('Residuals')
	P.legend(bbox_to_anchor=(0.95, 1.0),prop={'size':10})

def comp_rms(residuals, dof, err=True, variance=None):
	"""
	Compute the RMS or WRMS of a given distribution.
	:param 1D-array residuals: the residuals of the fit.
	:param int dof: the number of degree of freedom of the fit.
	:param bool err: return the error on the RMS (WRMS) if set to True.
	:param 1D-aray variance: variance of each point. If given,
               return the weighted RMS (WRMS).
	:return: rms or rms, rms_err
	"""
	if variance is None:                # RMS
		rms = float(N.sqrt(N.sum(residuals**2)/dof))
		rms_err = float(rms / N.sqrt(2*dof))
	else:	  # Weighted RMS
		assert len(residuals) == len(variance)
		rms = float(N.sqrt(N.sum((residuals**2)/variance) / N.sum(1./variance)))
		#rms_err = float(N.sqrt(1./N.sum(1./variance)))
		rms_err = N.sqrt(2.*len(residuals)) / (2*N.sum(1./variance)*rms)
	if err:
		return rms, rms_err
	return rms

def intfun(x,y):
	"""
	Function that build an array that contains theoretical mu (funtion of omgM and luminosity distance)
	imputs:
	-x: represent the redshift
	-y:represent the parameter omgM
	"""
	return 1/sqrt(((1+y*x)*(1+x)**2)-(1-y)*x*(2+x))

def fitfundL(zcmb,omgM):
	"""
	Function which create the distance modulus for each supernovae
	imputs:
	-zcmb is the redshift of each SNe (array)
	-omgM = 0.295.
	outputs:
	-MU is the array containing the distance modulus of each SNe
	"""
	MU=[]
	for i in range (len(zcmb)):
		zz=zcmb[i]
		MU.append(dL_z(zz,zz,omgM))
	return MU  

def dL_z(zcmb,zhel,omgM):
	"""
	Function that compute the integral for the comoving distance.
	imputs:
	-zcmb is the redshift in the cmb framework (array)
	-zhel the redshift in the heliocentric framework (array)
	-omgM = 0.295
	outputs:
	-mu_zz is the array contaning the distance for each SNe
	"""
	mu_zz = 5*log10((1+zcmb)*clight*(quad(intfun,0,zcmb,args=(omgM))[0]/(10*H)))
	#mu_zz = 5*log10((1+zhel)*clight*(quad(intfun,0,zcmb,args=(omgM))[0]/(10*H)))
	#mu_zz = 5*log10(sqrt((1+zhel)*(1+zcmb))*clight*(quad(intfun,0,zcmb,args=(omgM))[0]/(10*H)))
	return mu_zz

def muexp(mB,X1,C,alpha,beta,Mb,delta_M,M_stell):
	"""
	#Correction to muexp regarding the stellar mass of the host galaxy (cf Betoule et al 2014)
	imputs:
	-alpha: free parameter of the Hubble fit (factor of the stretch)
	-beta: free parameter of the Hubble fit (factor of the color)
	-delta_M is a free parameter of the Hubble fit (value of the step for the mass step correction (see Betoule et al 2014))
	-M_stell: the log_10 host stellar mass in units of solar mass
	-mB: B band peak magnitude of the SNs (array)
	-X1:  SALT2 shape parameter (array)
	-C: colour parameter (array)
	-M_stell : the stellar ;aass of each host (array)
	"""
	mu=[]
	#As M.Betoule (choose one or the other)
	for i in range(len(mB)):
		if M_stell[i]<10:
			mu.append(mB[i]-Mb+alpha*X1[i]-beta*C[i])
		else :
			mu.append(mB[i]-Mb-delta_M+alpha*X1[i]-beta*C[i])
	#With fixed delta_M (choose one or the other)
	#Without mass step correction
	'''
	for i in range(len(mB)):
		mu.append(mB[i]-Mb+alpha*X1[i]-beta*C[i])
	'''
	return mu

def dmuexp(dmB,dX1,dC,alpha,beta):
	"""
	Function that builds the list of dmuexp (uncertainties propagation) inputs:
	-dmB: uncertainty on mB
	-dX1: uncertainty on X1
	-dC: uncertainty on C	
	-alpha: free parameter of the Hubble fit (factor of the stretch)
	-beta: free parameter of the Hubble fit (factor of the color)
	"""
	dmu=[]
	for i in range(len(dmB)):
		dmu.append(sqrt(dmB[i]**2+(alpha*dX1[i])**2+(beta*dC[i])**2))
	return dmu

def Remove_Matrix(Tab,ID):
	"""
	function that remove from the 'Tab' matrix all the rows and colomns except those precised in 'ID'
	"""
	#create the list with the line to be removed
	try:
		tab = N.delete(N.arange(len(Tab[0])),ID,0)
	except:
		tab = N.delete(N.arange(len(Tab)),ID,0)

	#remove these lines to the original matrix
	Tab = N.delete(Tab,tab,0)
	try:
		Tab = N.delete(Tab,tab,1)
	except:
		'''nothing else to do'''
	return Tab

def mu_cov(alpha, beta, IDJLA):
	"""
	#Function that buil the covariance matrix as Betoule et al (2014) betzeen SNe to get the chi2
	imputs:
	-alpha: free parameter of the Hubble fit (factor of the stretch)
	-beta: free parameter of the Hubble fit (factor of the color)
	-IDJLA: is an array containing the indexing value of each SNe
	outuput:
	-Cmu : the covariance matrix (2 dimension array)
	"""
	#Assemble the full covariance matrix of distance modulus, See Betoule et al. (2014), Eq. 11-13 for reference
	#You have to acces the data which are in '/data/software/jla_likelihood_v6/covmat/C*.fits'. The C_hosts.fits has been removed from the analysis
	#Ceta = sum([pyfits.getdata(mat) for mat in glob.glob('/data/software/jla_likelihood_v6/covmat/C*.fits')])
	Ceta = sum([pyfits.getdata(mat) for mat in glob.glob('data/covmat/C*.fits')])
	Cmu = N.zeros_like(Ceta[::3,::3])
	for i, coef1 in enumerate([1., alpha, -beta]):
		for j, coef2 in enumerate([1., alpha, -beta]):
			Cmu += (coef1 * coef2) * Ceta[i::3,j::3]

	# Add diagonal term from Eq. 13
	sigma = N.loadtxt('data/covmat/sigma_mu.txt')
	sigma_pecvel = (5 * 150 / 3e5) / (N.log(10.) * sigma[:, 2])
	Cmu[N.diag_indices_from(Cmu)] += sigma[:, 0] ** 2 + sigma[:, 1] ** 2 + sigma_pecvel ** 2
	#Cmu = Remove_Matrix(Cmu,IDJLA)
	#print(len(Cmu))
	#print('coucou')
	return Cmu

def ecarts(zcmb,zhel,mu,omgM):
	"""
	Function that compute the difference between mu_exp and mu_theortical into a list (residuals)
	-zcmb is the redshift in the cmb framework (array)
	-zhel the redshift in the heliocentric framework (array)
	-omgM = 0.295
	-mu :  is the experimental value
	output:
	-ecarts5 is the array containing the residuals
	"""
	ecarts=[]
	for i in range(len(zcmb)):
		z=zcmb[i]
		zz=zhel[i]	
		ecarts.append(mu[i]-dL_z(z,zz,omgM))
	return ecarts

class Chi2: #Class definition
	'''Class for the chi2 '''

	def __init__(self,zcmb,zhel,mB,dmB,X1,dX1,C,dC,M_stell,IDJLA): # Construct method
		'''Attributes'''
		self.IDJLA = IDJLA
		self.chi2tot = 0.
		self.zcmb = zcmb
		self.zhel =  zhel
		self.mB =  mB
		self.dmB =  dmB
		self.X1 = X1
		self.dX1 = dX1
		self.C =  C
		self.dC =  dC
		self.M_stell =  M_stell
		self.dL = N.zeros(shape=(len(IDJLA)))

	def chi2(self,omgM,alpha,beta,Mb,delta_M):
		''' Funtion that calculate the chi2 '''
		result=0.
		Mat = inv(mu_cov(alpha,beta,self.IDJLA))
		mu_z=muexp(self.mB,self.X1,self.C,alpha,beta,Mb,delta_M,self.M_stell)
		#loop for matrix construction
		for i in range(len(self.zcmb)):
			zz = self.zcmb[i]
			zzz = self.zhel[i]
			self.dL[i] = dL_z(zz,zzz,omgM)

		#contruction of the chi2 by matrix product
		result =  P.dot( (mu_z-self.dL), P.dot((Mat),(mu_z-self.dL)))
		self.chi2tot = result
		return result

def Hubble_diagram(zcmb,zhel,mB,dmB,X1,dX1,C,dC,M_stell,IDJLA,exp,results,label='Hubble diagram',DisToCenter='',omgM=0.295,alpha=0.141,beta=3.101,Mb=-19.05,delta_M=-0.070,plot=''):
	"""
	Function which make the hubble diagram of the given compilation.
	inputs :
		-omgM: 1st free parameter to be fitted initialized to the 0.295 value if not precised
		-alpha: 2nd free parameter to be fitted initialized to the 0.141 value if not precised
		-beta: 3rd free parameter to be fitted initialized to the 3.101 value if not precised
		-Mb: 4th free parameter to be fitted initialized to the -19.05 value if not precised
		-delta_M: 5th free parameter to be fitted initialized to the -0.070 value if not precised
		-zcmb: an array which contains the redshifts of the SNs
		-mB: B band peak magnitude of the SNs
		-dmB: uncertainty on mB
		-X1:  SALT2 shape parameter
		-dX1: uncertainty on X1
		-dC: uncertainty on C	
		-C: colour parameter
		-M_stell: the log_10 host stellar mass in units of solar mass
		-IDJLA: index of the SNs from the 740 of jla used for the fit
		-label: string that will be the title of the hubble diagram plot
		-results : a file where some results wiil be written
		-plot : a subplot that need information to get fill. (used only by 'Cut.py')
	outputs:
		- The plot of the Hubble diagram
		- m : the iminuit object (see doc of iminuit for more information).
		-ecarts5 : the residuals of the fit
	"""
	#check : need to have at least 2 SNe
	'''f2,(test,P2,P3,P4,P5) = P.subplots(5, sharex=True, sharey=False, gridspec_kw=dict(height_ratios=[3,1,1,1,1]))'''
	if len(zcmb) == 1 or len(zcmb) == 0 :
		results.write('Not enough data \n')
		return 0
	print(len(zcmb))
	#Definition of the Chi2 object
	chi2mini=Chi2(zcmb,zhel,mB,dmB,X1,dX1,C,dC,M_stell,IDJLA)

	#minimisation of the chi2
	'''
	next instruction can be changed :
	1)the values of the free parameter can be initialized.
	2)the choice between fixed or free parameter can be made (end of the line)
	'''
	#m=Minuit(chi2mini.chi2,omgM=0.295,alpha=0.141,beta=3.101,Mb=-19.05,delta_M=-0.070,limit_omgM=(0.2,0.4),limit_alpha=(0.1,0.2),limit_beta=(2.0,4.0),limit_Mb=(-20.,-18.),limit_delta_M=(-0.1,-0.0),fix_omgM=True, fix_alpha=True, fix_beta=True, fix_Mb=True,fix_delta_M=True, print_level=1)
	m=Minuit(chi2mini.chi2,omgM=0.295,alpha=0.141,beta=3.101,Mb=-19.05,delta_M=-0.070,limit_omgM=(0.2,0.4),limit_alpha=(0.1,0.2),limit_beta=(2.0,4.0),limit_Mb=(-20.,-18.),limit_delta_M=(-0.1,-0.0),fix_omgM=False, fix_alpha=False, fix_beta=False, fix_Mb=False, fix_delta_M=False, print_level=1)
	
	m.migrad()
	#m.hesse()
	#print(m.values)
	#print(m.errors)
	omgM = m.args[0]
	alpha = m.args[1]
	beta = m.args[2]
	Mb = m.args[3]
	delta_M = m.args[4]

	#Pcov = P.subplots()
	#bins_X,bins_Y,contour_X_Y = m.contour('alpha', 'beta', bins=20, bound=1)
	#dic1 = {'alpha':bins_X,'beta':bins_Y,'contour_alpha_beta':contour_X_Y}
	#m.draw_contour('alpha', 'beta', bins=20, bound=1)
	#bins_X,bins_Y,contour_X_Y = m.contour('alpha', 'Mb', bins=20, bound=1)
	#dic1 = {'alpha':bins_X,'Mb':bins_Y,'contour_alpha_Mb':contour_X_Y}
	#m.draw_contour('alpha', 'Mb', bins=20, bound=1)
	#bins_X,bins_Y,contour_X_Y = m.contour('Mb', 'beta', bins=20, bound=1)
	#dic1 = {'Mb':bins_X,'beta':bins_Y,'contour_Mb_beta':contour_X_Y}
	#m.draw_contour('Mb', 'beta', bins=20, bound=1)

	#Save this data into pkl file for further uses
	if len(zcmb)==39:
		dic = {'param':m.values,'cov':m.covariance}
		savepkl(dic,'Ellip')
		#savepkl(dic1,'Ellip1')		
	elif len(zcmb)==88:
		dic = {'param':m.values,'cov':m.covariance}
		savepkl(dic,'Early')
		#savepkl(dic1,'Early1')
	elif len(zcmb)==65:
		dic = {'param':m.values,'cov':m.covariance}
		savepkl(dic,'Late')
		#savepkl(dic1,'Late1')

	#write some results in output
	results.write(str((chi2mini.chi2tot)/(len(IDJLA)-5))+'  ')
	results.write(str(chi2mini.chi2tot) + '  ')

	#Computation of the luminosity-distance modulus and its uncertainty for the minimized parameters
	mufinal=muexp(mB,X1,C,alpha,beta,Mb,delta_M,M_stell)
	dmufinal=dmuexp(dmB,dX1,dC,alpha,beta)

	#Computation of the difference between the measured luminosity-distance modulus and the theorical one (residuals)
	ecarts5=ecarts(zcmb,zhel,mufinal,omgM)

	#computation of the RMS,mean and their error of the distribution
	#std5,std5err = comp_rms(N.array(ecarts5),(len(ecarts5)))
	std5=N.std(ecarts5)
	std5err = RMSerr(ecarts5)
	mean5 = N.mean(ecarts5)
	mean5err = MEANerr(ecarts5)
	results.write(str(mean5)+'  ')
	sigma_mean5 = std5/N.sqrt(len(ecarts5))
	results.write(str(std5)+'  ')

	#mean = N.mean(ecarts5)

	#Computaion of the rms of the fit.
	rms,err_rms = comp_rms(N.array(ecarts5), len(mB), err=True, variance=None)
	results.write(str((omgM))+' ')
	results.write(str(alpha)+' ')
	results.write(str(beta)+' ')
	results.write(str(Mb)+' ')
	results.write(str(delta_M) + '\n')

	#creation of the curve that represents the function of the fit
	xfunc = N.linspace(0.001,2,1000)
	yfunc = N.zeros(len(xfunc))
	yfunc = fitfundL(xfunc,omgM)
	#Plotconfig()#zcmb,mufinal,'$Redshift$','$\\mu = m_b - (M - \\alpha X_1 - \\beta C)$'	)

	#Plot of the fit (a lot of configuration can be added here to make the plot visually better)
	if plot=='':
		f, (P1, P2) = P.subplots(2, sharex=True, sharey=False, gridspec_kw=dict(height_ratios=[3,1]))
		P1.set_title(label + ' with ' + str(len(IDJLA)) + ' SNe',fontsize = 25)
		P.xticks(fontsize=12)
		P1.scatter(zcmb,mufinal,label=('Chi2=' + str("%.2f" % chi2mini.chi2tot) + '\n' + 'Chi2/dof= ' + str("%.2f" % ((chi2mini.chi2tot)/(len(IDJLA)-3))) + '\n' 'Mean =' + str("%.3f" % mean5) + '\n' + 'RMS = ' + str("%.3f" % std5))) # + '$\pm$' + str("%.3f" % (std5/N.sqrt(740)))
		P1.set_ylabel('$\\mu = m_b^{*} - (M_{B} - \\alpha X_1 + \\beta C)$',fontsize=20)
		P2.set_xscale('log')
		P2.set_ylabel('$\\mu - \\mu_{\\Lambda {\\rm CDM}}$',fontsize = 20)
		P.yticks(fontsize=12)
		dz = 0
		P1.errorbar(zcmb,mufinal,linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
		P1.set_xscale('log')
		P1.set_ylim(30,48)
		P1.plot(xfunc,yfunc)
		P1.legend(bbox_to_anchor=(0.5,1.),fontsize=10)
		P2.scatter(zcmb,ecarts5,c='black',s=2)
		P2.errorbar(zcmb,ecarts5, linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)

		#The x and y axis limits have to be changed manually
		P2.set_xlim((min(zcmb)-0.02),(max(zcmb)+log10(1.14)))
		P2.set_ylim(-1,1)

		psfile  = './HubbleDiagram' + '(' + str(len(zcmb))+'SNe)' + '.eps'
		P.savefig(psfile)

		#second plot of the residuals (separately)
		res, (P1, P2) = P.subplots(1,2, sharex=False, sharey=True,gridspec_kw=dict(width_ratios=[3,1]))
		P1.scatter(zcmb,ecarts5,c='black',s=2)
		P1.errorbar(zcmb,ecarts5, linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
		P1.set_xscale('log')
		P1.set_xlabel('$Redshift$')
		P1.set_ylabel('$\\mu - \\mu_{\\Lambda {\\rm CDM}}$')
		P1.set_xlim((min(zcmb)-0.02),(max(zcmb)+log10(2)))
		P1.set_ylim(-1,1)
		P1.set_title(label)
		P1.plot(N.linspace(0.001,1000,1000),N.linspace(0,0,1000))
		res.subplots_adjust(hspace=0.30)
		Histo(ecarts5,mean5,sigma_mean5,std5,std5err,P2)

		psfile  = './residuals' + '(' + str(len(zcmb))+'SNe)' + '.eps'
		P.savefig(psfile)

	elif plot!='':
		plot.scatter(zcmb,mufinal,color='black',s=2,marker='.')#,label=('Chi2=' + str("%.2f" % chi2mini.chi2tot) + '\n' + 'Chi2/dof= ' + str("%.2f" % ((chi2mini.chi2tot)/(len(IDJLA)-3))) + '\n' 'Mean =' + str("%.3f" % mean5) + '\n' + 'RMS = ' + str("%.3f" % std5)),color='black') # + '$\pm$' + str("%.3f" % (std5/N.sqrt(740)))
		plot.yaxis.set_major_locator(MultipleLocator(base=5.))
		plot.yaxis.set_minor_locator(MultipleLocator(1.))
		dz=0
		plot.errorbar(zcmb,mufinal,linestyle='',xerr=dz,yerr=dmufinal,ecolor='black',capsize=2)
		plot.text(0.012,37,'$\chi^2$ = ' + str("%.2f" % chi2mini.chi2tot) + '\n' + '$\chi^2$/dof = ' + str("%.2f" % ((chi2mini.chi2tot)/(len(IDJLA)-3))) + '\n' '$<\Delta \mu >$ = ' + str("%.3f" % mean5) + '\n' + 'RMS = ' + str("%.3f" % std5) ,fontsize=12)
		plot.set_xlim((min(zcmb)-0.02),(max(zcmb)+log10(2)))
		plot.set_ylim(31,44)
		plot.plot(xfunc,yfunc,'k',lw=0.8)		

	if DisToCenter != '':
		#Same plot as previously but for another xaxis so the titles and label are different.
		dist, (P1, P2) = P.subplots(1,2, sharex=False, sharey=True,gridspec_kw=dict(width_ratios=[3,1]))
		P1.scatter(DisToCenter,ecarts5,c='black',s=2)
		P1.errorbar(DisToCenter,ecarts5, linestyle='',xerr=dz,yerr=dmufinal,ecolor='blue',alpha=1.0,zorder=0)
		P1.set_xscale('log')
		P1.set_xlabel('$Distance\ to\ the\ center$')
		P1.set_ylabel('$\\mu - \\mu_{\\Lambda {\\rm CDM}}$')
		P1.set_xlim(min(DisToCenter)-log10(2),max(DisToCenter)+log10(2))
		P1.set_ylim(-1,1)
		P1.set_title(label)
		P1.plot(N.linspace(0.001,1000,1000),N.linspace(0,0,1000))
		res.subplots_adjust(hspace=0.30)
		Histo(ecarts5,mean5,sigma_mean5,std5,std5err,P2)
		P.savefig('./'+label+' - residuals_DisToCenter'+'.png',dpi=300)
		pass
	return m,ecarts5

#######
#Main
#######
if __name__=='__main__':
	#First file is JLA data itself
	jlaData=N.loadtxt('./data/jla_lcparams.txt',dtype='str')
	#jlaData=N.loadtxt('/data/software/jla_likelihood_v6/data/jla_lcparams.txt',dtype='str')

	#Select the chosen data
	SnIDjla = list(jlaData[:,0])
	zcmb = map(float,list(jlaData[:,1]))
	zhel = map(float,list(jlaData[:,2]))
	mB = map(float,list(jlaData[:,4]))
	dmB = map(float,list(jlaData[:,5]))
	X1 = map(float,list(jlaData[:,6]))
	dX1 = map(float,list(jlaData[:,7]))
	C = map(float,list(jlaData[:,8]))
	dC = map(float,list(jlaData[:,9]))
	M_stell= map(float,list(jlaData[:,10]))
	exp=map(float,list(jlaData[:,17]))
	Rajla=map(float,list(jlaData[:,18]))
	Decjla=map(float,list(jlaData[:,19]))
	IDJLA = N.arange(len(jlaData))
	#IDJLA = N.loadtxt('/users/divers/lsst/henne/Desktop/index.csv',dtype='int')

	# Perform the JLA fit and write results
	results=open('./Hubble.txt','w') # file used to write some results (use the one you want instead of mine).
	#results=open('/users/divers/lsst/rosnet/lsst/JLA2014/Hubble.txt','w') # file used to write some results (use the one you want instead of mine).
	results.write('#Chi2/d.o.f. Chi2 mean_of_the_residuals rms Omega_M  alpha beta Mb delta_M \n')
	fit = Hubble_diagram(zcmb,zhel,mB,dmB,X1,dX1,C,dC,M_stell,IDJLA,exp,results,'Hubble diagram, made with JLA dataset,','',omgM,alpha,beta,Mb,delta_M)
	results.close()
	if False: P.show()
	pass
