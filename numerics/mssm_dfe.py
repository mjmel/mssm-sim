import numpy as np
import sys
import mpmath as mp
from scipy.special import airy
import scipy.integrate as integrate
import scipy.special

z0 = -2.33811
sys.stdout.flush()

def get_exp_factor(a, s, p):
  if p > 0:
    return np.exp(a*s)
  elif p == 0:
    return np.exp(a*s) - 1
  elif p == -1:
    return np.exp(a*s) - a*s - 1
  else: 
    return np.nan

class BaseDFE(object):
    '''
    Base DFE class. Defines helper functions applicable to all DFE types.
    '''
  def get_v(self, a):
    return self.get_Mp(a, 1)

  def get_b(self, a):
    return (self.get_Mp(a, 2) / 2)**(1/3.0)

  def get_c(self, a):
    return a*self.get_Mp(a, 1) - self.get_Mp(a, 0)

  def get_v_inf(self, a):
    return self.get_v(0) + a*self.get_Mp(0, 2)

  def get_b_inf(self, a):
    return self.get_b(0)

  def get_c_inf(self, a):
    return 1/2.0*a**2*self.get_Mp(0, 2)

class SingleBeneficialDFE(BaseDFE):
    '''
    Defines a class for DFE consisting of a single beneficial effect.
    '''
  def __init__(self, Ub, sb):
    '''
        Parameters:
        Ub (float): Beneficial mutation rate 
        sb (float): Effect size
    '''
    self.Ub = Ub
    self.sb = sb

  def get_Mp(self, a, p):
    return self.Ub*self.sb**p*get_exp_factor(a, self.sb, p)

class SingleDeleteriousDFE(BaseDFE):
  def __init__(self, Ud, sd):
    '''
        Parameters:
        Ud (float): Deleterious mutation rate 
        sd (float): Effect size.
    '''
    self.Ud = Ud
    self.sd = sd

  def get_Mp(self, a, p):
    return (-1)**p*self.Ud*self.sd**p*get_exp_factor(a, -self.sd, p)

class ExplikeBeneficialDFE(BaseDFE):
  def __init__(self, Ub, sb, beta_b):
    '''
        Parameters:
        Ub (float): Beneficial mutation rate 
        sb (float): Mean effect size.
        beta_b (float): Steepness parameter beta.
    '''
    self.Ub = Ub
    self.sb = sb
    self.beta_b = beta_b
    self.sigma_b = scipy.special.gamma(1.0/self.beta_b) / scipy.special.gamma(2.0/self.beta_b) * self.sb

  def get_Mp(self, a, p):
    if a == 0:
        if p > 0:
            res = self.Ub*scipy.special.gamma((1+p)/self.beta_b)/scipy.special.gamma(1/self.beta_b)*self.sigma_b**p
        elif p == 0:
            res = 0
        elif p == -1:
            res = 0
    else:
        if self.beta_b < 1:
            res = np.inf
        elif self.beta_b == 1:
            if a*self.sigma_b >= (1-10**-6):
                res = np.inf
            else:
                if p > 0:
                    res = scipy.special.gamma(p+1)*self.Ub*self.sb**p/(1 - a*self.sb)**(p+1)            
                elif p == 0:
                    res = self.Ub*a*self.sb/(1 - a*self.sb)
                elif p == -1:
                    res = self.Ub/self.sb*(np.log(1/(1 - a*self.sb)) - a*self.sb)

        elif self.beta_b > 1:
            if p > 0:
                res, acc_b = scipy.integrate.quad(lambda u: self.Ub*self.sigma_b**p/scipy.special.gamma(1 + 1/self.beta_b)/u**2*(1/u - 1)**p*np.exp(a*self.sigma_b*(1/u-1) - (1/u - 1)**self.beta_b)  ,0,1  )
            elif p == 0:
                res, acc_b = scipy.integrate.quad(lambda u: self.Ub*self.sigma_b**p/scipy.special.gamma(1 + 1/self.beta_b)/u**2*(1/u - 1)**p*np.exp(a*self.sigma_b*(1/u-1) - (1/u - 1)**self.beta_b)  ,0,1  )
                res -= self.Ub
            elif p == -1:
                res, acc_b = scipy.integrate.quad(lambda u: self.Ub*self.sigma_b**p/scipy.special.gamma(1 + 1/self.beta_b)/u**2*(1/u - 1)**p*(mp.exp(a*self.sigma_b*(1/u-1)) - a*self.sigma_b*(1/u-1) - 1)*mp.exp( - (1/u - 1)**self.beta_b)  ,0,1  )
    
    return res

class ExplikeDeleteriousDFE(BaseDFE):
  def __init__(self, Ud, sd, beta_d):
    '''
        Parameters:
        Ud (float): Deleterious mutation rate 
        sd (float): Mean effect size.
        beta_d (float): Steepness parameter beta.
    '''
    self.Ud = Ud
    self.sd = sd
    self.beta_d = beta_d
    self.sigma_d = scipy.special.gamma(1.0/self.beta_d) / scipy.special.gamma(2.0/self.beta_d) * self.sd

  def get_Mp(self, a, p):
    if a ==0:
        if p > 0:
            res = (-1)**p*self.Ud*scipy.special.gamma((1+p)/self.beta_d)/scipy.special.gamma(1/self.beta_d)*self.sigma_d**p
        elif p == 0:
            res = 0
        elif p == -1:
            res = 0
    else:
        if self.beta_d == 1:
            if p > 0:
                res = (-1)**p*scipy.special.factorial(p)*self.Ud*self.sd**p/(1+a*self.sd)**(p+1)
            if p == 0:
                res = -self.Ud*a*self.sd/(1+a*self.sd)
            if p == -1:
                res = -self.Ud/self.sd*(a*self.sd - np.log(1 + a*self.sd)) 
        else:
            if p > 0:
                res_int, acc_d = scipy.integrate.quad(lambda u: u**p*np.exp(-u - (u/a/self.sigma_d)**self.beta_d),0,np.inf  )
                res = (-1)**p*self.Ud/a**(p+1)/self.sigma_d/scipy.special.gamma(1 + 1/self.beta_d) * res_int
            elif p == 0:
                res_int, acc_d = scipy.integrate.quad(lambda u: np.exp(-u - (u/a/self.sigma_d)**self.beta_d),0,np.inf  )
                res = self.Ud/a/self.sigma_d/scipy.special.gamma(1 + 1/self.beta_d) * res_int - self.Ud
            elif p == -1:
                res_int, acc_d = scipy.integrate.quad(lambda u: u**(-1)*(np.exp(-u) - 1 + u)*np.exp( - (u/a/self.sigma_d)**self.beta_d),0,np.inf  )
                res = -self.Ud/self.sigma_d/scipy.special.gamma(1 + 1/self.beta_d) * res_int

    return res

class CompositeDFE(BaseDFE):
  def __init__(self, dfe_arr):
    '''
        Parameters:
        dfe_arr: List of DFEs that make up the composite DFE.
    '''
    self.dfe_arr = dfe_arr
    
  def get_Mp(self, a, p):
    return sum(dfe.get_Mp(a, p) for dfe in self.dfe_arr)