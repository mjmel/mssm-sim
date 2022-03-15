import numpy as np
import sys
import mpmath as mp
from scipy.special import airy
from scipy.optimize import fsolve
from scipy.optimize import brentq
import mssm_dfe

z0 = -2.33811
sys.stdout.flush()

def solve_for_xc(a: float, b: float, c: float, numind=40):
    """Compute xc given Tc = a, b, and c in MSSM approximation."""

    func2 = lambda xcv : airy((c-xcv)/b)[0]*(a*xcv - 1) - xcv/b*airy((c-xcv)/b)[1]
    
    xc = np.nan
    zc = np.nan
    xcv_arr = np.linspace(0, c-b*z0, numind) 
    f2arr = func2(xcv_arr)
        
    sign_change_up = np.where(np.diff(np.sign(f2arr)) == 2)[0]
    sign_change_down = np.where(np.diff(np.sign(f2arr)) == -2)[0]

    candidate_xc = []
    for sign_change_ind in sign_change_up:
        xc_dict = brentq(func2, xcv_arr[sign_change_ind], xcv_arr[sign_change_ind + 1], full_output=True)
        if xc_dict[-1].converged:
            candidate_xc.append(xc_dict[0])
    for sign_change_ind in sign_change_down:
        xc_dict = brentq(func2, xcv_arr[sign_change_ind], xcv_arr[sign_change_ind + 1], full_output=True)
        if xc_dict[-1].converged:
            candidate_xc.append(xc_dict[0])
    if len(candidate_xc)> 0:
        xc = np.max(np.array(candidate_xc))
        zc = (c-xc)/b - 1/a/b
    return xc, zc

def get_MSSM_prediction(N: int, dfe: mssm_dfe.BaseDFE, scheme: str):
    '''
    Compute theoretical predictions for the rate of adaptation v and other scales using infinitesimal approximation or MSSM approximation.

    Parameters:
        N (int): Number of individuals
        dfe (mssm_dfe.BaseDFE): DFE object (SingleBeneficialDFE, SingleDeleteriousDFE, ExplikeBeneficialDFE, ExplikeDeleteriousDFE, or CompositeDFE)
        scheme (str): Set to 'MSSM' to obtain predictions of MSSM approximation. Set to 'inf' to obtain predictions of infinitesimal approximation.

    Returns:
        Tuple with (v, Tc, b, c) corresponding to the specified approximation. If scheme = 'inf', returned Tc = sigma^2/2D, b = D^{1/3}, c = sigma^4/4D. If scheme = 'MSSM', returned Tc, b, c are as defined in MSSM approximation.
    '''
    func2 = lambda xcv, a, b, c : airy((c-xcv)/b)[0]*(a*xcv - 1) - xcv/b*airy((c-xcv)/b)[1]
    def func_master(a, dfe):
        if scheme == 'inf':
            b = dfe.get_b_inf(a)
            c = dfe.get_c_inf(a)
        elif scheme == 'MSSM':
            b = dfe.get_b(a)
            c = dfe.get_c(a)
            v = dfe.get_v(a)
            M_minus_one = dfe.get_Mp(a, -1)
        
        xc, zc = solve_for_xc(a, b, c, numind=40)

        if not np.isnan(xc):
            if scheme =='inf':
                return 1 - N*xc*np.exp(-(a*b)**3/3.0 + a*b*z0  + 1 ) / airy(zc + 1/a/b)[0] *( (airy(zc + 1/a/b)[1])**2 - (zc + 1/a/b)*(airy(zc+1/a/b)[0])**2 )
            elif scheme == 'MSSM':
                return 1 - N*xc*np.exp(a**2/2.0*v - M_minus_one - a*xc ) / airy(zc + 1/a/b)[0] *( (airy(zc + 1/a/b)[1])**2 - (zc + 1/a/b)*(airy(zc+1/a/b)[0])**2 )
            else:
                return np.nan
        else:
            return np.nan

    func_master_pyfunc = lambda a, dfe : float(np.frompyfunc(func_master, 3, 1)(a, dfe))    
        
    v, a, b, c = np.nan, np.nan, np.nan, np.nan
     
    num_a_guesses_arr = np.array([500])
    sign_change_up, sign_change_down = np.array([]), np.array([])
    for num_a_guesses in num_a_guesses_arr:
        ag_arr = np.logspace(0,np.log10(N),num_a_guesses)
        fmarr = []
        for i in range(len(ag_arr)):
            fmarr.append(func_master(ag_arr[i], dfe))
        fmarr = np.array(fmarr)

        sign_change_up = np.where(np.diff(np.sign(fmarr)) == 2)[0]
        sign_change_down = np.where(np.diff(np.sign(fmarr)) == -2)[0]
        if len(sign_change_up) + len(sign_change_down) > 0:
            break

    if len(np.append(sign_change_up, sign_change_down))>0:
        candidate_a, candidate_b, candidate_c = [],[],[]
        for sign_change_ind in np.append(sign_change_up, sign_change_down):
            a_dict = brentq(func_master, ag_arr[sign_change_ind], ag_arr[sign_change_ind + 1],args=(dfe),full_output=True)
            if a_dict[-1].converged:
                candidate_a.append(a_dict[0])
                if scheme=='inf':
                    candidate_b.append(dfe.get_b_inf(a_dict[0]))
                    candidate_c.append(dfe.get_c_inf(a_dict[0]))
                elif scheme=='MSSM':
                    candidate_b.append(dfe.get_b(a_dict[0]))
                    candidate_c.append(dfe.get_c(a_dict[0]))    
            else:
                print('brentq solution did not converge')

        candidate_a, candidate_b, candidate_c = np.array(candidate_a), np.array(candidate_b), np.array(candidate_c)

        if len(candidate_a) > 0:
            which_select = np.argmax(candidate_a*candidate_b)
            a, b, c = candidate_a[which_select], candidate_b[which_select], candidate_c[which_select]
        else:
            a, b, c = np.nan, np.nan, np.nan
    else:
        if np.sum(np.isfinite(fmarr))>0:
            a_guess = ag_arr[np.nanargmin(abs(fmarr))]
            a_sol, sol_dict, a_converged, msg = fsolve(func_master_pyfunc, a_guess,args=(dfe), full_output=True)        
            if a_converged:
                a = a_sol[0]
            else:
                print('fsolve solution did not converge')
        else:
            a, b, c = np.nan, np.nan, np.nan

    if scheme=='inf':
        b = dfe.get_b_inf(0)
        c = dfe.get_c_inf(a)
        v = dfe.get_v_inf(a)
    elif scheme=='MSSM':
        b = dfe.get_b(a)
        c = dfe.get_c(a)
        v = dfe.get_v(a)              

    return v, a, b, c

