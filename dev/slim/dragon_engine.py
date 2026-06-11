"""
DRAGON engine (self-contained): covariate-parameterized directed structured-
coalescent model for asymmetric gene flow, fit by maximum Wishart likelihood with
an analytic adjoint gradient through the Strobeck solve.

Used by run_scen4_recap.py for rung-3 (SLiM) validation. Depends only on
numpy/scipy. The graph is a DIM x DIM rook stepping stone with elevation =
standardized x-coordinate, matching the SLiM scenarios. Parameter vector
psi = (beta0, gdir, c0, kappa, ltau, lnug):
  M_ab    = exp(beta0 + gdir * (elev_a - elev_b))     directed migration
  gamma_i = exp(c0 + kappa * elev_i)                  coalescence (drift) rate
  Sigma_c = -1/2 exp(ltau) (Co T_obs Co') + exp(lnug) I   Gower + nugget
"""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy import optimize

def make_lattice(DIM):
    d=DIM*DIM
    gx=np.array([i%DIM for i in range(d)],float); gy=np.array([i//DIM for i in range(d)],float)
    elev=(gx-gx.mean())/gx.std()
    edges=[]
    for i in range(d):
        cx,cy=i%DIM,i//DIM
        nb=([i-1] if cx>0 else [])+([i+1] if cx<DIM-1 else [])+([i-DIM] if cy>0 else [])+([i+DIM] if cy<DIM-1 else [])
        edges+=[(i,j) for j in nb]
    ea=np.array([a for a,b in edges]); eb=np.array([b for a,b in edges]); dvec=elev[ea]-elev[eb]
    return dict(d=d,DIM=DIM,gx=gx,gy=gy,elev=elev,ea=ea,eb=eb,dvec=dvec)

def contrast(o):
    M0=np.eye(o)-np.ones((o,o))/o; w,V=np.linalg.eigh(M0)
    return V[:,np.argsort(w)[::-1][:o-1]].T

def build_M_L(beta0,gdir,g):
    d=g['d']; M=np.zeros((d,d)); M[g['ea'],g['eb']]=np.exp(beta0+gdir*g['dvec'])
    return M,np.diag(M.sum(1))-M

def assemble_A(L,gam):
    d=L.shape[0]; Ls=sp.csr_matrix(L); I=sp.identity(d,format='csr')
    di=np.arange(d)*d+np.arange(d); inj=np.zeros(d*d); inj[di]=gam
    return (sp.kron(Ls,I)+sp.kron(I,Ls)+sp.diags(inj)).tocsc()

def forwardg(psi,g,p,Shat_c,obs,Co,zc):
    beta0,gdir,c0,kappa,ltau,lnug=psi; d=g['d']
    tau,nug=np.exp(ltau),np.exp(lnug); gam=np.exp(c0+kappa*zc)
    M,L=build_M_L(beta0,gdir,g)
    try:
        lu=splu(assemble_A(L,gam)); T=lu.solve(np.ones(d*d)).reshape(d,d)
    except Exception:
        return 1e12,None
    T=0.5*(T+T.T); To=T[np.ix_(obs,obs)]; o=len(obs)
    Sig=-0.5*tau*(Co@To@Co.T)+nug*np.eye(o-1)
    sgn,ld=np.linalg.slogdet(Sig)
    if sgn<=0 or not np.isfinite(ld): return 1e12,None
    Sinv=np.linalg.inv(Sig); nll=0.5*p*(ld+np.trace(Sinv@Shat_c))
    return nll,dict(lu=lu,T=T,To=To,M=M,gam=gam,tau=tau,nug=nug,Sinv=Sinv,Co=Co,obs=obs,zc=zc)

def gradg(psi,g,p,Shat_c,obs,Co,zc):
    nll,c=forwardg(psi,g,p,Shat_c,obs,Co,zc)
    if c is None: return 1e12,np.zeros(6)
    d=g['d']; ea,eb,dvec=g['ea'],g['eb'],g['dvec']
    G=c['Sinv']-c['Sinv']@Shat_c@c['Sinv']
    Tbar=np.zeros((d,d)); Tbar[np.ix_(obs,obs)]=-0.25*p*c['tau']*(c['Co'].T@G@c['Co'])
    Lam=c['lu'].solve(Tbar.reshape(-1),trans='T').reshape(d,d)
    Phi=(Lam+Lam.T)@c['T']
    ge=c['M'][ea,eb]*(Phi[ea,ea]-Phi[ea,eb])
    dbeta0=-np.sum(ge); dgdir=-np.sum(ge*dvec)
    dgam=-np.diag(Lam)*np.diag(c['T'])
    dc0=np.sum(dgam*c['gam']); dkappa=np.sum(dgam*c['gam']*c['zc'])
    dtau=-0.25*p*np.trace(G@(c['Co']@c['To']@c['Co'].T)); dnug=0.5*p*np.trace(G)
    return nll,np.array([dbeta0,dgdir,dc0,dkappa,dtau*c['tau'],dnug*c['nug']])

def cov_from_YN(Y,N,maf=0.05):
    freq=Y/N[:,None]; pbar=freq.mean(0); keep=np.minimum(pbar,1-pbar)>=maf
    F=freq[:,keep]; p=F.shape[1]; pb=F.mean(0)
    Gs=(F-pb)/np.sqrt(pb*(1-pb)); return (Gs@Gs.T)/p, p

def _fit(g,p,Shat_c,obs,Co,zc,free,n_start=8,seed=0):
    # NULL defaults for every parameter; non-free entries stay at their null.
    null=np.array([0.0, 0.0, 0.0, 0.0, -1.0, -4.0])   # beta0,gdir,c0,kappa,ltau,lnug
    rng_lo={1:-1.2,2:-1.0,3:-0.8,4:-2.5,5:-7.0}; rng_hi={1:1.2,2:1.5,3:0.8,4:1.5,5:-1.0}
    rng=np.random.default_rng(seed); best=None
    def make(x,free=free):
        psi=null.copy()
        for j,i in enumerate(free): psi[i]=x[j]
        return psi
    for s in range(n_start):
        x0=np.array([rng.uniform(rng_lo[i],rng_hi[i]) for i in free])
        def obj(x): nll,gr=gradg(make(x),g,p,Shat_c,obs,Co,zc); return nll,gr[free]
        r=optimize.minimize(obj,x0,jac=True,method="L-BFGS-B",options=dict(maxiter=400,ftol=1e-11,gtol=1e-7))
        if r.fun<1e11 and (best is None or r.fun<best.fun): best=r
    return best, make

def fit_dragon(Y,N,demeids,DIM,elev_by_deme,maf=0.05):
    g=make_lattice(DIM)
    obs=np.asarray(demeids,int)               # deme ids are graph node indices
    Co=contrast(len(obs)); zc=g['elev']
    Shat_c0,p=cov_from_YN(Y,N,maf)
    Shat_c=Co@Shat_c0@Co.T
    out={"p":p,"n_obs":len(obs)}
    # symmetric (gamma=0, uniform coalescence): free c0,ltau,lnug
    rs,_=_fit(g,p,Shat_c,obs,Co,zc,[2,4,5]); out["nll_sym"]=rs.fun
    # direction only: free gdir,c0,ltau,lnug
    rd,md=_fit(g,p,Shat_c,obs,Co,zc,[1,2,4,5]); out["nll_dir"]=rd.fun; out["gdir_dironly"]=md(rd.x)[1]
    # direction + drift: free gdir,c0,kappa,ltau,lnug
    rj,mj=_fit(g,p,Shat_c,obs,Co,zc,[1,2,3,4,5]); out["nll_joint"]=rj.fun
    pj=mj(rj.x); out["gdir_joint"]=pj[1]; out["kappa_joint"]=pj[3]
    out["LRT_dir_vs_sym"]=2*(out["nll_sym"]-out["nll_dir"])
    out["LRT_joint_vs_dir"]=2*(out["nll_dir"]-out["nll_joint"])
    out["LRT_joint_vs_sym"]=2*(out["nll_sym"]-out["nll_joint"])
    return out
