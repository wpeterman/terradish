import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg, spsolve, factorized
from scipy.sparse import diags
np.random.seed(1)

# ---- Build a small synthetic raster grid graph (terradish conventions) ----
nr, nc = 6, 7
N = nr*nc
def cid(r,c): return r*nc + c
edges=[]
for r in range(nr):
    for c in range(nc):
        if c+1<nc: edges.append((cid(r,c),cid(r,c+1)))
        if r+1<nr: edges.append((cid(r,c),cid(r+1,c)))
        # 8-neighbour diagonals
        if r+1<nr and c+1<nc: edges.append((cid(r,c),cid(r+1,c+1)))
        if r+1<nr and c-1>=0: edges.append((cid(r,c),cid(r+1,c-1)))
edges=np.array(edges)  # (E,2), i<j not guaranteed but undirected
demes = np.array([0, 5, 17, 30, 41])   # focal vertices (0-based), m of them
m=len(demes)

# conductance as function of theta via log-linear on 2 covariates
X = np.random.randn(N,2)
def cond(theta):
    return np.exp(X@theta)            # per-vertex conductance, positive
def dcond_dtheta(theta):
    c=cond(theta)
    return c[:,None]*X                # N x p   (dc_i/dtheta_l)

def build_reduced(theta):
    c=cond(theta)
    w = c[edges[:,0]]+c[edges[:,1]]
    # full Laplacian
    import scipy.sparse as sp
    I=np.concatenate([edges[:,0],edges[:,1]])
    J=np.concatenate([edges[:,1],edges[:,0]])
    V=np.concatenate([-w,-w])
    Lf=sp.coo_matrix((V,(I,J)),shape=(N,N)).tolil()
    d=-np.asarray(Lf.sum(1)).ravel()
    Lf.setdiag(d)
    Lf=Lf.tocsr()
    red=np.arange(N-1)
    return Lf[red][:,red].tocsr(), w, c

# RHS Z (constant in theta): columns per deme
def make_Z():
    Z=np.full((N-1,m),-1.0/N)
    for k,dm in enumerate(demes):
        if dm<N-1: Z[dm,k]+=1.0
    return Z
Z=make_Z()

# Selection operator A so that E = A @ G   (m x (N-1))
def make_A():
    A=np.full((m,N-1),-1.0/N)   # baseline -mean over reduced rows
    for j,dm in enumerate(demes):
        if dm<N-1: A[j,dm]+=1.0
    return A
A=make_A()

def forward(theta):
    L,_,_=build_reduced(theta)
    G=spsolve(L,Z)
    if G.ndim==1: G=G[:,None]
    E=A@G
    return L,G,E

# scalar loss: a quadratic that mimics a measurement model on E
Btar=np.random.randn(m,m); Btar=Btar@Btar.T
def loss_and_dE(E):
    R=E-Btar
    l=0.5*np.sum(R*R)
    return l, R           # dl/dE = R

# ---- analytic adjoint gradient (1 forward + 1 adjoint solve, indep of p) ----
def grad_adjoint(theta):
    L,G,E=forward(theta)
    l,dE=loss_and_dE(E)
    dG=A.T@dE                       # (N-1)x m
    Lam=spsolve(L,dG)               # adjoint solve
    if Lam.ndim==1: Lam=Lam[:,None]
    # dl/dL = -Lam @ G^T  (we need contraction with dL/dtheta)
    # dL/dtheta via edges: each edge e=(i,j), w_e=c_i+c_j.
    # contribution of edge to reduced L: for a,b in reduced indices.
    c=cond(theta); dcdt=dcond_dtheta(theta); p=theta.size
    g=np.zeros(p)
    # Build per-edge dl/dw_e = sum over L-entries (dl/dL)*(dL/dw_e)
    # reduced L entry contributions from edge (i,j):
    #   L[i,i]+=w, L[j,j]+=w, L[i,j]-=w, L[j,i]-=w  (only if index<N-1)
    # dl/dL_{ab} = -(Lam G^T)_{ab}; symmetric loss => use M=-(Lam@G.T)
    M = -(Lam@G.T)                  # (N-1)x(N-1) dense (small here)
    red=N-1
    dl_dw=np.zeros(len(edges))
    for e,(i,j) in enumerate(edges):
        s=0.0
        ii=i<red; jj=j<red
        if ii: s+=M[i,i]
        if jj: s+=M[j,j]
        if ii and jj: s+=-M[i,j]-M[j,i]
        dl_dw[e]=s
    # dw_e/dc_i = 1 for both endpoints; accumulate to vertices
    dl_dc=np.zeros(N)
    for e,(i,j) in enumerate(edges):
        dl_dc[i]+=dl_dw[e]; dl_dc[j]+=dl_dw[e]
    g = dcdt.T@dl_dc
    return l,g

def grad_fd(theta,h=1e-6):
    g=np.zeros(theta.size)
    for k in range(theta.size):
        tp=theta.copy(); tp[k]+=h
        tm=theta.copy(); tm[k]-=h
        _,_,Ep=forward(tp); lp,_=loss_and_dE(Ep)
        _,_,Em=forward(tm); lm,_=loss_and_dE(Em)
        g[k]=(lp-lm)/(2*h)
    return g

theta=np.array([0.3,-0.4])
l,ga=grad_adjoint(theta)
gf=grad_fd(theta)
print("loss",round(l,6))
print("adjoint grad",ga)
print("fd grad     ",gf)
print("max rel err ", np.max(np.abs(ga-gf)/(np.abs(gf)+1e-9)))

# ---- block PCG (Jacobi) vs direct, multi-RHS ----
L,_,_=build_reduced(theta)
Gd=spsolve(L,Z); Gd=Gd if Gd.ndim>1 else Gd[:,None]
Mdiag=L.diagonal()
def pcg_block(L,B,tol=1e-10,maxit=5000):
    out=np.zeros_like(B)
    for k in range(B.shape[1]):
        x,info=cg(L,B[:,k],rtol=tol,maxiter=maxit,M=diags(1.0/Mdiag))
        out[:,k]=x
    return out
Gp=pcg_block(L,Z)
print("block-PCG vs direct max err", np.max(np.abs(Gp-Gd)))

# ---- Hutchinson trace estimate sanity: tr(L^{-1}) ----
true_tr=np.trace(np.linalg.inv(L.toarray()))
ests=[]
for t in range(2000):
    z=np.random.choice([-1.0,1.0],size=L.shape[0])
    x=spsolve(L,z)
    ests.append(z@x)
print("Hutchinson tr(Linv): est %.4f true %.4f rel %.3f"%(np.mean(ests),true_tr,abs(np.mean(ests)-true_tr)/true_tr))
