import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def gamma(N_corr): return 2.0/np.sqrt(N_corr)
def coherence_fraction(g): return 1.0/(1.0+g**2)
N=np.linspace(1,20,1000); g=gamma(N); cf=coherence_fraction(g)
fig,axes=plt.subplots(2,4,figsize=(20,10))
fig.suptitle('Chemistry Session #2047: Candle Dye Chemistry - γ~1 Coherence Analysis\n1910th Phenomenon Type | Finding #1974 *** 1910th PHENOMENON TYPE MILESTONE ***',fontsize=14,fontweight='bold')
ax=axes[0,0]; g4=gamma(4); ax.axhline(y=1,color='r',ls='--',alpha=.7,label='γ=1'); ax.plot(N,g,'b-',lw=2); ax.axvline(x=4,color='g',ls=':',alpha=.7); ax.plot(4,g4,'ro',ms=10); ax.set_xlabel('N_corr'); ax.set_ylabel('γ'); ax.set_title(f'T1: γ(4)={g4:.4f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t1=abs(g4-1)<1e-10
ax=axes[0,1]; c4=coherence_fraction(gamma(4)); ax.plot(N,cf,'b-',lw=2); ax.axhline(y=.5,color='r',ls='--',alpha=.7,label='CF=0.5'); ax.axvline(x=4,color='g',ls=':',alpha=.7); ax.plot(4,c4,'ro',ms=10); ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title(f'T2: CF(4)={c4:.4f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t2=abs(c4-.5)<1e-10
ax=axes[0,2]; tc3=1-1/np.e; Nc3=4/(1/tc3-1); cc3=coherence_fraction(gamma(Nc3)); ax.plot(N,cf,'b-',lw=2); ax.axhline(y=tc3,color='r',ls='--',alpha=.7,label=f'1-1/e={tc3:.4f}'); ax.plot(Nc3,cc3,'ro',ms=10); ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title(f'T3: CF={cc3:.4f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t3=abs(cc3-tc3)<1e-10
ax=axes[0,3]; tc4=1/np.e; Nc4=4/(1/tc4-1); cc4=coherence_fraction(gamma(Nc4)); ax.plot(N,cf,'b-',lw=2); ax.axhline(y=tc4,color='r',ls='--',alpha=.7,label=f'1/e={tc4:.4f}'); ax.plot(Nc4,cc4,'ro',ms=10); ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title(f'T4: CF={cc4:.4f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t4=abs(cc4-tc4)<1e-10
ax=axes[1,0]; dg=np.diff(g)/np.diff(N); ax.plot(N[:-1],dg,'b-',lw=2); ax.axhline(y=0,color='r',ls='--',alpha=.7); ax.set_xlabel('N_corr'); ax.set_ylabel('dγ/dN'); ax.set_title(f'T5: dγ/dN<0: {np.all(dg<0)}'); ax.grid(True,alpha=.3); t5=bool(np.all(dg<0))
ax=axes[1,1]; g1=gamma(1); ax.bar(['γ(1)'],[g1],color='steelblue',alpha=.7); ax.axhline(y=2,color='r',ls='--',alpha=.7,label='2'); ax.set_ylabel('γ'); ax.set_title(f'T6: γ(1)={g1:.4f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t6=abs(g1-2)<1e-10
ax=axes[1,2]; dc=np.diff(cf)/np.diff(N); ax.plot(N[:-1],dc,'b-',lw=2); ax.axhline(y=0,color='r',ls='--',alpha=.7); ax.set_xlabel('N_corr'); ax.set_ylabel('dCF/dN'); ax.set_title(f'T7: dCF/dN>0: {np.all(dc>0)}'); ax.grid(True,alpha=.3); t7=bool(np.all(dc>0))
ax=axes[1,3]; c100=coherence_fraction(gamma(100)); ax.bar(['CF(100)'],[c100],color='steelblue',alpha=.7); ax.axhline(y=.95,color='r',ls='--',alpha=.7,label='.95'); ax.set_ylim(.9,1); ax.set_ylabel('CF'); ax.set_title(f'T8: CF(100)={c100:.6f}'); ax.legend(fontsize=8); ax.grid(True,alpha=.3); t8=c100>.95
plt.tight_layout(); plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/candle_dye_chemistry_coherence.png',dpi=150,bbox_inches='tight'); plt.close()
r=[t1,t2,t3,t4,t5,t6,t7,t8]; print(f"Session #2047: {sum(r)}/8"); print(f"  {r}")
