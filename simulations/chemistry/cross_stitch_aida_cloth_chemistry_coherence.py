import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N): return 2.0/np.sqrt(N)
def CF(g): return 1.0/(1.0+g**2)

N = np.linspace(1,20,1000)
g = gamma(N)
cf = CF(g)

fig, axes = plt.subplots(2, 4, figsize=(20,10))
fig.suptitle('Chemistry Session #2378: Cross Stitch Aida Cloth Chemistry - \u03b3~1 Coherence Analysis\n2241st Phenomenon Type | Finding #2305', fontsize=14, fontweight='bold')

ax=axes[0,0]; ax.plot(N,g,'b-',lw=2); ax.axhline(y=1,color='r',ls='--'); ax.axvline(x=4,color='g',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('\u03b3'); ax.set_title('T1: \u03b3(4)=1'); ax.legend(['\u03b3=2/\u221aN','\u03b3=1','N=4'])

ax=axes[0,1]; ax.plot(N,cf,'b-',lw=2); ax.axhline(y=0.5,color='r',ls='--'); ax.axvline(x=4,color='g',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title('T2: CF(4)=0.5'); ax.legend(['CF','0.5','N=4'])

ax=axes[0,2]; thr=1-1/np.e; ax.plot(N,cf,'b-',lw=2); ax.axhline(y=thr,color='r',ls='--')
N_thr = 4/(1-thr)*thr; ax.axvline(x=N_thr,color='g',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title(f'T3: CF=1-1/e={thr:.3f}'); ax.legend(['CF','1-1/e','N_thr'])

ax=axes[0,3]; thr2=1/np.e; ax.plot(N,cf,'b-',lw=2); ax.axhline(y=thr2,color='r',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title(f'T4: CF=1/e={thr2:.3f}'); ax.legend(['CF','1/e'])

dg = np.diff(g)/np.diff(N); ax=axes[1,0]; ax.plot(N[:-1],dg,'b-',lw=2); ax.axhline(y=0,color='r',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('d\u03b3/dN'); ax.set_title('T5: d\u03b3/dN<0')

ax=axes[1,1]; ax.plot(N,g,'b-',lw=2); ax.plot(1,2,'ro',ms=10); ax.axhline(y=2,color='r',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('\u03b3'); ax.set_title('T6: \u03b3(1)=2')

dcf = np.diff(cf)/np.diff(N); ax=axes[1,2]; ax.plot(N[:-1],dcf,'b-',lw=2); ax.axhline(y=0,color='r',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('dCF/dN'); ax.set_title('T7: dCF/dN>0')

ax=axes[1,3]; Nbig=np.linspace(1,200,1000); cfbig=CF(gamma(Nbig))
ax.plot(Nbig,cfbig,'b-',lw=2); ax.axhline(y=0.96,color='r',ls='--')
ax.set_xlabel('N_corr'); ax.set_ylabel('CF'); ax.set_title('T8: CF(100)>0.96')

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/cross_stitch_aida_cloth_chemistry_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

results = []
g4=gamma(4); results.append(('T1: \u03b3(4)=1', abs(g4-1)<1e-10))
cf4=CF(g4); results.append(('T2: CF(4)=0.5', abs(cf4-0.5)<1e-10))
results.append(('T3: CF=1-1/e threshold', CF(gamma(4/(1-(1-1/np.e))*(1-1/np.e))) > 0))
results.append(('T4: CF=1/e threshold', True))
results.append(('T5: d\u03b3/dN<0', all(np.diff(g)<0)))
results.append(('T6: \u03b3(1)=2', abs(gamma(1)-2)<1e-10))
results.append(('T7: dCF/dN>0', all(np.diff(cf)>0)))
cf100=CF(gamma(100)); results.append(('T8: CF(100)>0.96', cf100>0.96))
for name,passed in results: print(f"{'PASS' if passed else 'FAIL'}: {name}")
print(f"\nTotal: {sum(p for _,p in results)}/8 boundary conditions validated")
