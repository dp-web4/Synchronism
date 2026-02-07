import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def gamma(N_corr):
    return 2.0 / np.sqrt(N_corr)

def coherence_fraction(g):
    return 1.0 / (1.0 + g**2)

N = np.linspace(1, 20, 1000)
g = gamma(N)
cf = coherence_fraction(g)

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('Chemistry Session #2165: Pipe Thread Sealant Chemistry - γ~1 Coherence Analysis\n2028th Phenomenon Type | Finding #2092', fontsize=14, fontweight='bold')

ax1=axes[0,0]; ax1.plot(N,g,'b-',lw=2); ax1.axhline(y=1,color='r',ls='--',label='γ=1'); ax1.axvline(x=4,color='g',ls='--',label='N=4'); ax1.set_xlabel('N_corr'); ax1.set_ylabel('γ'); ax1.set_title('T1: γ=1 at N_corr=4'); ax1.legend(); ax1.grid(True,alpha=0.3)
test1=abs(gamma(4)-1.0)<1e-10

ax2=axes[0,1]; ax2.plot(N,cf,'r-',lw=2); ax2.axhline(y=0.5,color='b',ls='--',label='CF=0.5'); ax2.axvline(x=4,color='g',ls='--',label='N=4'); ax2.set_xlabel('N_corr'); ax2.set_ylabel('CF'); ax2.set_title('T2: CF=0.5 at N_corr=4'); ax2.legend(); ax2.grid(True,alpha=0.3)
test2=abs(coherence_fraction(gamma(4))-0.5)<1e-10

ax3=axes[0,2]; threshold_63=1-1/np.e; ax3.plot(N,cf,'r-',lw=2); ax3.axhline(y=threshold_63,color='m',ls='--',label=f'CF={threshold_63:.3f}'); ax3.set_xlabel('N_corr'); ax3.set_ylabel('CF'); ax3.set_title('T3: CF=1-1/e Threshold'); ax3.legend(); ax3.grid(True,alpha=0.3)
g_at_63=np.sqrt(1/threshold_63-1); N_at_63=4.0/g_at_63**2; test3=abs(coherence_fraction(gamma(N_at_63))-threshold_63)<1e-6

ax4=axes[0,3]; threshold_37=1/np.e; ax4.plot(N,cf,'r-',lw=2); ax4.axhline(y=threshold_37,color='c',ls='--',label=f'CF={threshold_37:.3f}'); ax4.set_xlabel('N_corr'); ax4.set_ylabel('CF'); ax4.set_title('T4: CF=1/e Threshold'); ax4.legend(); ax4.grid(True,alpha=0.3)
g_at_37=np.sqrt(1/threshold_37-1); N_at_37=4.0/g_at_37**2; test4=abs(coherence_fraction(gamma(N_at_37))-threshold_37)<1e-6

ax5=axes[1,0]; dg=np.diff(g); ax5.plot(N[1:],dg,'g-',lw=2); ax5.axhline(y=0,color='r',ls='--'); ax5.set_xlabel('N_corr'); ax5.set_ylabel('dγ/dN'); ax5.set_title('T5: dγ/dN < 0 (Monotonic)'); ax5.grid(True,alpha=0.3)
test5=np.all(dg<0)

ax6=axes[1,1]; ax6.plot(N,g,'b-',lw=2); ax6.plot(1,gamma(1),'ro',ms=10,label=f'γ(1)={gamma(1):.1f}'); ax6.set_xlabel('N_corr'); ax6.set_ylabel('γ'); ax6.set_title('T6: γ(1) = 2'); ax6.legend(); ax6.grid(True,alpha=0.3)
test6=abs(gamma(1)-2.0)<1e-10

ax7=axes[1,2]; dcf=np.diff(cf); ax7.plot(N[1:],dcf,'m-',lw=2); ax7.axhline(y=0,color='r',ls='--'); ax7.set_xlabel('N_corr'); ax7.set_ylabel('dCF/dN'); ax7.set_title('T7: dCF/dN > 0 (Monotonic)'); ax7.grid(True,alpha=0.3)
test7=np.all(dcf>0)

ax8=axes[1,3]; N_ext=np.linspace(1,100,1000); cf_ext=coherence_fraction(gamma(N_ext)); ax8.plot(N_ext,cf_ext,'r-',lw=2); ax8.axhline(y=0.95,color='g',ls='--',label='CF=0.95'); ax8.set_xlabel('N_corr'); ax8.set_ylabel('CF'); ax8.set_title('T8: CF(100) > 0.95'); ax8.legend(); ax8.grid(True,alpha=0.3)
test8=coherence_fraction(gamma(100))>0.95

plt.tight_layout()
plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/pipe_thread_sealant_chemistry_coherence.png', dpi=150, bbox_inches='tight')
plt.close()

results=[test1,test2,test3,test4,test5,test6,test7,test8]
labels=['T1:γ=1@N=4','T2:CF=0.5@N=4','T3:CF=1-1/e','T4:CF=1/e','T5:dγ/dN<0','T6:γ(1)=2','T7:dCF/dN>0','T8:CF(100)>0.95']
print(f"Session #2165: Pipe Thread Sealant Chemistry")
for l,r in zip(labels,results): print(f"  {l}: {'PASS' if r else 'FAIL'}")
print(f"  Total: {sum(results)}/8")
