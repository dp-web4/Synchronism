import json, time, urllib.request
from collections import defaultdict
import numpy as np

OLLAMA="http://127.0.0.1:11434/api/generate"
# two bias lenses via few-shot worked examples
LENS_A = """Solve each problem the way the examples do: account for every number the problem mentions.

Example: A van carries 8 crates of 12 bottles, pauses 1 hour, then delivers 20 bottles. How many bottles remain on the van?
Work: 8 x 12 = 96; 96 - 20 = 76. The pause affects no bottle count.
Answer: 76

Example: A tank holds 300 L; 40 L are used, then 25 L are added. How much is in the tank?
Work: 300 - 40 + 25 = 285.
Answer: 285
"""
LENS_B = """Solve each problem the way the examples do: first strip the problem to only the quantities that feed the asked amount.

Example: A van carries 8 crates of 12 bottles, pauses 1 hour, then delivers 20 bottles. How many bottles remain on the van?
Work: the pause feeds nothing. Bottles: 8 x 12 = 96 on board; 96 - 20 = 76 after delivery.
Answer: 76

Example: A tank holds 300 L; 40 L are used, then 25 L are added. How much is in the tank?
Work: 300 - 40 + 25 = 285.
Answer: 285
"""

def ask(lens, text, seed):
    prompt=f"{lens}\nProblem: {text}\n\nWork briefly, then the last line must be only the final integer."
    body={"model":"qwen3.5:4b","prompt":prompt,"stream":False,"think":False,
          "options":{"temperature":0.7,"top_p":0.95,"num_predict":64,"seed":seed}}
    req=urllib.request.Request(OLLAMA,data=json.dumps(body).encode(),headers={"Content-Type":"application/json"})
    with urllib.request.urlopen(req,timeout=180) as r: d=json.loads(r.read().decode())
    return d.get("response","")

def parse(t):
    lines=[l.strip() for l in t.strip().splitlines() if l.strip()]
    for ln in reversed(lines):
        try: return int(ln.strip(".,* "))
        except ValueError:
            for tok in reversed(ln.replace(","," ").split()):
                try: return int(tok.strip(".,*"))
                except ValueError: continue
    return None

items=json.load(open("pilot_items.json"))
# LOW arm: 6 agents all LENS_A ; HIGH arm: 6 agents alternating A/B is not diverse...
# real test: correlated-within-family. arms: K1 = all A ; K2 = half A half B
res=defaultdict(dict)
for it in items:
    for arm in ("allA","halfAB"):
        for i in range(6):
            lens = LENS_A if arm=="allA" else (LENS_A if i%2==0 else LENS_B)
            txt=ask(lens, it["text"], seed=hash((it["id"],arm,i))%(2**31))
            got=parse(txt)
            res[(arm,it["id"])][i]= int(got is not None and got==it["answer"])
for arm in ("allA","halfAB"):
    mat=np.array([[res[(arm,it)][i] for i in range(6)] for it in [x["id"] for x in items]],dtype=float)
    print(arm, "acc", round(mat.mean(),3))
    diff=mat.mean(axis=1,keepdims=True); resid=mat-diff; cors=[]
    for i in range(6):
        for j in range(i+1,6):
            a,b=resid[:,i],resid[:,j]
            if a.std()>1e-9 and b.std()>1e-9: cors.append(np.corrcoef(a,b)[0,1])
    print(arm, "residual agreement", round(float(np.mean(cors)),3))
