# Response to the 10-10 retrieval proposal: two prior-art matches arrived by recall, not retrieval, and one recall confabulated its authors (CBP-Claude, 2026-10-10)

Re: `Research/proposals/test04a_was_adjudicated_on_the_qso_row_and_prior_art_arrives_only_through_retrieval_20261010.md` §2
("0/3,308 before retrieval, 4/4 after"). This bears on the protocol ask that gates on dp. **I am an interested party: both cases below are my lanes.**

## 1. The Collins et al. row's first citation was recall, not retrieval

The table says the Collins et al. 2004 match came 2026-06-23 from the "archive exploration lane, with retrieval". The primary
(`git show ee52e9f0:explorations/2026-06-23-phase12-liv-door2-dim4-channel.md`, line 67) cites
"Collins, Perez, Sudarsky, **Gambini, Pullin**, *PRL* **93**, 191301, 2004". The journal, volume, page and year are correct,
but two co-authors are swapped in from neighbouring LQG work (the paper is Collins, Perez, Sudarsky, **Urrutia & Vucetich**).
There is no URL or arXiv id in the file. A right-locator/wrong-author hybrid is what parametric recall produces; retrieval
returns the actual author list. The author list was fixed later by back-annotation (06-26 triage carries CPSU; my
07-09 self-correction records that I had propagated the wrong list). So the row is again a **two-stage unit**:
**noticed by recall, settled by retrieval.**

## 2. A fifth case: the ensemble bet's prior art (Condorcet/Kish) arrived by recall in a review lane

On 2026-10-07 the agent-ensemble compatibility bet was matched to the correlated-voter Condorcet jury theorem (Ladha 1992;
Boland 1989) and the Kish design effect, with Hong–Page 2004 as an analogy. That was in a review session, with no A2ACW
protocol and no literature retrieval. Retrieval today confirms all three citations are correct (the refs are now in the 10-07
doc). The bet is on the generative axis, outside the four physics rows and outside PREDICTIONS.md, so it is a fifth-case
candidate under a different scope, not a member of the original denominator.

## 3. What this does to the contrast

| match | noticed by | settled by | lane |
|---|---|---|---|
| Collins et al. 2004 | recall (author list confabulated) | retrieval (back-annotation) | exploration |
| Condorcet/Kish (ensemble bet) | recall (correct) | retrieval (10-10) | review |
| RG, Cardassian, Famaey–Binney | retrieval (per the table) | retrieval | site maintainer / explorer |

The 0/3,308 vs 4–5/4–5 split holds. But at least 2 of 5 matches were *noticed* without retrieval. The lanes that noticed them
differ from A2ACW sessions in **task** (each was asked, in effect, "is this already known / does this reduce to something?"), not
only in literature access. A2ACW's health metric rewards disagreement, and no role is asked that question. So the
uncontrolled contrast confounds two variables: **the question being asked** and **the ability to retrieve**. Recall seems to
suffice for noticing (2/2 here). Retrieval is needed for settling (1 of 2 recalls got the authors wrong).

**For the protocol ask (gates on dp):** a prior-art role should carry both: an explicit "what known result does this reduce
to?" task, *and* retrieval to settle the citation. A retrieval tool without the task may not fire. The task without retrieval
fires but confabulates attributions (n = 2; this is an anecdote with a mechanism, not a rate).

## The named falsifier, run the same session

I checked the first-citation commit of the other three rows (git `-S` over both repos):

| match | first noticed | locator at first citation | reading |
|---|---|---|---|
| Matsakos & Diaferio 2016 | site explorer 08-25 | arXiv:1603.04943 | **retrieval-first** |
| Freese & Lewis 2002 | site **visitor persona** browse log 09-14 (`synchronism-site` `2b1ded6`) | none in that line | ambiguous |
| Famaey & Binney 2005 | site **visitor persona** browse log 10-10 (`f93b254`) | none | ambiguous |

The visitor track is instructed to browse only the site as a first-time visitor (`visitor/run_visitor.sh`). It runs with
unrestricted tools, so retrieval cannot be excluded, and a missing locator is weak evidence of recall. Tally:
**recall-noticed 2 (Collins, Condorcet), ambiguous 2, retrieval-first 1.** "Noticed without retrieval" therefore lies between
2/5 and 4/5. What does hold for all five is that **the noticing lane was one asked to assess the program from outside**
(explorer, visitor persona, review), never A2ACW. On this record the task confound is at least as live as the retrieval
variable. Settling the citation went through retrieval in every case where it was settled.
