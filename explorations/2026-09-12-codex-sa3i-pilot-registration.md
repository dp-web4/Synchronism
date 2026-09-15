# SA-3I pilot 1: model, transport, and resource registration

Codex, 2026-09-12. dp approved executing one real-participant pilot after the
[instrument handoff](2026-09-12-codex-sa3i-instrument-handoff.md). This supersedes
the original delivery's "no model calls" stop, not its design or outcome rules.
Commit and push this supplement plus the tested runner **before any generation**.
No participant answer has been collected during model-availability discovery.

## Frozen participant and transport

- Installed model identifier: `qwen3.8-distill:4b`, Q6_K; reported parameter
  count 4,205,751,296; architecture `qwen35`.
- Model manifest digest:
  `5f0b0f885f97ea5aa382292093099fc8f472683e7acd1e23a0b625e309fb58aa`.
  This identifies the installed artifact; it is not a claim of official Qwen
  model provenance or equivalence to another quantization.
- Ollama version `0.20.7`, `/api/chat`, `stream=false`, `think=true`,
  `format="json"`; exactly one user message containing the frozen stage prompt.
  No system message, tools, agent wrapper, conversation history, or persistent
  participant memory. Verification uses the frozen verification prompt in a
  fresh request, carrying the original baseline and advice.
- Options: `num_ctx=16384`, `num_predict=4096`, `temperature=0`, `top_p=1`,
  `top_k=0`, `repeat_penalty=1`, `presence_penalty=0`, `frequency_penalty=0`,
  `seed=20260912`. No custom stop strings or response postprocessing.
- Keep the existing model lifecycle policy: no downloads, service reconfiguration,
  explicit unload, or keep-alive override. One request at a time.

The installed model metadata reports a 262,144-token training context and a
GPT-2-style BPE tokenizer (`tokenizer.ggml.model=gpt2`, pre-tokenizer `qwen35`).
We request the smaller 16,384-token runtime context explicitly. Before the first
generation, construct **all** possible stage prompts: four baselines, both
baseline choices for all 64 advice prompts, and both choices for all 64
verification prompts. Reject non-ASCII input for this pilot. Use the conservative
input bound **UTF-8 bytes + 512 tokens** for chat framing; require that bound
plus the 4,096-output-token reserve fit 16,384 for every prompt. BPE merges cannot
increase a byte-token count; the additional framing allowance is an explicit
adapter assumption, not an exact tokenization measurement. No prompt truncation
is implemented or permitted. Check observed prompt counts against this bound
and runtime context after every response; a violation terminates the pilot as
an instrument failure, not as a trust finding.

This verifies constructed-message completeness and conservative context capacity,
not private model attention or server-internal token equality. Save request
hashes and the literal request payloads; report the receipt checks separately.
Reject a changed model digest, server version, tokenizer class, or inadequate
reported context before generation. Recheck the model digest after the run;
if it changed, mark provenance unresolved and make no behavioral attribution.

## Hard resource envelope

- Four baseline generations, then the registered seeded order of 64 episodes.
- At most 132 generations total; at most 4,096 generated tokens per request.
- Total generated tokens at most **131,072**; total reported input plus generated
  tokens at most **400,000**. Reserve each call's conservative input bound and
  full output allowance before sending, then settle against returned counts.
- **1,800 seconds** from the first generation attempt; per-request deadline at
  most 120 seconds and no later than the overall deadline. No parallel requests.
- External inference spend cap: **$0**. Existing compute still has a cost;
  retain wall time and token counters rather than calling inference costless.

No automatic retries, replacement samples, temperature changes, or alternative
model fallback. Do not run capability-generation probes outside the budget.
No new calls after the deadline; a timed-out request is an indeterminate
server-side completion and terminates this pilot (no attempt to cancel another
user's work or restart the service). The request output cap remains in force;
client disconnect alone is not claimed to prove server cancellation.

Persist a start event before each request and a response/error event immediately
afterwards. If interrupted, do not rerun a previously started request. This runner
does not resume; preserve the partial ledger and report an incomplete pilot.
Missing or invalid baseline answers stop before interventions. For an episode,
invalid output is retained as a scored failure; valid `verify` actions buy one
reference follow-up. A malformed object still spelling `action="verify"` is
recorded by the frozen scorer as a requested purchase; preserve the distinction
between this scored request and an actually delivered verification event.

Transport failure, missing/invalid usage telemetry, model mismatch, or a context
bound violation stops new calls. Token counters absent after an error remain
unknown; retain the reserved worst-case charge, never record them as zero.
Budget exhaustion before a follow-up leaves that requested verification without
a valid final answer; all remaining episode IDs stay missing. No complete-grid
contrasts from an incomplete pilot. Syntax-invalid responses, including output
that reaches the generation cap before producing JSON, are never repaired.

## Frozen measurement and scope

Use the instrument at `c2d4ae33`, source SHA-256
`34293a6927bf90f48b3d3239742bc13090e4a457966c111cc5348237d420b34b`,
fixture SHA-256
`7ff57c5aa7227cdbbe2a9353164cfa5a47ada58b7bcda83aeee1cfd35083862e`.
Baseline order C1, C2, I1, I2; episode order is the registered seed 20260912.
Do not replace tasks if baselines hit a ceiling or floor.

Report the original condition-by-advice-correctness cells and per-task paired
contrasts, all protocol failures and receipts, verification requests/deliveries,
and actual run cost. Preserve raw responses, including any thinking field;
score only the untouched assistant content. Do not infer hidden reasoning from
the receipt or count a verification-mediated answer change as direct advice
uptake. Any model/provider/system wording found in raw output is checked before
public publication; operational identifiers must not be published.

JSON-mode and the synthetic histories are part of this intervention. Results
apply to this artifact under these settings on four tasks, not the live fleet,
other models, a population trust effect, or physical Synchronism. No inferential
significance claim or post-hoc challenge-mixture reweighting. Stop after one
pilot, its trace audit, a short result note, and push.

Transport references: [Ollama chat API](https://docs.ollama.com/api/chat) and
[parameter reference](https://docs.ollama.com/modelfile). The
[public distill model card](https://huggingface.co/empero-ai/Qwen3.8-4B-Distill-GGUF)
describes the model family; the installed manifest above is the run's identity.
