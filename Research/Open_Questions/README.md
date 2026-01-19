# Open Questions

**Status**: Active exploration queue
**Purpose**: Track unresolved theoretical questions flagged for formal investigation

---

## Current Questions

### 1. Coherence Backpropagation
**File**: `Coherence_Backpropagation.md`
**Origin**: LeCun 1988 Lagrangian backprop framework discussion

**Question**: The coherence function C(x) = tanh(γ × g(x)) is forward-only. Is there a feedback mechanism analogous to backpropagation?

**Key Insight**: Backprop doesn't require time reversal - error signals inform the *next* present state. Could coherence dynamics include similar error-based biasing?

**Status**: Exploration phase - needs formal mathematical treatment

---

### 2. Gnosis-EP Connection
**File**: `Gnosis_EP_Connection.md`

**Question**: How does Gnosis's self-awareness mechanism connect to Epistemic Proprioception (EP) in SAGE?

**Status**: Tracked for integration work

---

## Adding New Questions

When flagging a question for investigation:

1. Create a new `.md` file with descriptive name
2. Include:
   - Date and origin of question
   - Clear statement of the question
   - Current understanding/partial answers
   - What investigation would look like
   - Potential implications if resolved

---

## Resolution Process

Questions move from here to:
- **Session logs**: If they become formal research sessions
- **Framework docs**: If resolved and integrated into theory
- **Archive**: If determined to be outside scope or already answered
