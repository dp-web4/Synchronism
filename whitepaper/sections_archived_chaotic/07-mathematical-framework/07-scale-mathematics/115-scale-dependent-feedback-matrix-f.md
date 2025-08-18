# Scale-Dependent Feedback Matrix (F)

The feedback matrix `F` models how interaction outcomes influence
subsequent behavior. To include scale dependence:

    F(r, t, s) = g_s(Î(r, t, s), h_s(t - t'))

where:

-   `F(r, t, s)`: Feedback matrix at position `r`, time `t`, and scale
    `s`.

-   `g_s`: Scale-dependent function relating the interaction tensor `Î`
    to feedback at scale `s`.

-   `h_s(t - t')`: Scale-dependent time-delay function accounting for
    the lag in feedback effects at scale `s`.

**Example:**

At the cellular scale (`s = cellular`), `g_s` and `h_s` could be:

`g_cellular(``Î``, ``Ä``t) = ``â``_c * ``Î`₃` * exp(-``Ä``t / ``ô``_c)`

    h_cellular(Ät) = Ät 

where `â``_c` is a coupling constant, `ô``_c` is a characteristic
cellular feedback time, and `Î`₃ is the alteration component of the
interaction tensor. This implies that feedback strength is proportional
to the alteration caused by the interaction, decays exponentially with
time, and has no significant time delay at the cellular level.