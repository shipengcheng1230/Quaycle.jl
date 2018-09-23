# Quasi-dynamic Simulation

## Basic Theory
The governing equation is that at every time step, shear stress across the fault plane equals to frictional force plus a radiation damping term for approximating wave propagation effect:

```math
τ = σf + ηv
```

Here ``μ`` is shear stress across the fault plain. Using [Okada's dislocation theory](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html), it can be shown as:

```math
τ = K ⊗ δ
```

where ``K`` is the so-called stiffness tensor, depicting relationship between displacements at one position regarding to dislocations somewhere else. ``δ`` is the dislocation, i.e. displacement at everywhere on the fault. ``⊗`` denotes tensor contraction.

Back to ``f``, we use rate-and-state frictional law to calculate its value, specifically as below:

```math
f(v, θ) = f_0 + a \ln{\frac{V}{V_0}} + b \ln{\frac{V_0 θ}{L}}
```

where ``f_0`` and ``v_0`` are reference friction coefficient and velocity, ``v`` and ``θ`` are velocity and state variable based on which frictional force is. ``a`` and ``b`` are two frictional parameters denoting contributions each of which comes from velocity and state variable respectively.

Sometimes people use regularized form to avoid infinity when ``v ≈ 0``, namely:

```math
f(v, θ) = a \sinh ^{-1}{[\frac{v}{2v_0} \exp{\frac{f_0 + b \ln{v_0 θ/L}}{a}}]}
```

There are many state evolution law that describes how state variable ``θ`` changes with time, one of which that most widely used is Dieterich law:

```math
\frac{\mathrm{d}θ}{\mathrm{d}t} = 1 - \frac{v θ}{L}
```

Further, ``η`` is a damping coefficient whose value is often chosen as ``μ / 2vs`` where ``μ`` is shear modulus and ``vs`` shear wave velocity. ``σ`` is the effective normal stress.

To simulate how fault evolves with time, we then take the derivative of the governing equation:

```math
\frac{\mathrm{d} τ}{\mathrm{d} t} = \frac{\mathrm{d} f(v, θ)}{\mathrm{d} t} + η \frac{\mathrm{d} v}{\mathrm{d} t}
= \frac{\mathrm{d} f}{\mathrm{d} v} \frac{\mathrm{d} v}{\mathrm{d} t} + \frac{\mathrm{d} f}{\mathrm{d} θ} \frac{\mathrm{d} θ}{\mathrm{d} t} + η \frac{\mathrm{d} v}{\mathrm{d} t}
```

Thus we arrive at:
```math
\frac{\mathrm{d} v}{\mathrm{d} t} = \frac{\frac{\mathrm{d} τ}{\mathrm{d} t} - \frac{\mathrm{d} τ}{\mathrm{d} θ} \frac{\mathrm{d} θ}{\mathrm{d} t}}{\frac{\mathrm{d} τ}{\mathrm{d} v} + η}
```

where ``\frac{\mathrm{d} τ}{\mathrm{d} t} = K ⊗ (vpl - v)`` where ``vpl`` is the plate rate.

Hence, with both derivatives of velocity ``v`` and state variable ``θ``, we are able to discover how fault evolves with various parameters settings.
