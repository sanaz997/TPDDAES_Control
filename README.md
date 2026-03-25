## TPDDAEs_Stabilization

This MATLAB software package provides stability analysis and stabilization tools for linear time-periodic delay differential algebraic equations (TPDDAEs). The presented toolbox can deal with the time delay systems of the retarded type. These systems are assumed to be given in state-space form. In other words, we study the software toolbox accepts the systems of the form:
 
 $$
\begin{equation}
  \begin{cases}
    E \dot x(t)=\underset{i=0}{\overset{m}\sum}A_i(t)x(t-\tau_{i})+ B u(t),~~t\geq t_0,\\
     y(t)=Cx(t),
\end{cases}
\end{equation}\tag{1}
$$
where:

 - The vectors $`x(t)\in \mathbb{C}^d,~~y(t)\in \mathbb{C}^p,~~\text{and}~~ u(t)\in \mathbb{C}^{r}`$ denote the system state, the measured output, and the control input vectors at time $`t`$. 
 
 - The variable $`t_0`$ represents the initial time. The matrix-valued functions $`A_i:\mathbb{R}\rightarrow \mathbb{R}^{d\times d},~~\text{for}~~i=0,\dots,m,
`$ are real-valued, **time-periodic**, piecewise continuous functions with period $`T>0`$.

- The matrices $`B \in \mathbb{R}^{d\times r},~~ C\in \mathbb{R}^{p\times d}`$ represent the **constant input and output matrices**, respectively.

- The **discrete delay components** are given by the nonnegative parameters $`0=\tau_0<\tau_1 < \dots < \tau_m.`$

- The leading matrix $`E\in \mathbb{R}^{d\times d} `$ is allowed to be **singular**, satisfying $`rank(E)=d-\nu,~~\text{where}~~ 0\leq \nu <d`$.
  
## Assumptions
Note that if $`\nu=0`$, the TPDDAEs **(1)** scales down to time-periodic delay differential equations (TPDDEs). Let us define the columns of $`{Q},{R} \in \mathbb{R}^{d\times \nu}`$ as orthogonal bases for the right and left nullspace of $`E`$, respectively. This implies that $` EQ=0`$ and $`R^T E=0.`$ We have the following assumptions required for our stability analysis approach:

### Assumption 1
Matrices $`A_i(t)`$ satisfy $`R^TA_i(t)Q=0`$ for $`i=1,\dots,m`$, and for all $`t`$.

### Assumption 2
  Matrix $`R^TA_0(t) Q`$ is invertible for all $`t`$, and either $`R^TB=0`$ or $`CQ=0`$.

As illustrated in our paper, these assumptions imply that the introduced TPDDAEs are semi-explicit and of the retarded type.

## Controller Design
In this article, our objective is to design the T-periodic feedback gain $`t\mapsto \mathcal{K}(t)\in \mathbb{R}^{r\times p}`$ such that the 
corresponding closed-loop system obtained from steering **(1)** with a static output feedback controller of the form    
   $$ u(t)=\mathcal{K}(t)y(t) $$ 
is uniformly exponentially stable. At first glance, this structure might seem limiting as the delay components appear exclusively in the state vector, and the matrices $`B`$ and $`C`$are constant; however, the singularity of the matrix $`E`$ allows us to have a rich framework in addressing many control structures, which are demonstrated in our paper.

## Overview
The reader can refer to the comments in the MATLAB file, which demonstrate how the developed toolbox can be used; however, to get a general overview, the following functions are provided:
- `tpddae_create` makes a struct from the provided TPDDAEs.

- `FM_computation` calculates a finite number of dominant Floquet multipliers.
- `tpddae_stabopt_static` designs a static controller for the provided TPDDAEs. The feedback gains can be time-periodic or time-invariant based on the provided input.

- `tpddae_stabopt_dynamic` designs a static controller for the provided TPDDAEs. The feedback gains can be time-periodic or time-invariant based on the provided input.

- `plot_K` plots the temporal profile of the designed optimal feedback gain across a period.

The other functions along with several case studies can be found in the provided toolbox. 

## Application
### 1. Time-periodic systems with time-periodic feedback gain
The provided software provides stability analysis and stabilization methods for TPDDAEs based on the spectrum analysis; one of the most famous example in the study of systems with linear time-periodic dynamics is milling machine. The milling machine is subject to self-regenerative vibrations; therefore, with a proper design of time-periodic feedback gain, we wish to reduce such chatter.

### 2. Time-invariant systems with time-periodic feedback gain
The other application of our software is when you have a system with time-invariant dynamics; however, in order to get better stabilization, you wish to use a time-periodic controller rather than a time-invariant one. It is shown that in many cases time-periodic feedback gain can serve better in system stabilization, for instance, smaller spectral radius, finite-time stabilization, etc. Act-and-wait controller is one of the examples. For a more detailed explanation, you can refer to our paper.

### 3. Time-invariant systems with time-invariant feedback gain
Last but not least, the provided software can as well be used to design time-invariant feedback gain as well. Since time-invariant systems  can be considered periodic with any period, the proposed approach in this paper for their stability analysis works as well. 

## Examples
We provide four case studies involving the design of structured and unstructured dynamic controllers: 
1.  We design static and unstructured dynamic controllers to stabilize the delayed Mathieu equation.
2. The second case study examines the suppression of regenerative chatter during the milling operation using a delay-based feedback controller. 
3. The third case study focuses on the design of a PIR controller in the control problem of a fuel cell system, to provide the required output voltage independent of a pulsating load.
4. The final case study involves the stabilization of an unstable periodic orbit (UPO) of a nonlinear autonomous system via the design of a Pyragas-type feedback controller beyond the odd number limitation. 
You can find these case studies in the subfolder with the name **Case-study**.

## Dependencies
This software uses HANSO to solve non-smooth, non-convex optimization problems arising in the context of controller design. HANSO is available for download from https://cs.nyu.edu/~overton/software/hanso/. HANSO is based on the methods described in [1][2].

## Theory
Interested in reading more and learning about the theory? You can currently refer to [3].

## How to cite?
You can currently cite:

S. Akbarisisi, W. Michiels. **Controller design for systems governed by time-periodic delay differential equations: A spectrum optimization approach.** In Proceedings of the 63rd IEEE Conference on Decision and Control, 2024.

## Authors
This software was written by Sanaz Akbarisisi under the supervision of Wim Michiels.

## Future Extension 
In an extension of our work, not only the feedback gain but also the time delay will serve as controller parameters, leading to better stabilization of TPDDAEs.

## References
[1] A. S. Lewis and M. L. Overton, Nonsmooth optimization via quasi-Newton methods, Mathematical Programming, vol. 141, no. 1, pp. 135–163, 2013. 


[2] J. V. Burke, F. E. Curtis, A. S. Lewis, M. L. Overton, and L. E. Simões, Gradient sampling methods for nonsmooth optimization, Numerical Nonsmooth Optimization, pp. 201–225, 2020.

[3] S. Akbarisisi, W. Michiels. Controller design for systems governed by time-periodic delay differential equations: A spectrum optimization approach. In Proceedings of the 63rd IEEE Conference on Decision and Control, 2024.
