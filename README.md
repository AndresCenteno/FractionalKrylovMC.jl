# FractionalKrylovMC

[![Build Status](https://github.com/AndresCenteno/FractionalKrylovMC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AndresCenteno/FractionalKrylovMC.jl/actions/workflows/CI.yml?query=branch%3Amain)

# Quadrature based methods for fractional in space-time equations

We are trying to solve the fractional in space-time heat equation

$$^C_0 D_t^\alpha u = J^{1-\alpha} Du=\int_0^t d\tau\frac{(t-\tau)^{-\alpha}}{\gamma(1-\alpha)} \frac{\partial u}{\partial t}(\tau) = -A^{\beta}u,\qquad u(0)=u_0$$
where $A\in\mathbb{R}^{n\times n}$ is a really big sparse matrix with non-negative eigenvalues (so positive semidefinite) and $\alpha\in(0,1)$, $\beta\in(0,1)$ although it is generalizable to $\alpha\in(0,2)$, $\beta\in(0,\infty)$. The power $\beta$ of a matrix is defined by diagonalizing $A=V\Lambda V^{-1}$ where $\Lambda$ is a diagonal matrix with the diagonal being the spectrum of the matrix $A$ and the columns of $V$ are the eigenvectors
$$AV=\Lambda V$$
and raising the spectrum to $\beta$
$$A^\beta = V \Lambda^\beta V^{-1}$$
as we are dealing with eigenvalues $\lambda\in\mathbb{R}_{\geq 0}$ we don't need to define $a^\beta$ for complex number $a$.

**Important:** this is the part in which the _heat_ is present, all eigenvalues of $-A$ and of $-A^\beta$ are non-positive (non-negative if you don't put a minus) meaning that the value $u(t)$ of the vector will approach a steady state $u(\infty)$ when $t\rightarrow \infty$ (typically either $u(\infty)=0$, $u(\infty)=\bar{u(0)}$ the mean of the initial values or a combination of the eigenvectors $v_i$ of the eigenvalues $\lambda_i$) meaning that there is no periodicity in the solution but total damping. 

There's no need to worry about all this mathematical stuff, the solution is given by 
$$u(t) = E_\alpha(-t^\alpha A^\beta)u(0)$$
where
$$E_\alpha(z) = \sum_{k=0}^\infty \frac{z^k}{\Gamma(\alpha k+1)}$$
is the Mittag-Leffler function (and $\Gamma$ is some sort of continuation of the factorial to a bigger class of numbers than integers).

When dealing with small upper-Hessenberg matrices, or whichever matrix, you either
- Do the Schur-Parlett algorithm, which I don't understand why you would do that.
- Use quadrature rules.
- Use Monte-Carlo methods for the quadrature rules above.
and I'm going to explain you how to do the last two, which are essentially the same.

There is a quadrature rule here https://arxiv.org/abs/0805.3823 (best paper in FractionalDiffEq ever) that says that you can actually compute the Mittag-Leffler function of a matrix by averaging the values of the exponential of the matrix at different times (or sampling from the solution in the classical regime at different end times, however you want to think about it)
$$
E_\alpha(-\lambda A^\alpha)u(0)=\int_0^\infty dr \exp(-rA)u(0)K_\alpha(r,\lambda)
$$
it is just that easy! This function $K_\alpha(r,t^\alpha)$ is called a spectral kernel and I've written more about it in spectral_kernel.pdf. It has these two properties
$$K_\alpha(r,\lambda)\geq 0, \int_0^\infty dr K_\alpha(r,\lambda)=1$$
meaning that it is a probability density function (pdf), and that the integral can be casted as the expected value of 
$$E_\alpha(-t^\alpha A^{\beta})u(0)=\mathbb{E}_r (\exp(-rA^{\beta/\alpha})),\quad r\sim K_\alpha(r,t^\alpha)$$
We still don't want (or can't) compute the power $A^\beta$, or $A^{\beta/\alpha}$, as it is dense (the power $A^n$ of a sparse $A$ will typically be sparse only if $n$ is an integer, this is the equivalent to losing locality in space), so we use another cool property
$$\exp(-A^\beta t)=\int_0^t dr\exp(-Ar) p_\beta(t,r)$$
with $T^\beta_t\sim p_\beta(t,r)$ a non-negative random variable called $\alpha$-stable Lévy subordinator (we are using $\beta$ because $\alpha$ is taken in space). What does this mean? This again means that you only need to care about
- Computing the exponential of you upper Hessenberg matrix
- Adding up the exponential at different times with different weights
        - The weight can be deterministic: meaning that you think of a quadrature rule and add the values up yourself.
        - You obtain weights by sampling, so: Monte-Carlo.

First option means numerically evaluating this integral
$$\int_0^\infty \int_0^\infty dr d\tau \exp(-r^{1/\gamma}\tau A^n)K_\alpha(r,t^\alpha)p_\beta(\tau)$$
where $n$ is the smallest integer such that $\beta/(\alpha n)<1=\gamma$. Second option is just taking the expected value
$$\mathbb{E}_t (\exp(-At))$$ where $t$ follows a weird RNG. This little idea is what makes this worth of publishing.

**This might seem daunting** but I fully understand it, and I left a Julia file fully commented so that you can see the code. There is nothing much to write left! Just experiments. Doing Monte-Carlo is kind of stupid to be completely honest but I don't know a quadrature rule for the Lévy stable subordinator, we can talk about this :P, I do Monte-Carlo because I can't be bothered, but you should do quadratures and use the semigroup property of the exponential
$$\exp((t+s)A)=\exp(tA)\exp(sA)=\exp(sA)\exp(tA)$$
to advance in time and quadratize, or maybe just sample a bunch of points in time in each thread and advance accordingly to those sampled times! This is a good idea.
