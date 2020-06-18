# Algorithm 675: Fortran subroutines for computing the square root covariance filter and square root information filter in dense or Hessenberg forms

For the original description please refer to
* Vanbegin, M., & Verhaegen, M. (1989). Algorithm 675: FORTRAN subroutines for computing the square root covariance filter and square root information filter in dense or Hessenberg forms. _ACM Transactions on Mathematical Software (TOMS)_, 15(3), 243-256. doi:[10.1145/66888.69647](https://doi.org/10.1145/66888.69647)

A PDF version of the original article free of charge is available [here](https://dl.acm.org/doi/pdf/10.1145/66888.69647). The original source code can be downloaded from [Netlib](http://www.netlib.org/toms/).

> **_NOTE:_** The code available here is **_not_** the original Software and contains modifications (see differences below). Use of this code is still subject to the [ACM Software License Agreement](https://www.acm.org/publications/policies/software-copyright-notice). This agreement grants the right to execute, copy, modify and distribute both the binary and source code solely for academic, research, and other similar noncommercial uses. See the [LICENSE](./LICENSE) document for the remaining conditions. New files (i.e. those which are not present in the original package from Netlib) are provided under terms of the permissive MIT License.

## Routines

### Square Root Covariance Filter

The calling sequence is
```fortran
call srcf(s,lds,a,lda,b,ldb,q,ldq,c,ldc,r,ldr,n,m,p,k,ldk,wrk,ldw,multbq,withk,tol)
```
with parameters as follows:
* `s`: A two-dimensional real array containing the lower triangular square root of the `n`✕`n` state covariance matrix $P_{k|k-1}$. Upon return it contains the corresponding square root of $P_{k+1|k}$.
* `lds`: The first (integer) dimension of `s` as declared in the calling program.
* `a`: A two-dimensional real array containing the `n`✕`n` state transition matrix $A_k$ of the discrete-time system.
* `lda`: The first (integer) dimension of `a` as declared in the calling program.
* `b`: A two-dimensional real array containing the `n`✕`m` input weight matrix $B_k$ (or its product with Q if `multbq=.true.`).
* `ldb`: The first (integer) dimension of `b` as declared in the calling program.
* `q`: A two-dimensional real array containing the `m`✕`m` lower triangular Cholesky factor $Q_k^L$ of the process noise covariance matrix.
* `ldq`: The first (integer) dimension of `q` as declared in the calling program.
* `c`: A two-dimensional real array containing the `p`✕`n` output weight matrix $C_k$ of the discrete-time system.
* `ldc`: The first (integer) dimension of `c` as declared in the calling program.
* `r`: A two-dimensional real array containing the `p`✕`p` lower triangular Cholesky factor $R_k^L$ of the measurement noise covariance matrix.
* `ldr`: The first (integer) dimension of `r` as declared in the calling program.
* `n`: Dimension (integer) of the state (actual value used in `srcf`).
* `m`: Dimension (integer) of the input (actual value used in `srcf`).
* `p`: Dimension (integer) of the output (actual value used in `srcf`).
* `k`: A two-dimensional real array containing the `n`✕`p` Kalman gain $K_k$ upon return (if `withk=.true.`).
* `ldk`: The first (integer) dimension of `k` as declared in the calling program.
* `wrk`: A real working array of dimensions `ldw`✕`nw` containing the postarray upon return. `nw` must be at least `(m+n+p)`.
* `ldw`: The first (integer) dimension of `wrk`, which must be at least `(n+p)`.
* `multbq`: Logical variable indicating how input matrices `b` and `q` are transferred. If `multbq =.true.`, then the product `b*q` is transferred via `b`, and `q` is not used. If `multbq =.false.`, then `b` and `q` are transferred via the parameters `b` and `q`, respectively.
* `withk`: Logical variable indicating whether or not `k` is requested. If a nearly singular matrix is encountered during the calculation of `k`, `withk` is set to `.false.` upon return.
* `tol`: Real indicating the tolerance for the reciprocal condition number of the computed $R_{inov,k}^L$ when solving for $K_k$.

### Square Root Information Filter

The calling sequence for `srif` is

```fortran
call srif(t,ldt,ainv,lda,b,ldb,rinv,ldr,c,ldc,qinv,ldq,x,rinvy,w,n,m,p,wrk,ldw,multab,multrc,withx,tol)
```
with parameters as follows:
* `t`: A two dimensional real array containing the `n`✕`n` upper triangular square root $T_k$ of the inverse of the state covariance matrix $P_{k|k}$. Upon return, it contains the corresponding factor of $P_{k+1|k+1}$.
* `ldt`: The first (integer) dimension of `t` as declared in the calling program.
* `ainv`: The two-dimensional real array containing the inverse $A_k^{-1}$ of the `n`✕`n` state transition matrix of the discrete-time system (if `multab=.false.`).
* `lda`: The first (integer) dimension of `ainv` as declared in the calling program.
* `b`: A two-dimensional real array containing the `n`✕`m` input weight matrix $B_k$ (or its product with `ainv` if `multab=.true.`) of the discrete-time system.
* `ldb`: The first (integer) dimension of `b` as declared in the calling program.
* `rinv`: A two-dimensional real array containing the `p`✕`p` upper triangular Cholesky factor $R_{k+1}^U$ of the inverse of the measurement noise covariance (if `multrc=.false.`).
* `ldr` The first (integer) dimension of `rinv` as declared in the calling program.
* `c`: A two-dimensional real array containing the `p`✕`n` output weight matrix $C_{k+1}$ (or its product with `rinv` if `multrc=.true.`).
* `ldc`: The first (integer) dimension of `c` as declared in the calling program.
* `qinv`: A two-dimensional real array containing the `m`✕`m` upper triangular Cholesky factor $Q_{k}^U$ of the inverse of the process noise covariance.
* `ldq`: The first (integer) dimension of `q` as declared in the calling program.
* `x`: The real vector containing the estimated state $`\hat{x}_{k|k}`$ upon input and containing the estimated state $`\hat{x}_{k+1|k+1}`$ upon return (if `withx=.true.`). 
* `rinvy`: The real product of `rinv` with the measured output $y_{k+1}$.
* `w`: The real mean value of the state process nosie $\bar{w}_k$
* `n`: Dimension (integer) of the state (actual value used in `srif`).
* `m`: Dimension (integer) of the input (actual value used in `srif`).
* `p`: Dimension (integer) of the output (actual value used in `srif`).
* `wrk`: A real working array of dimensions `ldw`x`nw` containing the postarray upon return. `nw` must be at least `(m+n+1)`.
* `ldw`: The first (integer) dimension of `wrk`, which must be at least `(m+n+p)`.
* `multab`: Logical variable indicating how input matrices `ainv` and `b` are transferred. If `multab =.true.`, then the product `ainv*b` is transferred via `b`, and `ainv` is not used. If `multab =.false.`, then `ainv` and `b` are transferred via their respective parameters.
* `multrc`: Logical variable indicating how input matrices `rinv` and `c` are transferred. If `multrc =.true.`, then the product `rinv*c` is transferred via `c`, and `rinv` is not used. If `multrc =.false.`, then `rinv` and `c` are transferred via their respective parameters.
* `withx`: Logical variable indicating whether or not `x` is requested. If a nearly singular matrix is encountered during the calculation of `x`, `withx` is set to `.false.` upon return.
* `tol`: Real indicating the tolerance for the reciprocal condition number of the computed $T_{k+1}$ when solving for $`\hat{x}_{k+1|k+1}`$.

### Condensed forms

For _time invariant_ systems, considerable savings can be obtained using so-called _condensed forms_.

#### Observer Hessenberg form

The calling sequence for `obhess` is
```fortran
call obhess(a,lda,n,c,ldc,p,u,ldu,withu,upper)
```
with parameters as follows (those preceded by an asterisk are altered by the routine):

* `*a`: A two-dimensional real array containing the `n`✕`n` state transition matrix $A$ to be transformed.
* `lda`: The first (integer) dimension of `a` as declared in the calling program.
* `n`: Dimension (integer) of the state (actual value used in `obhess`).
* `*c`: A two-dimensional real array containing the `p`✕`n` input matrix $C$ to be transformed.
* `ldc`: The first (integer) dimension of `c` as declared in the calling program.
* `p`: Dimension (integer) of the input (actual value used in `obhess`).
* `*u`: A two-dimensional real array which, upon return, is multiplied (if `withu=.true.`) with the `n`✕`n` state space transformation $U$ reducing the given pair to observer Hessenberg form.
* `ldu`: The first (integer) dimension of `u` as declared in the calling program. 
* `withu`: Logical variable equal to `.true.` when `u` is requested and `.false.` otherwise.
* `upper`: Logical variable equal to `.true.` when upper Hessenberg is requested and `.false.` when lower Hessenberg is requested.


#### Controller Hessenberg form

The calling sequence for `cohess` is
```fortran
call cohess(a,lda,n,b,ldb,m,u,ldu,withu,upper)
```
with parameters as follows (those preceded by an asterisk are altered by the routine):

* `*a`: A two-dimensional real array containing the `n`✕`n` state transition matrix $A$ to be transformed.
* `lda`: The first (integer) dimension of `a` as declared in the calling program.
* `n`: Dimension (integer) of the state (actual value used in `cohess`).
* `*b`: A two-dimensional real array containing the `n`✕`m` input matrix $B$ to be transformed.
* `ldb`: The first (integer) dimension of `b` as declared in the calling program.
* `m`: Dimension (integer) of the input (actual value used in `cohess`).
* `*u`: A two-dimensional real array which, upon return, is multiplied (if `withu=.true.`) with the `n`✕`n` state space transformation $U$ reducing the given pair to controller Hessenberg form.
* `ldu`: The first (integer) dimension of `u` as declared in the calling program. 
* `withu`: Logical variable equal to `.true.` when `u` is requested and `.false.` otherwise.
* `upper`: Logical variable equal to `.true.` when upper Hessenberg is requested and `.false.` when lower Hessenberg is requested.


## Differences from the original software

* In routines `srcf`, `srcfob`, `srif`, and `srifco` the calls to the LINPACK routine `dtrco` for calculating the reciprocal condition number have been replaced with calls to the LAPACK routine `dtrcon`.