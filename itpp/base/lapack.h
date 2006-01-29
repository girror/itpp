/*!
 * \file
 * \brief Lapack header functions. For internal use only.
 * \author Tony Ottosson
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#ifndef LAPACK_H
#define LAPACK_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#if defined(HAVE_LAPACK) || defined(HAVE_MKL)

#ifndef DOXYGEN_SHOULD_SKIP_THIS

extern "C" {

// Fix for MKL Windows version so that naming is consistent with the 5.x and 7.x MKL LAPACK libraries
#ifdef HAVE_MKL

#define dgetrf_ dgetrf
#define zgetrf_ zgetrf

#define dgetri_ dgetri
#define zgetri_ zgetri

#define dgesvd_ dgesvd
#define zgesvd_ zgesvd

#define dsyev_ dsyev
#define zheev_ zheev

#define dgeev_ dgeev
#define zgeev_ zgeev

#define dpotrf_ dpotrf
#define zpotrf_ zpotrf

#define dgeqrf_ dgeqrf
#define zgeqrf_ zgeqrf

#define dgeqp3_ dgeqp3
#define zgeqp3_ zgeqp3

#define dorgqr_ dorgqr
#define zungqr_ zungqr
 
#define dormqr_ dormqr
#define zunmqr_ zunmqr

#define dgesv_ dgesv
#define zgesv_ zgesv

#define dposv_ dposv
#define zposv_ zposv

#define dtrtrs_ dtrtrs
#define ztrtrs_ ztrtrs

#define dgels_ dgels
#define zgels_ zgels

#endif // #ifdef HAVE_MKL


  // Exists in ATLAS
/* LU factorization
   a is of size m*n and with lda rows.
   ipiv is the permutation vector of rows. Row i should be replaced by row ipiv(i).
   info=0 if OK. info=-i if ith value is illegal. info=i factorization OK but the system is singular if solved.
*/
void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
void zgetrf_(int *m, int *n, std::complex<double> *a, int *lda, int *ipiv, int *info);

  // In ATLAS
/* Inverting a matrix of an LU-factored general matrix (first call xGETRF)
   a is of square size n*n with lda rows containing the factorization as returned by xGETRF
   ipiv is vector as returned by xGETRF
   lwork >= n
   output: a is overwritten by the inverse
   info=0 if OK. info=-i if ith parameter is illegal. info=i the ith diagonal element = 0 and U is singular.
*/
void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void zgetri_(int *n, std::complex<double> *a, int *lda, int *ipiv, std::complex<double> *work, int *lwork, int *info);

  
/* SVD of a general rectangular matrix A = U S V^H
   a is of size m*n and with lda rows.
   Output: s with sorted singular values (vector)
           u, and vt (for U and V^H). U is m*m, and V^H is n*n
   jobu='A','S','O','N'. Different versions. 'A' = all columns of U calculated and returned in u.
   jobvt='A','S','O','N'. Different versions. 'A' = all columns of V^H calculated and returned in vt.
   ldu = no rows in U
   ldvt = no rows in V^H
   info = 0 successful, =-i ith parameter is illegal, =i did not converge

   work is a workspace vector of size lwork.
   lwork >= max(3*min(m,n)+max(m,n), 5*min(m,n)) for double
   lwork >= 2*min(m,n)+max(m,n) for std::complex<double>
   Good performance. Make lwork larger!
   rwork is a workspace array for complex version. Size max(1, 5*min(m,n)).
 
*/
void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, int *info);
void zgesvd_(char *jobu, char *jobvt, int *m, int *n, std::complex<double> *a, int *lda, double *s, std::complex<double> *u, int *ldu, std::complex<double> *vt, int *ldvt, std::complex<double> *work, int *lwork, double *rwork, int *info);

/* Eigenvalues and eigenvectors of a symmetric/hermitian matrix A */
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);
void zheev_(char *jobz, char *uplo, int *n, std::complex<double> *a, int *lda, double *w, std::complex<double> *work, int *lwork, double *rwork, int *info);


/* Eigenvalues and eigenvectors of a general matrix A */
void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
void zgeev_(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda, std::complex<double> *w, std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr, std::complex<double> *work, int *lwork, double *rwork, int *info);

  // In ATLAS
/* Cholesky factorization */
void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
void zpotrf_(char *uplo, int *n, std::complex<double> *a, int *lda, int *info);

/* QR factorization of a general matrix A  */
void dgeqrf_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
void zgeqrf_(int *m, int *n, std::complex<double> *a, int *lda, std::complex<double> *tau, std::complex<double> *work, int *lwork, int *info);

/* QR factorization of a general matrix A with pivoting */
void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
void zgeqp3_(int *m, int *n, std::complex<double> *a, int *lda, int *jpvt, std::complex<double> *tau, std::complex<double> *work, int *lwork, double *rwork, int *info);

/* Calculation of Q matrix from QR-factorization */
void dorgqr_(int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info);
void zungqr_(int *m, int *n, int *k, std::complex<double> *a, int *lda, std::complex<double> *tau, std::complex<double> *work, int *lwork, int *info);
 

/* Multiplies a real matrix by the orthogonal matix Q of the QR factorization formed by dgeqp3_() */
void dormqr_(char *side, char *trans, int *m, int *n, int *k, double *a, int *lda, double *tau, double *c, int *ldc, double *work, int *lwork, int *info);
/* Multiplies a complex matrix by the unitary matix Q of the QR factorization formed by zgeqp3_() */
void zunmqr_(char *side, char *trans, int *m, int *n, int *k, std::complex<double> *a, int *lda, std::complex<double> *tau, std::complex<double> *c, int *ldc, std::complex<double> *work, int *lwork, int *info);

  // In ATLAS
/* Solves a system on linear equations, Ax=b, with a square matrix A, Using LU-factorization */
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
void zgesv_(int *n, int *nrhs, std::complex<double> *a, int *lda, int *ipiv, std::complex<double> *b, int *ldb, int *info);

  // In ATLAS
/* Solves a system on linear equations, Ax=b, with a square symmetric/hermitian positive definite matrix A, Using Cholesky-factorization */
void dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
void zposv_(char *uplo, int *n, int *nrhs, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, int *info);

/* Solves a system of linear equations with a triangular matrix with multiple right-hand sides */
void dtrtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, int *info);
void ztrtrs_(char *uplo, char *trans, char *diag, int *n, int *nrhs, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, int *info);

/* Solves a linear over/underdetermined system using QR or LQ factorization. Assumes a full rank matrix  */
void dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);
void zgels_(char *trans, int *m, int *n, int *nrhs, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, std::complex<double> *work, int *lwork, int *info);

} // extern C

#endif // #ifndef DOXYGEN_SHOULD_SKIP_THIS

#endif // #if defined(HAVE_LAPACK) || defined(HAVE_MKL)

#endif // #ifndef LAPACK_H
