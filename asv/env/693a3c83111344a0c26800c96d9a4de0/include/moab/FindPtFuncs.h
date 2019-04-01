#ifndef FINDPTFUNCS_H
#define FINDPTFUNCS_H

#include "float.h"
#include "math.h"

/*======================================================
/ from types.h
/======================================================
*/

#define MOAB_POLY_EPS (128*DBL_EPSILON)
#define MOAB_POLY_PI  3.1415926535897932384626433832795028841971693993751058209749445923

#define mbsqrt sqrt
#define mbabs fabs
#define mbcos cos
#define mbsin sin
#define mbfloor floor
#define mbceil ceil

/* integer type to use for everything */
#if   defined(USE_LONG)
#  define INTEGER long
#elif defined(USE_LONG_LONG)
#  define INTEGER long long
#elif defined(USE_SHORT)
#  define INTEGER short
#else
#  define INTEGER int
#endif

/* when defined, use the given type for global indices instead of INTEGER */
#if   defined(USE_GLOBAL_LONG_LONG)
#  define GLOBAL_INT long long
#elif defined(USE_GLOBAL_LONG)
#  define GLOBAL_INT long
#else
#  define GLOBAL_INT long
#endif

/* floating point type to use for everything */
#if   defined(USE_FLOAT)
   typedef float realType;
#elif defined(USE_LONG_DOUBLE)
   typedef long double realType;
#else
   typedef double realType;
#endif

/* apparently uint and ulong can be defined already in standard headers */
#ifndef uint
typedef   signed INTEGER sint;
#endif

#ifndef uint
typedef unsigned INTEGER uint;
#endif
#undef INTEGER

#ifndef slong
#ifdef GLOBAL_INT
  typedef   signed GLOBAL_INT slong;
#else
  typedef sint slong;
#endif
#endif

#ifndef ulong
#ifdef GLOBAL_INT
  typedef unsigned GLOBAL_INT ulong;
#else
  typedef uint ulong;
#endif
#endif

/*======================================================
/ from poly.h
/======================================================
*/

/* 
  For brevity's sake, some names have been shortened
  Quadrature rules
    Gauss   -> Gauss-Legendre quadrature (open)
    Lobatto -> Gauss-Lobatto-Legendre quadrature (closed at both ends)
  Polynomial bases
    Legendre -> Legendre basis
    Gauss    -> Lagrangian basis using Gauss   quadrature nodes
    Lobatto  -> Lagrangian basis using Lobatto quadrature nodes
*/

/*--------------------------------------------------------------------------
   Legendre Polynomial Matrix Computation
   (compute P_i(x_j) for i = 0, ..., n and a given set of x)
  --------------------------------------------------------------------------*/

/* precondition: n >= 1
   inner index is x index (0 ... m-1);
   outer index is Legendre polynomial number (0 ... n)
 */
void legendre_matrix(const realType *x, int m, realType *P, int n);

/* precondition: n >= 1
   inner index is Legendre polynomial number (0 ... n)
   outer index is x index (0 ... m-1);
 */
void legendre_matrix_t(const realType *x, int m, realType *P, int n);

/* precondition: n >= 1
   compute P_i(x) with i = 0 ... n
 */
void legendre_row(realType x, realType *P, int n);


/*--------------------------------------------------------------------------
   Quadrature Nodes and Weights Calculation
   
   call the _nodes function before calling the _weights function
  --------------------------------------------------------------------------*/

void gauss_nodes(realType *z, int n);   /* n nodes (order = 2n-1) */
void lobatto_nodes(realType *z, int n); /* n nodes (order = 2n-3) */

void gauss_weights(const realType *z, realType *w, int n);
void lobatto_weights(const realType *z, realType *w, int n);

/*--------------------------------------------------------------------------
   Lagrangian to Legendre Change-of-basis Matrix
  --------------------------------------------------------------------------*/

/* precondition: n >= 2
   given the Gauss quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Gauss basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_gauss(j)

   computes J   = .5 (2i+1) w  P (z )
             ij              j  i  j
             
   in column major format (inner index is i, the Legendre index)
 */
void gauss_to_legendre(const realType *z, const realType *w, int n, realType *J);

/* precondition: n >= 2
   same as above, but
   in row major format (inner index is j, the Gauss index)
 */
void gauss_to_legendre_t(const realType *z, const realType *w, int n, realType *J);

/* precondition: n >= 3
   given the Lobatto quadrature rule (z,w,n), compute the square matrix J
   for transforming from the Lobatto basis to the Legendre basis:
   
      u_legendre(i) = sum_j J(i,j) u_lobatto(j)

   in column major format (inner index is i, the Legendre index)
 */
void lobatto_to_legendre(const realType *z, const realType *w, int n, realType *J);

/*--------------------------------------------------------------------------
   Lagrangian basis function evaluation
  --------------------------------------------------------------------------*/

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions at all points x
   
   inner index of output J is the basis function index (row-major format)
   provide work array with space for 4*n doubles
 */
void lagrange_weights(const realType *z, unsigned n,
                      const realType *x, unsigned m,
                      realType *J, realType *work);

/* given the Lagrangian nodes (z,n) and evaluation points (x,m)
   evaluate all Lagrangian basis functions and their derivatives
   
   inner index of outputs J,D is the basis function index (row-major format)
   provide work array with space for 6*n doubles
 */
void lagrange_weights_deriv(const realType *z, unsigned n,
                            const realType *x, unsigned m,
                            realType *J, realType *D, realType *work);

/*--------------------------------------------------------------------------
   Speedy Lagrangian Interpolation
   
   Usage:
   
     lagrange_data p;
     lagrange_setup(&p,z,n);    *  setup for nodes z[0 ... n-1] *
     
     the weights
       p->J [0 ... n-1]     interpolation weights
       p->D [0 ... n-1]     1st derivative weights
       p->D2[0 ... n-1]     2nd derivative weights
     are computed for a given x with:
       lagrange_0(p,x);  *  compute p->J *
       lagrange_1(p,x);  *  compute p->J, p->D *
       lagrange_2(p,x);  *  compute p->J, p->D, p->D2 *
       lagrange_2u(p);   *  compute p->D2 after call of lagrange_1(p,x); *
     These functions use the z array supplied to setup
       (that pointer should not be freed between calls)
     Weights for x=z[0] and x=z[n-1] are computed during setup; access as:
       p->J_z0, etc. and p->J_zn, etc.

     lagrange_free(&p);  *  deallocate memory allocated by setup *
  --------------------------------------------------------------------------*/

typedef struct {
  unsigned n;                /* number of Lagrange nodes            */
  const realType *z;             /* Lagrange nodes (user-supplied)      */
  realType *J, *D, *D2;          /* weights for 0th,1st,2nd derivatives */
  realType *J_z0, *D_z0, *D2_z0; /* ditto at z[0]   (computed at setup) */
  realType *J_zn, *D_zn, *D2_zn; /* ditto at z[n-1] (computed at setup) */
  realType *w, *d, *u0, *v0, *u1, *v1, *u2, *v2; /* work data           */
} lagrange_data;

void lagrange_setup(lagrange_data *p, const realType *z, unsigned n);
void lagrange_free(lagrange_data *p);

void lagrange_0(lagrange_data *p, realType x) ;
void lagrange_1(lagrange_data *p, realType x) ;
void lagrange_2(lagrange_data *p, realType x) ;
void lagrange_2u(lagrange_data *p) ;

/*======================================================
/ from tensor.h
/======================================================
*/

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application
   
   the 3d case:
   tensor_f3(R,mr,nr, S,ms,ns, T,mt,nt, u,v, work1,work2)
     gives v = [ R (x) S (x) T ] u
     where R is mr x nr, S is ms x ns, T is mt x nt,
       each in row- or column-major format according to f := r | c
     u is nr x ns x nt in column-major format (inner index is r)
     v is mr x ms x mt in column-major format (inner index is r)
  --------------------------------------------------------------------------*/

void tensor_c1(const realType *R, unsigned mr, unsigned nr, 
               const realType *u, realType *v);
void tensor_r1(const realType *R, unsigned mr, unsigned nr, 
               const realType *u, realType *v);

/* work holds mr*ns realTypes */
void tensor_c2(const realType *R, unsigned mr, unsigned nr,
               const realType *S, unsigned ms, unsigned ns,
               const realType *u, realType *v, realType *work);
void tensor_r2(const realType *R, unsigned mr, unsigned nr,
               const realType *S, unsigned ms, unsigned ns,
               const realType *u, realType *v, realType *work);

/* work1 holds mr*ns*nt realTypes,
   work2 holds mr*ms*nt realTypes */
void tensor_c3(const realType *R, unsigned mr, unsigned nr,
               const realType *S, unsigned ms, unsigned ns,
               const realType *T, unsigned mt, unsigned nt,
               const realType *u, realType *v, realType *work1, realType *work2);
void tensor_r3(const realType *R, unsigned mr, unsigned nr,
               const realType *S, unsigned ms, unsigned ns,
               const realType *T, unsigned mt, unsigned nt,
               const realType *u, realType *v, realType *work1, realType *work2);

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors (for Interpolation)
   
   the 3d case:
   v = tensor_i3(Jr,nr, Js,ns, Jt,nt, u, work)
   same effect as tensor_r3(Jr,1,nr, Js,1,ns, Jt,1,nt, u,&v, work1,work2):
     gives v = [ Jr (x) Js (x) Jt ] u
     where Jr, Js, Jt are row vectors (interpolation weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar
  --------------------------------------------------------------------------*/

realType tensor_i1(const realType *Jr, unsigned nr, const realType *u);

/* work holds ns realTypes */
realType tensor_i2(const realType *Jr, unsigned nr,
               const realType *Js, unsigned ns,
               const realType *u, realType *work);

/* work holds ns*nt + nt realTypes */
realType tensor_i3(const realType *Jr, unsigned nr,
               const realType *Js, unsigned ns,
               const realType *Jt, unsigned nt,
               const realType *u, realType *work);

/*--------------------------------------------------------------------------
   1-,2-,3-d Tensor Application of Row Vectors
             for simultaneous Interpolation and Gradient computation
   
   the 3d case:
   v = tensor_ig3(Jr,Dr,nr, Js,Ds,ns, Jt,Dt,nt, u,g, work)
     gives v   = [ Jr (x) Js (x) Jt ] u
           g_0 = [ Dr (x) Js (x) Jt ] u
           g_1 = [ Jr (x) Ds (x) Jt ] u
           g_2 = [ Jr (x) Js (x) Dt ] u
     where Jr,Dr,Js,Ds,Jt,Dt are row vectors
       (interpolation & derivative weights)
     u is nr x ns x nt in column-major format (inner index is r)
     v is a scalar, g is an array of 3 realTypes
  --------------------------------------------------------------------------*/

realType tensor_ig1(const realType *Jr, const realType *Dr, unsigned nr,
                const realType *u, realType *g);

/* work holds 2*ns realTypes */
realType tensor_ig2(const realType *Jr, const realType *Dr, unsigned nr,
                const realType *Js, const realType *Ds, unsigned ns,
                const realType *u, realType *g, realType *work);

/* work holds 2*ns*nt + 3*ns realTypes */
realType tensor_ig3(const realType *Jr, const realType *Dr, unsigned nr,
                const realType *Js, const realType *Ds, unsigned ns,
                const realType *Jt, const realType *Dt, unsigned nt,
                const realType *u, realType *g, realType *work);

/*======================================================
/ from findpt.h
/======================================================
*/

typedef struct {
  realType x[2], A[4], axis_bnd[4];
} obbox_2;

typedef struct {
  realType x[3], A[9], axis_bnd[6];
} obbox_3;

typedef struct {
  unsigned hash_n;
  realType bnd[4]; /* bounds for all elements */
  realType fac[2]; /* fac[i] = hash_n / (bnd[2*i+1]-bnd[2*i]) */
  obbox_2 *obb; /* obb[nel] -- bounding box info for each element */
  uint *offset; /* hash table -- for cell i,j:
                         uint index = j*hash_n+i,
                                  b = offset[index  ],
                                  e = offset[index+1];
                         elements in cell are
                           offset[b], offset[b+1], ..., offset[e-1] */
  unsigned max; /* maximum # of elements in any cell */
} hash_data_2;

typedef struct {
  unsigned hash_n;
  realType bnd[6]; /* bounds for all elements */
  realType fac[3]; /* fac[i] = hash_n / (bnd[2*i+1]-bnd[2*i]) */
  obbox_3 *obb; /* obb[nel] -- bounding box info for each element */
  uint *offset; /* hash table -- for cell i,j,k:
                         uint index = (k*hash_n+j)*hash_n+i,
                                  b = offset[index  ],
                                  e = offset[index+1];
                         elements in cell are
                           offset[b], offset[b+1], ..., offset[e-1] */
  unsigned max; /* maximum # of elements in any cell */
} hash_data_3;

typedef struct {
  uint el;
  realType r[3];
  realType dist;
} findpt_listel;

typedef struct {
  unsigned constraints;
  unsigned de, d1;
  realType *x[2], *fd1[2];
} opt_edge_data_2;

typedef struct {
  unsigned constraints;
  realType x[2], jac[4];
} opt_point_data_2;

typedef struct {
  lagrange_data *ld;
  unsigned size[3];
  const realType *elx[2];
  opt_edge_data_2 ed;
  opt_point_data_2 pd;
  realType *work;
  realType x[2], jac[4];
} opt_data_2;

typedef struct {
  const realType *xw[2];   /* geometry data */
  realType *z[2];          /* lobatto nodes */
  lagrange_data ld[2]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  hash_data_2 *hash;   /* geometric hashing data */
  findpt_listel *list, **sorted, **end;
                                        /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  opt_data_2 *od; /* data for the optimization algorithm */
  realType *od_work;
} findpt_data_2;

typedef struct {
  unsigned constraints;
  unsigned dn, d1, d2;
  realType *x[3], *fdn[3];
} opt_face_data_3;

typedef struct {
  unsigned constraints;
  unsigned de, d1, d2;
  realType *x[3], *fd1[3], *fd2[3];
} opt_edge_data_3;

typedef struct {
  unsigned constraints;
  realType x[3], jac[9];
} opt_point_data_3;

typedef struct {
  lagrange_data *ld;
  unsigned size[4];
  const realType *elx[3];
  opt_face_data_3 fd;
  opt_edge_data_3 ed;
  opt_point_data_3 pd;
  realType *work;
  realType x[3], jac[9];
} opt_data_3;

typedef struct {
  const realType *xw[3];   /* geometry data */
  realType *z[3];          /* lobatto nodes */
  lagrange_data ld[3]; /* interpolation, derivative weights & data */
  unsigned nptel;      /* nodes per element */
  hash_data_3 *hash;   /* geometric hashing data */
  findpt_listel *list, **sorted, **end;
                                        /* pre-allocated list of elements to
                                           check (found by hashing), and
                                           pre-allocated list of pointers into
                                           the first list for sorting */
  opt_data_3 *od; /* data for the optimization algorithm */
  realType *od_work;
} findpt_data_3;

findpt_data_2 *findpt_setup_2(
          const realType *const xw[2], const unsigned n[2], uint nel,
          uint max_hash_size, realType bbox_tol);
findpt_data_3 *findpt_setup_3(
          const realType *const xw[3], const unsigned n[3], uint nel,
          uint max_hash_size, realType bbox_tol);

void findpt_free_2(findpt_data_2 *p);
void findpt_free_3(findpt_data_3 *p);

const realType *findpt_allbnd_2(const findpt_data_2 *p);
const realType *findpt_allbnd_3(const findpt_data_3 *p);

typedef int (*findpt_func)(void *, const realType *, int, uint *, realType *, realType *);
int findpt_2(findpt_data_2 *p, const realType x[2], int guess,
             uint *el, realType r[2], realType *dist);
int findpt_3(findpt_data_3 *p, const realType x[3], int guess,
             uint *el, realType r[3], realType *dist);

void findpt_weights_2(findpt_data_2 *p, const realType r[2]);

void findpt_weights_3(findpt_data_3 *p, const realType r[3]);

double findpt_eval_2(findpt_data_2 *p, const realType *u);

double findpt_eval_3(findpt_data_3 *p, const realType *u);

/*======================================================
/ from extrafindpt.h
/======================================================
*/

void opt_alloc_3(opt_data_3 *p, lagrange_data *ld);
void opt_free_3(opt_data_3 *p);
double opt_findpt_3(opt_data_3 *p, const realType *const elx[3],
                           const realType xstar[3], realType r[3], unsigned *constr);
void opt_vol_set_intp_3(opt_data_3 *p, const realType r[3]);

/* for 2d spectralQuad */
/*--------------------------------------------------------------------------

   2 - D

  --------------------------------------------------------------------------*/

void opt_alloc_2(opt_data_2 *p, lagrange_data *ld);
void opt_free_2(opt_data_2 *p);
double opt_findpt_2(opt_data_2 *p, const realType *const elx[2],
                           const realType xstar[2], realType r[2], unsigned *constr);

extern const unsigned opt_no_constraints_2;
extern const unsigned opt_no_constraints_3;

/*======================================================
/ from errmem.h
/======================================================
*/

/* requires:
     <stdlib.h> for malloc, calloc, realloc, free
*/

/*--------------------------------------------------------------------------
   Error Reporting
   Memory Allocation Wrappers to Catch Out-of-memory
  --------------------------------------------------------------------------*/

/* #include "malloc.h" */
#include <stdlib.h>

#ifdef __GNUC__
void fail(const char *fmt, ...) __attribute__ ((noreturn));
#define MAYBE_UNUSED __attribute__ ((unused))
#else
void fail(const char *fmt, ...);
#define MAYBE_UNUSED
#endif

#if 0
{}
#endif

static void *smalloc(size_t size, const char *file) MAYBE_UNUSED;
static void *smalloc(size_t size, const char *file)
{
  void *res = malloc(size);
  if(!res && size) fail("%s: allocation of %d bytes failed\n",file,(int)size);
  return res;
}

static void *scalloc(size_t nmemb, size_t size, const char *file) MAYBE_UNUSED;
static void *scalloc(size_t nmemb, size_t size, const char *file)
{
  void *res = calloc(nmemb, size);
  if(!res && nmemb)
    fail("%s: allocation of %d bytes failed\n",file,(int)size*nmemb);
  return res;
}

static void *srealloc(void *ptr, size_t size, const char *file) MAYBE_UNUSED;
static void *srealloc(void *ptr, size_t size, const char *file)
{
  void *res = realloc(ptr, size);
  if(!res && size) fail("%s: allocation of %d bytes failed\n",file,(int)size);
  return res;
}

#define tmalloc(type, count) \
  ((type*) smalloc((count)*sizeof(type),__FILE__) )
#define tcalloc(type, count) \
  ((type*) scalloc((count),sizeof(type),__FILE__) )
#define trealloc(type, ptr, count) \
  ((type*) srealloc((ptr),(count)*sizeof(type),__FILE__) )
/*
typedef struct { size_t size; void *ptr; } buffer;

static void buffer_init_(buffer *b, size_t size, const char *file) MAYBE_UNUSED;
static void buffer_init_(buffer *b, size_t size, const char *file)
{
  b->size=size, b->ptr=smalloc(size,file);
}
static void buffer_reserve_(buffer *b, size_t min, const char *file)
  MAYBE_UNUSED;
static void buffer_reserve_(buffer *b, size_t min, const char *file)
{
  size_t size = b->size;
  if(size<min) {
    size+=size/2+1;
    if(size<min) size=min;
    b->ptr=srealloc(b->ptr,size,file);
  }
}
static void buffer_free(buffer *b) MAYBE_UNUSED;
static void buffer_free(buffer *b) { free(b->ptr); }

#define buffer_init(b,size) buffer_init_(b,size,__FILE__)
#define buffer_reserve(b,min) buffer_reserve_(b,min,__FILE__)
*/


/*======================================================
/ from minmax.h
/======================================================
*/

/*--------------------------------------------------------------------------
   Min, Max, Norm
  --------------------------------------------------------------------------*/

#ifdef __GNUC__
#define MAYBE_UNUSED __attribute__ ((unused))
#else
#define MAYBE_UNUSED
#endif

#define DECLMINMAX(type, prefix) \
  static type prefix##min_2(type a, type b) MAYBE_UNUSED;               \
  static type prefix##min_2(type a, type b) {                           \
    return b<a?b:a;                                                     \
  }                                                                     \
  static type prefix##max_2(type a, type b) MAYBE_UNUSED;               \
  static type prefix##max_2(type a, type b) {                           \
    return a>b?a:b;                                                     \
  }                                                                     \
  static void prefix##minmax_2(type *min, type *max, type a,            \
                               type b) MAYBE_UNUSED;                    \
  static void prefix##minmax_2(type *min, type *max, type a, type b) {  \
    if(b<a) *min=b, *max=a;                                             \
    else *min=a, *max=b;                                                \
  }                                                                     \
  static type prefix##min_3(type a, type b, type c) MAYBE_UNUSED;       \
  static type prefix##min_3(type a, type b, type c) {                   \
    return b<a?(c<b?c:b):(c<a?c:a);                                     \
  }                                                                     \
  static type prefix##max_3(type a, type b, type c) MAYBE_UNUSED;       \
  static type prefix##max_3(type a, type b, type c) {                   \
    return a>b?(a>c?a:c):(b>c?b:c);                                     \
  }                                                                     \
  static void prefix##minmax_3(type *min, type *max, type a, type b,    \
                               type c) MAYBE_UNUSED;                    \
  static void prefix##minmax_3(type *min, type *max, type a, type b,    \
                               type c) {                                \
    if(b<a) *min=prefix##min_2(b,c), *max=prefix##max_2(a,c);           \
    else    *min=prefix##min_2(a,c), *max=prefix##max_2(b,c);           \
  }

DECLMINMAX(int, i)
DECLMINMAX(unsigned, u)
DECLMINMAX(realType, r)
#undef DECLMINMAX

static realType r1norm_1(realType a) MAYBE_UNUSED;
static realType r1norm_1(realType a) {
  return mbabs(a);
}
static realType r1norm_2(realType a, realType b) MAYBE_UNUSED;
static realType r1norm_2(realType a, realType b) {
  return mbabs(a)+mbabs(b);
}
static realType r1norm_3(realType a, realType b, realType c) MAYBE_UNUSED;
static realType r1norm_3(realType a, realType b, realType c) {
  return mbabs(a)+mbabs(b)+mbabs(c);
}
static realType r2norm_1(realType a) MAYBE_UNUSED;
static realType r2norm_1(realType a) {
  return mbsqrt(a*a);
}
static realType r2norm_2(realType a, realType b) MAYBE_UNUSED;
static realType r2norm_2(realType a, realType b) {
  return mbsqrt(a*a+b*b);
}
static realType r2norm_3(realType a, realType b, realType c) MAYBE_UNUSED;
static realType r2norm_3(realType a, realType b, realType c) {
  return mbsqrt(a*a+b*b+c*c);
}

#endif


