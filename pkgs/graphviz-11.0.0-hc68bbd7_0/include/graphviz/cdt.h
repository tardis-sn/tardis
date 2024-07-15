/**
 * @file
 * @brief container data types API
 * @ingroup public_apis
 *
 * **CDT** manages run-time dictionaries using standard container data types:
 * unordered set/multiset, ordered set/multiset, list, stack, and queue.
 *
 * [man 3 cdt](https://graphviz.org/pdf/cdt.3.pdf)
 *
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*	Public interface for the dictionary library
**
**      Written by Kiem-Phong Vo
*/

#define CDT_VERSION	20050420L

#include <stddef.h>	/* size_t */
#include <string.h>

#ifdef GVDLL
#ifdef EXPORT_CDT
#define CDT_API __declspec(dllexport)
#else
#define CDT_API __declspec(dllimport)
#endif
#endif

#ifndef CDT_API
#define CDT_API /* nothing */
#endif

typedef struct _dtlink_s	Dtlink_t;
typedef struct _dthold_s	Dthold_t;
typedef struct _dtdisc_s	Dtdisc_t;
typedef struct _dtmethod_s	Dtmethod_t;
typedef struct _dtdata_s	Dtdata_t;
typedef struct _dt_s		Dt_t;
typedef struct _dt_s		Dict_t;	/* for libdict compatibility */
typedef struct _dtstat_s	Dtstat_t;
typedef void*			(*Dtsearch_f)(Dt_t*,void*,int);
typedef void* 		(*Dtmake_f)(void*,Dtdisc_t*);
typedef void 			(*Dtfree_f)(void*,Dtdisc_t*);
typedef int			(*Dtcompar_f)(Dt_t*,void*,void*,Dtdisc_t*);

struct _dtlink_s
{	Dtlink_t*	right;	/* right child		*/
	union
	{ unsigned int	_hash;	/* hash value		*/
	  Dtlink_t*	_left;	/* left child		*/
	} hl;
};

/* private structure to hold an object */
struct _dthold_s
{	Dtlink_t	hdr;	/* header		*/
	void*		obj;	/* user object		*/
};

/* method to manipulate dictionary structure */
struct _dtmethod_s
{	Dtsearch_f	searchf; /* search function	*/
	int		type;	/* type of operation	*/
};

/* stuff that may be in shared memory */
struct _dtdata_s
{	int		type;	/* type of dictionary			*/
	Dtlink_t*	here;	/* finger to last search element	*/
	union
	{ Dtlink_t**	_htab;	/* hash table				*/
	  Dtlink_t*	_head;	/* linked list				*/
	} hh;
	int		ntab;	/* number of hash slots			*/
	int		size;	/* number of objects			*/
	int		loop;	/* number of nested loops		*/
};

/* structure to hold methods that manipulate an object */
struct _dtdisc_s
{	int		key;	/* where the key begins in an object	*/
	int		size;	/* key size and type			*/
	int		link;	/* offset to Dtlink_t field		*/
	Dtmake_f	makef;	/* object constructor			*/
	Dtfree_f	freef;	/* object destructor			*/
	Dtcompar_f	comparf;/* to compare two objects		*/
};

#define DTDISC(dc, ky, sz, lk, mkf, frf, cmpf) \
	( (dc)->key = (ky), (dc)->size = (sz), (dc)->link = (lk), \
	  (dc)->makef = (mkf), (dc)->freef = (frf), \
	  (dc)->comparf = (cmpf) )

/* the dictionary structure itself */
struct _dt_s
{	Dtsearch_f	searchf;/* search function			*/
	Dtdisc_t*	disc;	/* method to manipulate objs		*/
	Dtdata_t*	data;	/* sharable data			*/
	Dtmethod_t*	meth;	/* dictionary method			*/
	int		nview;	/* number of parent view dictionaries	*/
	Dt_t*		view;	/* next on viewpath			*/
	Dt_t*		walk;	/* dictionary being walked		*/
	void*		user;	/* for user's usage			*/
};

/* structure to get status of a dictionary */
struct _dtstat_s
{	int	dt_meth;	/* method type				*/
	int	dt_size;	/* number of elements			*/
	size_t dt_n; // number of chains or levels
	size_t dt_max; // max size of a chain or a level
	size_t* dt_count; // counts of chains or levels by size
};

/* supported storage methods */
#define DT_SET		0000001	/* set with unique elements		*/
#define DT_OSET		0000004	/* ordered set (self-adjusting tree)	*/
#define DT_OBAG		0000010	/* ordered multiset			*/
#define DT_QUEUE	0000100	/* queue: insert at top, delete at tail	*/
#define DT_METHODS	0000377	/* all currently supported methods	*/

/* types of search */
#define DT_INSERT	0000001	/* insert object if not found		*/
#define DT_DELETE	0000002	/* delete object if found		*/
#define DT_SEARCH	0000004	/* look for an object			*/
#define DT_NEXT		0000010	/* look for next element		*/
#define DT_PREV		0000020	/* find previous element		*/
#define DT_RENEW	0000040	/* renewing an object			*/
#define DT_CLEAR	0000100	/* clearing all objects			*/
#define DT_FIRST	0000200	/* get first object			*/
#define DT_LAST		0000400	/* get last object			*/
#define DT_MATCH	0001000	/* find object matching key		*/
#define DT_VSEARCH	0002000	/* search using internal representation	*/
#define DT_DETACH	0010000	/* detach an object from the dictionary	*/

CDT_API extern Dtmethod_t* 	Dtset; ///< set with unique elements
CDT_API extern Dtmethod_t* 	Dtoset; ///< ordered set (self-adjusting tree)
CDT_API extern Dtmethod_t* 	Dtobag; ///< ordered multiset
CDT_API extern Dtmethod_t*	Dtqueue; ///< queue: insert at top, delete at tail

CDT_API extern Dtmethod_t*	Dttree;
CDT_API extern Dtmethod_t	_Dttree;
CDT_API extern Dtmethod_t	_Dtqueue;

CDT_API Dt_t*		dtopen(Dtdisc_t*, Dtmethod_t*);
CDT_API int		dtclose(Dt_t*);
CDT_API Dt_t*		dtview(Dt_t*, Dt_t*);
CDT_API Dtdisc_t*	dtdisc(Dt_t *dt, Dtdisc_t*);
CDT_API Dtmethod_t*	dtmethod(Dt_t*, Dtmethod_t*);

CDT_API Dtlink_t*	dtflatten(Dt_t*);
CDT_API Dtlink_t*	dtextract(Dt_t*);
CDT_API int		dtrestore(Dt_t*, Dtlink_t*);

CDT_API int		dtwalk(Dt_t*, int(*)(void*,void*), void*);

CDT_API void*		dtrenew(Dt_t*, void*);

CDT_API int		dtsize(Dt_t*);
CDT_API int		dtstat(Dt_t*, Dtstat_t*, int);
CDT_API unsigned int dtstrhash(void*, int);

/* internal functions for translating among holder, object and key */
#define _DT(dt)		((Dt_t*)(dt))
#define _DTDSC(dc,ky,sz,lk,cmpf) \
			(ky = dc->key, sz = dc->size, lk = dc->link, cmpf = dc->comparf)
#define _DTLNK(o,lk)	((Dtlink_t*)((char*)(o) + lk) )
#define _DTOBJ(e,lk)	(lk < 0 ? ((Dthold_t*)(e))->obj : (void*)((char*)(e) - lk) )
#define _DTKEY(o,ky,sz)	(void*)(sz < 0 ? *((char**)((char*)(o)+ky)) : ((char*)(o)+ky))

#define _DTCMP(dt,k1,k2,dc,cmpf,sz) \
			(cmpf ? (*cmpf)(dt,k1,k2,dc) : \
			 (sz <= 0 ? strcmp(k1,k2) : memcmp(k1,k2,(size_t)sz)) )

#define dtlink(d,e)	(((Dtlink_t*)(e))->right)
#define dtobj(d,e)	_DTOBJ((e), _DT(d)->disc->link)
#define dtfinger(d)	(_DT(d)->data->here ? dtobj((d),_DT(d)->data->here):(void*)(0))

#define dtfirst(d)	(*(_DT(d)->searchf))((d),(void*)(0),DT_FIRST)
#define dtnext(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_NEXT)
#define dtlast(d)	(*(_DT(d)->searchf))((d),(void*)(0),DT_LAST)
#define dtprev(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_PREV)
#define dtsearch(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_SEARCH)
#define dtmatch(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_MATCH)
#define dtinsert(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_INSERT)
#define dtdelete(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_DELETE)
#define dtdetach(d,o)	(*(_DT(d)->searchf))((d),(void*)(o),DT_DETACH)
#define dtclear(d)	(*(_DT(d)->searchf))((d),(void*)(0),DT_CLEAR)

#define DT_PRIME	17109811 /* 2#00000001 00000101 00010011 00110011 */

/**
 * @dir lib/cdt
 * @brief container data types, API cdt.h
 *
 * [man 3 cdt](https://graphviz.org/pdf/cdt.3.pdf)
 *
 */

#ifdef __cplusplus
}
#endif
