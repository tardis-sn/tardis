/**
 * @file
 * @brief abstract graph C library, @ref cgraph_api
 * @ingroup cgraph_api
 *
 * **Libcgraph** supports graph programming by maintaining graphs
 * in memory and reading and writing graph files.
 * Graphs are composed of nodes, edges, and nested subgraphs.
 * These graph objects may be attributed with string name-value pairs
 * and programmer-defined records (see Attributes).
 * All of Libcgraph’s global symbols have the prefix **ag** (case varying).
 * In the following, if a function has a parameter `int createflag` and
 * the object does not exist, the function will create the specified object
 * if `createflag` is non-zero; otherwise, it will return NULL.
 *
 * [man 3 cgraph](https://graphviz.org/pdf/cgraph.3.pdf)
 *
 * @defgroup cgraph_api Cgraph API
 * @brief Abstract graph C library. API cgraph.h
 * @ingroup public_apis
 *
 * [man 3 cgraph](https://graphviz.org/pdf/cgraph.3.pdf)
 *
 * Main types @ref Agraph_t, @ref Agnode_t, @ref Agedge_t.
 * @{
 */

/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: Details at https://graphviz.org
 *************************************************************************/

#pragma once

#include <inttypes.h>
#include "cdt.h"

#ifdef __cplusplus
extern "C" {
#endif

/// @cond

#ifdef GVDLL
#ifdef EXPORT_CGRAPH
#define CGRAPH_API __declspec(dllexport)
#else
#define CGRAPH_API __declspec(dllimport)
#endif
#endif

#ifndef CGRAPH_API
#define CGRAPH_API /* nothing */
#endif

#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (!FALSE)
#endif

/// @endcond

/// @defgroup cgraph_other other
/// @{
typedef uint64_t IDTYPE;

/* forward struct type declarations */
typedef struct Agtag_s Agtag_t;
typedef struct Agobj_s Agobj_t;         ///< generic object header
/// @}
/// @ingroup cgraph_graph
typedef struct Agraph_s Agraph_t;       ///< graph, subgraph (or hyperedge)
/// @ingroup cgraph_node
typedef struct Agnode_s Agnode_t;       ///< node (atom)
/// @ingroup cgraph_edge
typedef struct Agedge_s Agedge_t;       ///< node pair
/// @ingroup cgraph_graph
typedef struct Agdesc_s Agdesc_t;       ///< graph descriptor
/// @addtogroup cgraph_other
/// @{
typedef struct Agiddisc_s Agiddisc_t;   ///< object ID allocator
typedef struct Agiodisc_s Agiodisc_t;   ///< IO services
typedef struct Agdisc_s Agdisc_t;       ///< union of client discipline methods
typedef struct Agdstate_s Agdstate_t;   ///< client state (closures)
typedef struct Agsym_s Agsym_t;         ///< string attribute descriptors
typedef struct Agattr_s Agattr_t;       ///< string attribute container
typedef struct Agcbdisc_s Agcbdisc_t;   ///< client event callbacks
typedef struct Agcbstack_s Agcbstack_t; ///< enclosing state for cbdisc
typedef struct Agclos_s Agclos_t;       ///< common fields for graph/subgs
typedef struct Agdatadict_s Agdatadict_t; ///< set of dictionaries per graph
typedef struct Agedgepair_s Agedgepair_t; ///< the edge object
typedef struct Agsubnode_s Agsubnode_t;
/// @}

/** @addtogroup cgraph_attr
 *  @{
 *
 *  @defgroup cgraph_rec records
 *  @brief These records are attached by client programs dynamically at runtime.
 *  @{
 *
 *  Uninterpreted records may be attached to graphs, subgraphs, nodes,
 *  and edges for efficient operations on values such as marks, weights,
 *  counts, and pointers needed by algorithms.
 *  Application programmers define the fields of these records,
 *  but they must be declared with a common record header @ref Agrec_t.
 *
 *  A unique string ID (stored in @ref Agrec_s.name) must be given
 *  to each record attached to the same object.
 *  Cgraph has functions to create, search for, and delete these records.
 *  The records are maintained in a circular list,
 *  with @ref Agobj_s.data pointing somewhere in the list.
 *  The search function @ref aggetrec has an option to lock this pointer on a given record.
 *  The application must be written so only one such lock is outstanding at a time.
 *
 *  Records are created and managed by Libcgraph.
 *  A programmer must explicitly attach them to the objects in a graph,
 *  either to individual objects one at a time via @ref agbindrec,
 *  or to all the objects of the same class in a graph via @ref aginit.
 *  The `name` argument of a record distinguishes various types of records,
 *  and is programmer defined.
 *
 *  Libcgraph reserves the prefix "_AG_" (in
 *  @ref DataDictName,
 *  @ref AgDataRecName,
 *  @ref DRName).
 *
 *  To allow referencing application-dependent data without function calls or search,
 *  Libcgraph allows setting and locking the list pointer
 *  of a graph, node, or edge on a particular record
 *  (see @ref Agtag_s.mtflock and @ref Agobj_s.data).
 *  This pointer can be obtained with the macro @ref AGDATA(obj).
 *  A cast, generally within a macro or inline function,
 *  is usually applied to convert the list pointer to
 *  an appropriate programmer-defined type (eg. @ref GD_parent).
 *
 *  To control the setting of this pointer,
 *  the `move_to_front` flag may be TRUE or FALSE.
 *  If `move_to_front` is TRUE, the record will be
 *  locked @ref Agtag_s.mtflock at the
 *  head of the list @ref Agobj_s.data,
 *  so it can be accessed directly by @ref AGDATA(obj).
 *
 *  The lock protects the data pointer from being moved.
 *  Function @ref aggetrec reports error when data pointer and lock
 *  are reassigned.
 *
 *  The lock can be released or reset by a call to @ref agdelrec.
 *
 */

typedef struct Agrec_s Agrec_t;
///< generic header of @ref Agattr_s, @ref Agdatadict_s and user records

/// implementation of @ref Agrec_t
struct Agrec_s {
    char *name;
    Agrec_t *next;
    /* following this would be any programmer-defined data */
};
/// @}

/** @brief Object tag for graphs, nodes, and edges.

While there may be several structs
for a given node or edges, there is only one unique ID (per main graph).  */

struct Agtag_s {
    unsigned objtype:2;		/* see literal tags below */
    unsigned mtflock:1;		/* move-to-front lock, see above */
    unsigned attrwf:1;		/* attrs written (parity, write.c) */
    unsigned seq:(sizeof(unsigned) * 8 - 4);	/* sequence no. */
    IDTYPE id;		        /* client  ID */
};

	/* object tags */
#define AGRAPH				0	/* can't exceed 2 bits. see Agtag_t. */
#define AGNODE				1
#define AGOUTEDGE			2
#define AGINEDGE			3	/* (1 << 1) indicates an edge tag.   */
#define AGEDGE 				AGOUTEDGE	/* synonym in object kind args */

/// a generic graph/node/edge header
struct Agobj_s {
    Agtag_t tag;
    Agrec_t *data;
};

#define AGTAG(obj)		(((Agobj_t*)(obj))->tag)
#define AGTYPE(obj)		(AGTAG(obj).objtype)
#define AGID(obj)		(AGTAG(obj).id)
#define AGSEQ(obj)		(AGTAG(obj).seq)
#define AGATTRWF(obj)		(AGTAG(obj).attrwf)
#define AGDATA(obj)		(((Agobj_t*)(obj))->data)
/// @}

/// @addtogroup cgraph_node
/// @{

/** @brief This is the node struct allocated per graph (or subgraph).

It resides in the n_dict of the graph.
The node set is maintained by libdict, but transparently to libgraph callers.
Every node may be given an optional string name at its time of creation,
or it is permissible to pass NIL(char*) for the name. */

struct Agsubnode_s {		/* the node-per-graph-or-subgraph record */
    Dtlink_t seq_link;		/* must be first */
    Dtlink_t id_link;
    Agnode_t *node;		/* the object */
    Dtlink_t *in_id, *out_id;	/* by node/ID for random access */
    Dtlink_t *in_seq, *out_seq;	/* by node/sequence for serial access */
};

struct Agnode_s {
    Agobj_t base;
    Agraph_t *root;
    Agsubnode_t mainsub;	/* embedded for main graph */
};
/// @}

/// @addtogroup cgraph_edge
/// @{
struct Agedge_s {
    Agobj_t base;
    Dtlink_t id_link;		/* main graph only */
    Dtlink_t seq_link;
    Agnode_t *node;		/* the endpoint node */
};

struct Agedgepair_s {
    Agedge_t out, in;
};
/// @}

/// @addtogroup cgraph_graph
/// @{
struct Agdesc_s {		/* graph descriptor */
    unsigned directed:1;	/* if edges are asymmetric */
    unsigned strict:1;		/* if multi-edges forbidden */
    unsigned no_loop:1;		/* if no loops */
    unsigned maingraph:1;	/* if this is the top level graph */
    unsigned no_write:1;	/* if a temporary subgraph */
    unsigned has_attrs:1;	/* if string attr tables should be initialized */
    unsigned has_cmpnd:1;	/* if may contain collapsed nodes */
};
/// @}

/** @defgroup cgraph_disc disciplines
 *  @brief disciplines for external resources needed by libgraph
 *  @{
 */

/// object ID allocator discipline

struct Agiddisc_s {
    void *(*open) (Agraph_t * g, Agdisc_t*);	/* associated with a graph */
    long (*map) (void *state, int objtype, char *str, IDTYPE *id,
		 int createflag);
    long (*alloc) (void *state, int objtype, IDTYPE id);
    void (*free) (void *state, int objtype, IDTYPE id);
    char *(*print) (void *state, int objtype, IDTYPE id);
    void (*close) (void *state);
    void (*idregister) (void *state, int objtype, void *obj);
};

struct Agiodisc_s {
    int (*afread) (void *chan, char *buf, int bufsize);
    int (*putstr) (void *chan, const char *str);
    int (*flush) (void *chan);	/* sync */
    /* error messages? */
};

/// user's discipline

struct Agdisc_s {
    Agiddisc_t *id;
    Agiodisc_t *io;
};

	/* default resource disciplines */

CGRAPH_API extern Agiddisc_t AgIdDisc;
CGRAPH_API extern Agiodisc_t AgIoDisc;

CGRAPH_API extern Agdisc_t AgDefaultDisc;
/// @}

/** @defgroup cgraph_graph graphs
 *  @{
 */
struct Agdstate_s {
    void *id;
    /* IO must be initialized and finalized outside Cgraph,
     * and channels (FILES) are passed as void* arguments. */
};

typedef void (*agobjfn_t) (Agraph_t * g, Agobj_t * obj, void *arg);
typedef void (*agobjupdfn_t) (Agraph_t * g, Agobj_t * obj, void *arg,
			      Agsym_t * sym);

struct Agcbdisc_s {
    struct {
	agobjfn_t ins;
	agobjupdfn_t mod;
	agobjfn_t del;
    } graph, node, edge;
};

/// object event callbacks

struct Agcbstack_s {
    Agcbdisc_t *f;		/* methods */
    void *state;		/* closure */
    Agcbstack_t *prev;		/* kept in a stack, unlike other disciplines */
};

struct Agclos_s {
    Agdisc_t disc;		/* resource discipline functions */
    Agdstate_t state;		/* resource closures */
    Dict_t *strdict;		/* shared string dict */
    uint64_t seq[3];	/* local object sequence number counter */
    Agcbstack_t *cb;		/* user and system callback function stacks */
    Dict_t *lookup_by_name[3];
    Dict_t *lookup_by_id[3];
};

struct Agraph_s {
    Agobj_t base;
    Agdesc_t desc;
    Dtlink_t link;
    Dict_t *n_seq;		/* the node set in sequence */
    Dict_t *n_id;		/* the node set indexed by ID */
    Dict_t *e_seq, *e_id;	/* holders for edge sets */
    Dict_t *g_dict;		/* subgraphs - descendants */
    Agraph_t *parent, *root;	/* subgraphs - ancestors */
    Agclos_t *clos;		/* shared resources */
};

CGRAPH_API void agpushdisc(Agraph_t * g, Agcbdisc_t * disc, void *state);
CGRAPH_API int agpopdisc(Agraph_t * g, Agcbdisc_t * disc);

/* graphs */
CGRAPH_API Agraph_t *agopen(char *name, Agdesc_t desc, Agdisc_t * disc);
CGRAPH_API int agclose(Agraph_t * g);
CGRAPH_API Agraph_t *agread(void *chan, Agdisc_t * disc);
CGRAPH_API Agraph_t *agmemread(const char *cp);
CGRAPH_API Agraph_t *agmemconcat(Agraph_t *g, const char *cp);
CGRAPH_API void agreadline(int);
CGRAPH_API void agsetfile(const char *);
CGRAPH_API Agraph_t *agconcat(Agraph_t * g, void *chan, Agdisc_t * disc);
CGRAPH_API int agwrite(Agraph_t * g, void *chan);
CGRAPH_API int agisdirected(Agraph_t * g);
CGRAPH_API int agisundirected(Agraph_t * g);
CGRAPH_API int agisstrict(Agraph_t * g);
CGRAPH_API int agissimple(Agraph_t * g);
/// @}

/// @defgroup cgraph_node nodes
/// @{
CGRAPH_API Agnode_t *agnode(Agraph_t * g, char *name, int createflag);
CGRAPH_API Agnode_t *agidnode(Agraph_t * g, IDTYPE id, int createflag);
CGRAPH_API Agnode_t *agsubnode(Agraph_t * g, Agnode_t * n, int createflag);
CGRAPH_API Agnode_t *agfstnode(Agraph_t * g);
CGRAPH_API Agnode_t *agnxtnode(Agraph_t * g, Agnode_t * n);
CGRAPH_API Agnode_t *aglstnode(Agraph_t * g);
CGRAPH_API Agnode_t *agprvnode(Agraph_t * g, Agnode_t * n);

CGRAPH_API Agsubnode_t *agsubrep(Agraph_t * g, Agnode_t * n);
CGRAPH_API int agnodebefore(Agnode_t *u, Agnode_t *v); /* we have no shame */
/// @}

/** @defgroup cgraph_edge edges
 *
 * An abstract edge has two endpoint nodes called tail and head
 * where all outedges of the same node have it as the tail
 * value and similarly all inedges have it as the head.
 * In an undirected graph, head and tail are interchangeable.
 * If a graph has multi-edges between the same pair of nodes,
 * the edge's string name behaves as a secondary key.
 *
 * Note that an abstract edge has two distinct concrete
 * representations: as an in-edge and as an out-edge.
 * In particular, the pointer as an out-edge is different
 * from the pointer as an in-edge.
 * The function @ref ageqedge canonicalizes the pointers before
 * doing a comparison and so can be used to test edge equality.
 * The sense of an edge can be flipped using @ref agopp.
 *
 * @{
 */

CGRAPH_API Agedge_t *agedge(Agraph_t * g, Agnode_t * t, Agnode_t * h,
			char *name, int createflag);
CGRAPH_API Agedge_t *agidedge(Agraph_t * g, Agnode_t * t, Agnode_t * h,
              IDTYPE id, int createflag);
CGRAPH_API Agedge_t *agsubedge(Agraph_t * g, Agedge_t * e, int createflag);
CGRAPH_API Agedge_t *agfstin(Agraph_t * g, Agnode_t * n);
CGRAPH_API Agedge_t *agnxtin(Agraph_t * g, Agedge_t * e);
CGRAPH_API Agedge_t *agfstout(Agraph_t * g, Agnode_t * n);
CGRAPH_API Agedge_t *agnxtout(Agraph_t * g, Agedge_t * e);
CGRAPH_API Agedge_t *agfstedge(Agraph_t * g, Agnode_t * n);
CGRAPH_API Agedge_t *agnxtedge(Agraph_t * g, Agedge_t * e, Agnode_t * n);
/// @}

/// @defgroup cgraph_generic generic
/// @{
CGRAPH_API Agraph_t *agraphof(void* obj);
CGRAPH_API Agraph_t *agroot(void* obj);
CGRAPH_API int agcontains(Agraph_t *, void *);
CGRAPH_API char *agnameof(void *);
CGRAPH_API int agrelabel_node(Agnode_t * n, char *newname);
CGRAPH_API int agdelete(Agraph_t * g, void *obj);
CGRAPH_API int agdelsubg(Agraph_t * g, Agraph_t * sub);	/* could be agclose */
CGRAPH_API int agdelnode(Agraph_t * g, Agnode_t * arg_n);
CGRAPH_API int agdeledge(Agraph_t * g, Agedge_t * arg_e);
CGRAPH_API int agobjkind(void *);
/// @}

/** @defgroup cgraph_attr attributes
 *  @brief strings, symbols, and @ref cgraph_rec
 *  @ingroup cgraph_api
 *
 * Programmer-defined values may be dynamically
 * attached to graphs, subgraphs, nodes, and edges.
 * Such values are either character string data (for I/O)
 * or uninterpreted binary @ref cgraph_rec (for implementing algorithms efficiently).
 *
 * @{
 */

CGRAPH_API char *agstrdup(Agraph_t *, const char *);
CGRAPH_API char *agstrdup_html(Agraph_t *, const char *);
CGRAPH_API int aghtmlstr(const char *);
CGRAPH_API char *agstrbind(Agraph_t * g, const char *);
CGRAPH_API int agstrfree(Agraph_t *, const char *);
CGRAPH_API char *agcanon(char *, int);
CGRAPH_API char *agstrcanon(char *, char *);
CGRAPH_API char *agcanonStr(char *str);  /* manages its own buf */

/// definitions for dynamic string attributes

struct Agattr_s {		/* dynamic string attributes */
    Agrec_t h;			/* common data header */
    Dict_t *dict;		/* shared dict to interpret attr field */
    char **str;			/* the attribute string values */
};

/// symbol in one of the above dictionaries

struct Agsym_s {
    Dtlink_t link;
    char *name;			/* attribute's name */
    char *defval;		/* its default value for initialization */
    int id;			/* its index in attr[] */
    unsigned char kind;		/* referent object type */
    unsigned char fixed;	/* immutable value */
    unsigned char print;	/* always print */
};

struct Agdatadict_s {		/* set of dictionaries per graph */
    Agrec_t h;			/* installed in list of graph recs */
    struct {
	Dict_t *n, *e, *g;
    } dict;
};

CGRAPH_API Agsym_t *agattr(Agraph_t * g, int kind, char *name,
                           const char *value);
CGRAPH_API Agsym_t *agattrsym(void *obj, char *name);
CGRAPH_API Agsym_t *agnxtattr(Agraph_t * g, int kind, Agsym_t * attr);
CGRAPH_API int      agcopyattr(void *oldobj, void *newobj);

/// @addtogroup cgraph_rec
/// @{
CGRAPH_API void *agbindrec(void *obj, const char *name, unsigned int recsize,
		       int move_to_front);
///< attach a new record of the given size to the object.

CGRAPH_API Agrec_t *aggetrec(void *obj, const char *name, int move_to_front);
///< find record in circular list and do optional move-to-front and lock

CGRAPH_API int agdelrec(void *obj, const char *name);
CGRAPH_API void aginit(Agraph_t * g, int kind, const char *rec_name,
                       int rec_size, int move_to_front);
CGRAPH_API void agclean(Agraph_t * g, int kind, char *rec_name);
/// @}

CGRAPH_API char *agget(void *obj, char *name);
CGRAPH_API char *agxget(void *obj, Agsym_t * sym);
CGRAPH_API int agset(void *obj, char *name, const char *value);
CGRAPH_API int agxset(void *obj, Agsym_t * sym, const char *value);
CGRAPH_API int agsafeset(void* obj, char* name, const char* value,
                         const char* def);
/// @}

/// @defgroup cgraph_subgraph definitions for subgraphs
/// @{
CGRAPH_API Agraph_t *agsubg(Agraph_t * g, char *name, int cflag);	/* constructor */
CGRAPH_API Agraph_t *agidsubg(Agraph_t * g, IDTYPE id, int cflag);	/* constructor */
CGRAPH_API Agraph_t *agfstsubg(Agraph_t * g);
CGRAPH_API Agraph_t *agnxtsubg(Agraph_t * subg);
CGRAPH_API Agraph_t *agparent(Agraph_t * g);
/// @}

/// @defgroup card set cardinality
/// @{
CGRAPH_API int agnnodes(Agraph_t * g);
CGRAPH_API int agnedges(Agraph_t * g);
CGRAPH_API int agnsubg(Agraph_t * g);
CGRAPH_API int agdegree(Agraph_t * g, Agnode_t * n, int in, int out);
CGRAPH_API int agcountuniqedges(Agraph_t * g, Agnode_t * n, int in, int out);
/// @}

/// @defgroup cgmem memory
/// @{
CGRAPH_API void *agalloc(Agraph_t * g, size_t size);
CGRAPH_API void *agrealloc(Agraph_t * g, void *ptr, size_t oldsize,
		       size_t size);
CGRAPH_API void agfree(Agraph_t * g, void *ptr);

/* an engineering compromise is a joy forever */
CGRAPH_API void aginternalmapclearlocalnames(Agraph_t * g);

#define agnew(g,t)		((t*)agalloc(g,sizeof(t)))
#define agnnew(g,n,t)	((t*)agalloc(g,(n)*sizeof(t)))
/// @}

/// @cond

/* support for extra API misuse warnings if available */
#ifdef __GNUC__
  #define PRINTF_LIKE(index, first) __attribute__((format(printf, index, first)))
#else
  #define PRINTF_LIKE(index, first) /* nothing */
#endif

/// @endcond

/// @defgroup cgraph_err error handling
/// @{
typedef enum { AGWARN, AGERR, AGMAX, AGPREV } agerrlevel_t;
typedef int (*agusererrf) (char*);
CGRAPH_API agerrlevel_t agseterr(agerrlevel_t);
CGRAPH_API char *aglasterr(void);
CGRAPH_API int agerr(agerrlevel_t level, const char *fmt, ...)
  PRINTF_LIKE(2, 3);
CGRAPH_API void agerrorf(const char *fmt, ...) PRINTF_LIKE(1, 2);
CGRAPH_API void agwarningf(const char *fmt, ...) PRINTF_LIKE(1, 2);
CGRAPH_API int agerrors(void);
CGRAPH_API int agreseterrors(void);
CGRAPH_API agusererrf agseterrf(agusererrf);
/// @}

#undef PRINTF_LIKE

/// @addtogroup cgraph_other
/// @{
/* data access macros */
/* this assumes that e[0] is out and e[1] is inedge, see @ref Agedgepair_s  */
#define AGIN2OUT(inedge)		((inedge)-1) ///< Agedgepair_s.in -> Agedgepair_s.out
#define AGOUT2IN(outedge)		((outedge)+1) ///< Agedgepair_s.out -> Agedgepair_s.in
#define AGOPP(e)		((AGTYPE(e)==AGINEDGE)?AGIN2OUT(e):AGOUT2IN(e))
#define AGMKOUT(e)		(AGTYPE(e) == AGOUTEDGE? (e): AGIN2OUT(e))
#define AGMKIN(e)		(AGTYPE(e) == AGINEDGE?  (e): AGOUT2IN(e))
#define AGTAIL(e)		(AGMKIN(e)->node)
#define AGHEAD(e)		(AGMKOUT(e)->node)
#define AGEQEDGE(e,f)		(AGMKOUT(e) == AGMKOUT(f))
/* These macros are also exposed as functions, so they can be linked against. */
#define agtail(e)		AGTAIL(e)
#define aghead(e)		AGHEAD(e)
#define agopp(e)		AGOPP(e) ///< opposite edge: flip Agedgepair_s.out ⇄ Agedgepair_s.in
#define ageqedge(e,f)		AGEQEDGE(e,f) ///< edges are equal

#define TAILPORT_ID		"tailport"
#define HEADPORT_ID		"headport"
/// @}

/// @addtogroup cgraph_graph
/// @{
CGRAPH_API extern Agdesc_t Agdirected;
CGRAPH_API extern Agdesc_t Agstrictdirected;
CGRAPH_API extern Agdesc_t Agundirected;
CGRAPH_API extern Agdesc_t Agstrictundirected;
/// @}

/// @defgroup cgraph_fast fast graphs
/// @{

/* this is expedient but a bit slimey because it "knows" that dict entries of both nodes
and edges are embedded in main graph objects but allocated separately in subgraphs */
#define AGSNMAIN(sn)        ((sn)==(&((sn)->node->mainsub)))
#define EDGEOF(sn,rep)		(AGSNMAIN(sn)?((Agedge_t*)((unsigned char*)(rep) - offsetof(Agedge_t,seq_link))) : ((Dthold_t*)(rep))->obj)
/// @}

#ifdef __cplusplus
}
#endif
/// @}
