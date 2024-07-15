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
 */

/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * https://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: Details at https://graphviz.org
 *************************************************************************/

#pragma once

#include "cdt.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup cgraph_api Cgraph API
 * @brief Abstract graph C library. API cgraph.h
 * @ingroup public_apis
 *
 * [man 3 cgraph](https://graphviz.org/pdf/cgraph.3.pdf)
 *
 * Main types @ref Agraph_t, @ref Agnode_t, @ref Agedge_t.
 * @{
 */

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

/// @endcond

/* forward struct type declarations */
/// @addtogroup cgraph_object
/// @{
typedef uint64_t IDTYPE; ///< unique per main graph ID
typedef struct Agtag_s Agtag_t;
typedef struct Agobj_s Agobj_t; ///< generic object header
/// @}
/// @addtogroup cgraph_graph
/// @{
typedef struct Agraph_s Agraph_t;     ///< graph, subgraph (or hyperedge)
typedef struct Agdesc_s Agdesc_t;     ///< graph descriptor
typedef struct Agdstate_s Agdstate_t; ///< client state (closures)
typedef struct Agclos_s Agclos_t;     ///< common fields for graph/subgs
/// @}
/// @addtogroup cgraph_node
/// @{
typedef struct Agnode_s Agnode_t; ///< node (atom)
typedef struct Agsubnode_s Agsubnode_t;
/// @}
/// @addtogroup cgraph_edge
/// @{
typedef struct Agedge_s Agedge_t;         ///< node pair
typedef struct Agedgepair_s Agedgepair_t; ///< the edge object
/// @}
/// @addtogroup cgraph_disc
/// @{
typedef struct Agiddisc_s Agiddisc_t; ///< object ID allocator
typedef struct Agiodisc_s Agiodisc_t; ///< IO services
typedef struct Agdisc_s Agdisc_t;     ///< union of client discipline methods
/// @}
/// @addtogroup cgraph_callback
/// @{
typedef struct Agcbdisc_s Agcbdisc_t;   ///< client event callbacks
typedef struct Agcbstack_s Agcbstack_t; ///< enclosing state for @ref Agcbdisc_t
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
 *  The search function @ref aggetrec has an option to lock this pointer on a
 * given record. The application must be written so only one such lock is
 * outstanding at a time.
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
 *  Internally, records Agrec_s are maintained in circular linked lists
 *  attached to graph objects Agobj_s.
 *  To allow referencing application-dependent data without function calls or
 * search, Libcgraph allows setting and locking the list pointer of a graph,
 * node, or edge on a particular record (see @ref Agtag_s.mtflock and @ref
 * Agobj_s.data). This pointer can be obtained with the macro @ref AGDATA(obj).
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

typedef struct Agsym_s Agsym_t;           ///< string attribute descriptors
typedef struct Agattr_s Agattr_t;         ///< string attribute container
typedef struct Agdatadict_s Agdatadict_t; ///< set of dictionaries per graph
typedef struct Agrec_s Agrec_t;
///< generic header of @ref Agattr_s, @ref Agdatadict_s and user records

/// implementation of @ref Agrec_t
struct Agrec_s {
  char *name;
  Agrec_t *next; ///< **circular** linked list of records
                 /* following this would be any programmer-defined data */
};
/// @}
/// @}

/** @defgroup cgraph_object objects
 *  @brief parent for @ref cgraph_graph, @ref cgraph_node, and @ref cgraph_edge.
 *  @ingroup cgraph_api
 *
 * Common parameter for functions **obj** is generic pointer
 * to @ref Agraph_t, @ref Agnode_t, or @ref Agedge_t
 *
 * @ref AGDATA, @ref AGID, @ref AGTYPE, and others are macros returning
 * the specified fields of the argument object.
 * @{
 */

/** @brief tag in @ref Agobj_s for graphs, nodes, and edges.

While there may be several structs
for a given node or edges, there is only one unique ID (per main graph).  */

struct Agtag_s {
  /// access with @ref AGTYPE
  unsigned objtype : 2; /* see enum below */
  unsigned mtflock : 1; ///< @brief move-to-front lock, guards @ref Agobj_s.data
  unsigned attrwf : 1;  /* attrs written (parity, write.c) */
  unsigned seq : (sizeof(unsigned) * 8 - 4); /* sequence no. */
  IDTYPE id;                                 /* client  ID */
};

/// Object tags. Can't exceed 2 bits. See Agtag_s.
enum { AGRAPH, AGNODE, AGEDGE, AGOUTEDGE = AGEDGE, AGINEDGE };

/// a generic header of @ref Agraph_s, @ref Agnode_s and @ref Agedge_s
struct Agobj_s {
  Agtag_t tag;   ///< access with @ref AGTAG
  Agrec_t *data; ///< stores programmer-defined data, access with @ref AGDATA
};

#define AGTAG(obj) (((Agobj_t *)(obj))->tag)
#define AGTYPE(obj) (AGTAG(obj).objtype)
///< @brief returns @ref AGRAPH, @ref AGNODE, or @ref AGEDGE depending
/// on the type of the object

#define AGID(obj) (AGTAG(obj).id)
///< returns the unique integer ID associated with the object

#define AGSEQ(obj) (AGTAG(obj).seq)
#define AGATTRWF(obj) (AGTAG(obj).attrwf)
#define AGDATA(obj) (((Agobj_t *)(obj))->data)
///< returns @ref Agrec_t

/// @}
/// @} cgraph_api

/** @defgroup cgraph_node nodes
 *  @ingroup cgraph_object
 *
 * A node is created by giving a unique string name or programmer
 * defined integer ID, and is represented by a unique internal object.
 * (%Node equality can checked by pointer comparison.)
 *
 * @{
 */

/** @brief This is the node struct allocated per graph (or subgraph).

It resides in the n_dict of the graph.
The node set is maintained by libdict, but transparently to libgraph callers.
Every node may be given an optional string name at its time of creation,
or it is permissible to pass NIL(char*) for the name. */

struct Agsubnode_s { /* the node-per-graph-or-subgraph record */
  Dtlink_t seq_link; /* must be first */
  Dtlink_t id_link;
  Agnode_t *node;             /* the object */
  Dtlink_t *in_id, *out_id;   /* by node/ID for random access */
  Dtlink_t *in_seq, *out_seq; /* by node/sequence for serial access */
};

struct Agnode_s {
  Agobj_t base;
  Agraph_t *root;
  Agsubnode_t mainsub; /* embedded for main graph */
};
/// @}

/// @addtogroup cgraph_edge
/// @{
struct Agedge_s {
  Agobj_t base;
  Dtlink_t id_link; /* main graph only */
  Dtlink_t seq_link;
  Agnode_t *node; /* the endpoint node */
};

struct Agedgepair_s {
  Agedge_t out, in;
};
/// @}

/// @addtogroup cgraph_graph
/// @{

/// graph descriptor
struct Agdesc_s {         /* graph descriptor */
  unsigned directed : 1;  /* if edges are asymmetric */
  unsigned strict : 1;    /* if multi-edges forbidden */
  unsigned no_loop : 1;   /* if no loops */
  unsigned maingraph : 1; /* if this is the top level graph */
  unsigned no_write : 1;  /* if a temporary subgraph */
  unsigned has_attrs : 1; /* if string attr tables should be initialized */
  unsigned has_cmpnd : 1; /* if may contain collapsed nodes */
};
/// @}

/** @defgroup cgraph_disc disciplines
 *  @ingroup cgraph_misc
 *  @brief disciplines for external resources needed by libgraph
 *
 *  (This section is not intended for casual users.)
 *
 *  Programmer-defined disciplines customize certain resources:
 *  ID namespace, memory, and I/O - needed by Libcgraph.
 *  A discipline struct (or NULL) is passed at graph creation time.
 *  @{
 */

/**
 * @brief object ID allocator discipline
 *
 * An ID allocator discipline allows a client to control assignment
 * of IDs (uninterpreted integer values) to objects, and possibly how
 * they are mapped to and from strings.
 */

/// object ID allocator
struct Agiddisc_s {
  void *(*open)(Agraph_t *g, Agdisc_t *); /* associated with a graph */
  long (*map)(void *state, int objtype, char *str, IDTYPE *id, int createflag);
  long (*alloc)(void *state, int objtype, IDTYPE id);
  void (*free)(void *state, int objtype, IDTYPE id);
  char *(*print)(void *state, int objtype, IDTYPE id);
  void (*close)(void *state);
  void (*idregister)(void *state, int objtype, void *obj);
};

/// IO services
struct Agiodisc_s {
  int (*afread)(void *chan, char *buf, int bufsize);
  int (*putstr)(void *chan, const char *str);
  int (*flush)(void *chan); /* sync */
                            /* error messages? */
};

/// @brief user's discipline
///
/// A default discipline is supplied when NULL is given for any of these fields.
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
 *  @ingroup cgraph_object
 *
 * The functions @ref agisdirected, @ref agisundirected,
 * @ref agisstrict, and @ref agissimple
 * can be used to query if a graph is directed, undirected,
 * strict (at most one edge with a given tail and head),
 * or simple (strict with no loops), respectively.
 *
 *  @{
 */

/// client state (closures)
struct Agdstate_s {
  void *id;
  /* IO must be initialized and finalized outside Cgraph,
   * and channels (FILES) are passed as void* arguments. */
};

typedef void (*agobjfn_t)(Agraph_t *g, Agobj_t *obj, void *arg);
typedef void (*agobjupdfn_t)(Agraph_t *g, Agobj_t *obj, void *arg,
                             Agsym_t *sym);

/** @defgroup cgraph_callback callbacks
 *  @brief virtual methods of initialization, modification, and finalization of
 * graph objects
 *
 * An @ref Agcbdisc_t defines callbacks to be invoked by Libcgraph when
 * initializing, modifying, or finalizing graph objects.
 * Disciplines are kept on a stack.
 * Libcgraph automatically calls the methods on the stack, top-down.
 * Callbacks are installed with @ref agpushdisc, uninstalled with @ref
 * agpopdisc.
 *
 * @{
 */

/// client event callbacks, used in Agcbstack_s
struct Agcbdisc_s {
  struct {
    agobjfn_t ins;
    agobjupdfn_t mod;
    agobjfn_t del;
  } graph, node, edge;
};

/// object event callbacks

/// enclosing state for Agcbdisc_s, used in Agclos_s
struct Agcbstack_s {
  Agcbdisc_t *f;     /* methods */
  void *state;       /* closure */
  Agcbstack_t *prev; /* kept in a stack, unlike other disciplines */
};

CGRAPH_API void agpushdisc(Agraph_t *g, Agcbdisc_t *disc, void *state);
CGRAPH_API int agpopdisc(Agraph_t *g, Agcbdisc_t *disc);

/// @}

/// shared resources for Agraph_s
struct Agclos_s {
  Agdisc_t disc;    /* resource discipline functions */
  Agdstate_t state; /* resource closures */
  Dict_t *strdict;  /* shared string dict */
  uint64_t seq[3];  /* local object sequence number counter */
  Agcbstack_t *cb;  /* user and system callback function stacks */
  Dict_t *lookup_by_name[3];
  Dict_t *lookup_by_id[3];
};

/// graph or subgraph
struct Agraph_s {
  Agobj_t base;
  Agdesc_t desc;
  Dtlink_t seq_link;
  Dtlink_t id_link;
  Dict_t *n_seq;           /* the node set in sequence */
  Dict_t *n_id;            /* the node set indexed by ID */
  Dict_t *e_seq, *e_id;    /* holders for edge sets */
  Dict_t *g_seq, *g_id;    /* subgraphs - descendants */
  Agraph_t *parent, *root; /* subgraphs - ancestors */
  Agclos_t *clos;          /* shared resources */
};

/* graphs */
CGRAPH_API Agraph_t *agopen(char *name, Agdesc_t desc, Agdisc_t *disc);
/**<
 * @brief creates a new graph with the given name and kind
 *
 * @param desc - graph kind, can be @ref Agdirected, @ref Agundirected,
 * @ref Agstrictdirected or @ref Agstrictundirected.
 * A strict graph cannot have multi-edges or self-arcs.
 *
 * @param disc - discipline structure which can be used
 * to tailor I/O, memory allocation, and ID allocation. Typically, a NULL
 * value will be used to indicate the default discipline @ref AgDefaultDisc.
 */

CGRAPH_API int agclose(Agraph_t *g);
///< deletes a graph, freeing its associated storage
CGRAPH_API Agraph_t *agread(void *chan, Agdisc_t *disc);
///< constructs a new graph

CGRAPH_API Agraph_t *agmemread(const char *cp);
///< reads a graph from the input string

CGRAPH_API Agraph_t *agmemconcat(Agraph_t *g, const char *cp);
CGRAPH_API void agreadline(int);
///< sets input line number for subsequent error reporting

CGRAPH_API void agsetfile(const char *);
///< sets the current file name for subsequent error reporting

CGRAPH_API Agraph_t *agconcat(Agraph_t *g, void *chan, Agdisc_t *disc);
/**< @brief merges the file contents with a pre-existing graph
 *
 * Though I/O methods may be overridden, the default is that
 * the channel argument is a stdio FILE pointer.
 * In that case, if any of the streams are wide-oriented,
 * the behavior is undefined.
 */

CGRAPH_API int agwrite(Agraph_t *g, void *chan);
CGRAPH_API int agisdirected(Agraph_t *g);
CGRAPH_API int agisundirected(Agraph_t *g);
CGRAPH_API int agisstrict(Agraph_t *g);
CGRAPH_API int agissimple(Agraph_t *g);
/// @}

/// @addtogroup cgraph_node
/// @{
CGRAPH_API Agnode_t *agnode(Agraph_t *g, char *name, int createflag);
CGRAPH_API Agnode_t *agidnode(Agraph_t *g, IDTYPE id, int createflag);
CGRAPH_API Agnode_t *agsubnode(Agraph_t *g, Agnode_t *n, int createflag);
CGRAPH_API Agnode_t *agfstnode(Agraph_t *g);
CGRAPH_API Agnode_t *agnxtnode(Agraph_t *g, Agnode_t *n);
CGRAPH_API Agnode_t *aglstnode(Agraph_t *g);
CGRAPH_API Agnode_t *agprvnode(Agraph_t *g, Agnode_t *n);

CGRAPH_API Agsubnode_t *agsubrep(Agraph_t *g, Agnode_t *n);
CGRAPH_API int agnodebefore(Agnode_t *u, Agnode_t *v); /* we have no shame */
CGRAPH_API int agdelnode(Agraph_t *g, Agnode_t *arg_n);
///< removes a node from a graph or subgraph.
CGRAPH_API int agrelabel_node(Agnode_t *n, char *newname);
/// @}

/** @defgroup cgraph_edge edges
 *  @ingroup cgraph_object
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

CGRAPH_API Agedge_t *agedge(Agraph_t *g, Agnode_t *t, Agnode_t *h, char *name,
                            int createflag);
CGRAPH_API Agedge_t *agidedge(Agraph_t *g, Agnode_t *t, Agnode_t *h, IDTYPE id,
                              int createflag);
CGRAPH_API Agedge_t *agsubedge(Agraph_t *g, Agedge_t *e, int createflag);
CGRAPH_API Agedge_t *agfstin(Agraph_t *g, Agnode_t *n);
CGRAPH_API Agedge_t *agnxtin(Agraph_t *g, Agedge_t *e);
CGRAPH_API Agedge_t *agfstout(Agraph_t *g, Agnode_t *n);
CGRAPH_API Agedge_t *agnxtout(Agraph_t *g, Agedge_t *e);
CGRAPH_API Agedge_t *agfstedge(Agraph_t *g, Agnode_t *n);
CGRAPH_API Agedge_t *agnxtedge(Agraph_t *g, Agedge_t *e, Agnode_t *n);
CGRAPH_API int agdeledge(Agraph_t *g, Agedge_t *arg_e);
/// @}

/// @addtogroup cgraph_object
/// @{

// Generic object (graphs, nodes and edges) functions

CGRAPH_API Agraph_t *agraphof(void *obj);
CGRAPH_API Agraph_t *agroot(void *obj);
///< takes any graph object (graph, subgraph, node, edge) and returns the root
///< graph in which it lives

CGRAPH_API int agcontains(Agraph_t *, void *obj);
///< returns non-zero if **obj** is a member of (sub)graph

CGRAPH_API char *agnameof(void *);
///< returns a string descriptor for the object.

CGRAPH_API int agdelete(Agraph_t *g, void *obj);
/**< @brief deletes object.
 * Equivalent to @ref agclose, @ref agdelnode, and @ref agdeledge
 * for **obj** being a graph, node or edge, respectively.
 * @returns -1 if **obj** does not belong to graph **g**.
 */

CGRAPH_API int agobjkind(void *obj);
///< returns @ref AGRAPH, @ref AGNODE, or @ref AGEDGE depending on the type of
///< the object. Synonym for @ref AGTYPE.
/// @}

/** @defgroup cgraph_string string utilities
 *  @brief reference-counted strings
 *  @ingroup cgraph_misc
 *
 *  Storage management of strings as reference-counted strings.
 *  The caller does not need to dynamically allocate storage.
 *
 * All uses of cgraph strings need to be freed using @ref agstrfree
 * in order to correctly maintain the reference count.
 *
 * @ref agcanonStr returns a pointer to a version of the input string
 * canonicalized for output for later re-parsing.
 * This includes quoting special characters and keywords.
 * It uses its own internal buffer, so the value will be lost on
 * the next call to @ref agcanonStr.
 * @ref agcanon is identical with @ref agcanonStr
 * except it can be used with any character string.
 * The second argument indicates whether or not the string
 * should be canonicalized as an HTML-like string.
 *
 * @{
 */
CGRAPH_API char *agstrdup(Agraph_t *, const char *);
///< @brief returns a pointer to a reference-counted copy of the argument
///< string,
/// creating one if necessary

CGRAPH_API char *agstrdup_html(Agraph_t *, const char *);
///< create an HTML-like string
///
CGRAPH_API int aghtmlstr(const char *);
///< query if a string is an ordinary string or an HTML-like string
///
CGRAPH_API int aghtmlstr(const char *);
CGRAPH_API char *agstrbind(Agraph_t *g, const char *);
///< returns a pointer to a reference-counted string if it exists, or NULL if
///< not

CGRAPH_API int agstrfree(Agraph_t *, const char *);
CGRAPH_API char *agcanon(char *str, int html);
CGRAPH_API char *agstrcanon(char *, char *);
CGRAPH_API char *agcanonStr(char *str); /* manages its own buf */
/// @}

/** @defgroup cgraph_attr attributes
 *  @brief symbols, and @ref cgraph_rec
 *  @ingroup cgraph_api
 *
 * Programmer-defined values may be dynamically
 * attached to graphs, subgraphs, nodes, and edges.
 * Such values are either character string data (see @ref agattr) (for I/O)
 * or uninterpreted binary @ref cgraph_rec (for implementing algorithms
 * efficiently).
 *
 * *String attributes* are handled automatically in reading and writing graph
 * files. A string attribute is identified by name and by an internal symbol
 * table entry (@ref Agsym_t) created by Libcgraph. Attributes of nodes, edges,
 * and graphs (with their subgraphs) have separate namespaces. The contents of
 * an @ref Agsym_t have a char* *name* for the attribute's name, a char*
 * *defval* field for the attribute's default value, and an int *id* field
 * containing the index of the attribute's specific value for an object in the
 * object's array of attribute values.
 *
 * @{
 */

// definitions for dynamic string attributes

/// string attribute container
struct Agattr_s { /* dynamic string attributes */
  Agrec_t h;      /* common data header */
  Dict_t *dict;   ///< shared dict of Agsym_s to interpret Agattr_s.str
  char **str;     ///< the attribute string values indexed by Agsym_s.id
};

/// @brief string attribute descriptor
/// symbol in Agattr_s.dict
struct Agsym_s {
  Dtlink_t link;
  char *name;          /* attribute's name */
  char *defval;        /* its default value for initialization */
  int id;              ///< index in Agattr_s.str
  unsigned char kind;  /* referent object type */
  unsigned char fixed; /* immutable value */
  unsigned char print; /* always print */
};

struct Agdatadict_s { ///< set of dictionaries per graph
  Agrec_t h;          /* installed in list of graph recs */
  struct {
    Dict_t *n, *e, *g;
  } dict;
};

CGRAPH_API Agsym_t *agattr(Agraph_t *g, int kind, char *name,
                           const char *value);
/**< @brief creates or looks up attributes of a graph
 * @param g graph. When is NULL, the default is set for all graphs created
 * subsequently.
 * @param kind may be @ref AGRAPH, @ref AGNODE, or @ref AGEDGE.
 * @param value default value. When is (char*)0, the request is to search
 * for an existing attribute of the given kind and name.
 *
 * If the attribute already exists, its default
 * for creating new objects is set to the given **value**;
 * if it does not exist, a new attribute is created with the
 * given default **value**, and the default is applied to all pre-existing
 * objects of the given **kind**
 */

CGRAPH_API Agsym_t *agattrsym(void *obj, char *name);
///< looks up a string attribute for a graph object given as an argument

CGRAPH_API Agsym_t *agnxtattr(Agraph_t *g, int kind, Agsym_t *attr);
///< @brief permits traversing the list of attributes of a given type
/// @param attr	if `NULL` the function returns the first attribute
/// @returns the next one in succession or `NULL` at the end of the list.

CGRAPH_API int agcopyattr(void *oldobj, void *newobj);
/**< @brief copies all of the attributes from one object to another
 * @return fails and returns non-zero if argument objects are different kinds,
 * or if all of the attributes of the source object have not been declared
 * for the target object
 */

/// @addtogroup cgraph_rec
/// @{
CGRAPH_API void *agbindrec(void *obj, const char *name, unsigned int recsize,
                           int move_to_front);
///< @brief attaches a new record of the given size to the object
/// @param recsize if 0, the call to @ref agbindrec is simply a lookup
/// @returns pointer to `Agrec_t` and user data

CGRAPH_API Agrec_t *aggetrec(void *obj, const char *name, int move_to_front);
///< find record in circular list and do optional move-to-front and lock

CGRAPH_API int agdelrec(void *obj, const char *name);
///<  deletes a named record from one object

CGRAPH_API void aginit(Agraph_t *g, int kind, const char *rec_name,
                       int rec_size, int move_to_front);
/**< @brief attach new records to objects of specified kind
 *
 * @param kind may be @ref AGRAPH, @ref AGNODE, or @ref AGEDGE
 * @param rec_size if is negative (of the actual rec_size) for graphs,
 * @ref aginit is applied recursively to the graph and its subgraphs
 */

CGRAPH_API void agclean(Agraph_t *g, int kind, char *rec_name);
///< @brief calls @ref agdelrec for all objects
/// of the same class in an entire graph

/// @}

CGRAPH_API char *agget(void *obj, char *name);
CGRAPH_API char *agxget(void *obj, Agsym_t *sym);
CGRAPH_API int agset(void *obj, char *name, const char *value);
CGRAPH_API int agxset(void *obj, Agsym_t *sym, const char *value);
CGRAPH_API int agsafeset(void *obj, char *name, const char *value,
                         const char *def);
///< @brief ensures the given attribute is declared
///  before setting it locally on an object

/// @}

/** @defgroup cgraph_subgraph subgraphs
 *  @ingroup cgraph_graph
 *
 * A "main" or "root" graph defines a namespace for a collection of
 * graph objects (subgraphs, nodes, edges) and their attributes.
 * Objects may be named by unique strings or by integer IDs.
 *
 * @ref agsubg finds or creates a subgraph by name.
 *
 * @ref agidsubg allows a programmer to specify the subgraph by a unique integer
 * ID.
 *
 * A new subgraph is initially empty and is of the same kind as its parent.
 * Nested subgraph trees may be created.
 * A subgraph's name is only interpreted relative to its parent.
 *
 * A program can scan subgraphs under a given graph
 * using @ref agfstsubg and @ref agnxtsubg.
 *
 * A subgraph is deleted with @ref agdelsubg (or @ref agclose).
 *
 * The @ref agparent function returns the immediate parent graph of a subgraph,
 * or itself if the graph is already a root graph.
 *
 * @{
 */

CGRAPH_API Agraph_t *agsubg(Agraph_t *g, char *name,
                            int cflag); /* constructor */
CGRAPH_API Agraph_t *agidsubg(Agraph_t *g, IDTYPE id,
                              int cflag); /* constructor */
CGRAPH_API Agraph_t *agfstsubg(Agraph_t *g);
CGRAPH_API Agraph_t *agnxtsubg(Agraph_t *subg);
CGRAPH_API Agraph_t *agparent(Agraph_t *g);
CGRAPH_API int agdelsubg(Agraph_t *g, Agraph_t *sub); /* could be agclose */
/// @}

/** @defgroup cgraph_misc miscellaneous
 *  @ingroup cgraph_api
 *  @{
 */

/** @defgroup card set cardinality
 *
 * By default, nodes are stored in ordered sets for
 * efficient random access to insert, find, and delete nodes.
 *
 * @ref agnnodes, @ref agnedges, and @ref agnsubg return the
 * sizes of node, edge and subgraph sets of a graph.
 *
 * The function @ref agdegree returns the size of a node’s edge set,
 * and takes flags to select in-edges, out-edges, or both.
 *
 * The function @ref agcountuniqedges returns
 * the size of a node’s edge set, and takes flags
 * to select in-edges, out-edges, or both.
 * Unlike @ref agdegree, each loop is only counted once.
 *
 * @{
 */
CGRAPH_API int agnnodes(Agraph_t *g);
CGRAPH_API int agnedges(Agraph_t *g);
CGRAPH_API int agnsubg(Agraph_t *g);
CGRAPH_API int agdegree(Agraph_t *g, Agnode_t *n, int in, int out);
CGRAPH_API int agcountuniqedges(Agraph_t *g, Agnode_t *n, int in, int out);
/// @}

/// @defgroup cgmem memory
/// @{
CGRAPH_API void *agalloc(Agraph_t *g, size_t size);
CGRAPH_API void *agrealloc(Agraph_t *g, void *ptr, size_t oldsize, size_t size);
CGRAPH_API void agfree(Agraph_t *g, void *ptr);

/* an engineering compromise is a joy forever */
CGRAPH_API void aginternalmapclearlocalnames(Agraph_t *g);

#define agnew(g, t) ((t *)agalloc(g, sizeof(t)))
#define agnnew(g, n, t) ((t *)agalloc(g, (n) * sizeof(t)))
/// @}

/// @cond

/* support for extra API misuse warnings if available */
#ifdef __GNUC__
#define PRINTF_LIKE(index, first) __attribute__((format(printf, index, first)))
#else
#define PRINTF_LIKE(index, first) /* nothing */
#endif

/// @endcond

/** @defgroup cgraph_err error handling
 *
 * The library provides a variety of mechanisms to control
 * the reporting of errors and warnings.
 * A message is only written if its type has higher priority than
 * a programmer-controlled minimum, which is @ref AGWARN by default.
 * The programmer can set this value using @ref agseterr,
 * which returns the previous value.
 * Calling `agseterr(AGMAX)` turns off the writing of messages.
 *
 * The function @ref agerr is the main entry point for reporting an anomaly.
 * The first argument indicates the type of message.
 * Usually, the first argument is @ref AGWARN or @ref AGERR
 * to indicate warnings and errors, respectively.
 * Sometimes additional context information is only available in functions
 * calling the function where the error is actually caught.
 * In this case, the calling function can indicate that it is continuing
 * the current error by using @ref AGPREV as the first argument.
 * The remaining arguments to @ref agerr are the same as
 * the arguments to `printf`.
 *
 * The functions @ref agwarningf and @ref agerrorf are shorthand for
 * `agerr(AGWARN,...)` and `agerr(AGERR,...)`, respectively.
 *
 * Some applications desire to directly control the writing of messages.
 * Such an application can use the function @ref agseterrf to register
 * the function that the library should call to actually write the message.
 * The previous error function is returned.
 * By default, the message is written to `stderr`.
 *
 * Errors not written are stored in a log file.
 * The last recorded error can be retrieved by calling @ref aglasterr.
 * Unless the printing of error messages has been completely disabled
 * by a call to `agseterr(AGMAX)`, standard error must not be wide-oriented,
 * even if a user-provided error printing function is provided.
 *
 * The function @ref agerrors returns non-zero if errors have been reported.
 *
 * @{
 */
typedef enum { AGWARN, AGERR, AGMAX, AGPREV } agerrlevel_t;
typedef int (*agusererrf)(char *);
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
/// @}

/// @addtogroup cgraph_edge
/// @{
/* data access macros */
/* this assumes that e[0] is out and e[1] is inedge, see @ref Agedgepair_s  */
#define AGIN2OUT(inedge) ((inedge)-1) ///< Agedgepair_s.in -> Agedgepair_s.out
#define AGOUT2IN(outedge)                                                      \
  ((outedge) + 1) ///< Agedgepair_s.out -> Agedgepair_s.in
#define AGOPP(e) ((AGTYPE(e) == AGINEDGE) ? AGIN2OUT(e) : AGOUT2IN(e))
#define AGMKOUT(e) (AGTYPE(e) == AGOUTEDGE ? (e) : AGIN2OUT(e))
#define AGMKIN(e) (AGTYPE(e) == AGINEDGE ? (e) : AGOUT2IN(e))
#define AGTAIL(e) (AGMKIN(e)->node)
#define AGHEAD(e) (AGMKOUT(e)->node)
#define AGEQEDGE(e, f) (AGMKOUT(e) == AGMKOUT(f))
/* These macros are also exposed as functions, so they can be linked against. */
#define agtail(e) AGTAIL(e)
#define aghead(e) AGHEAD(e)
#define agopp(e)                                                               \
  AGOPP(e) ///< opposite edge: flip Agedgepair_s.out ⇄ Agedgepair_s.in
#define ageqedge(e, f) AGEQEDGE(e, f) ///< edges are equal

#define TAILPORT_ID "tailport"
#define HEADPORT_ID "headport"
/// @}

/// @addtogroup cgraph_graph
/// @{
CGRAPH_API extern Agdesc_t Agdirected; ///< directed
CGRAPH_API extern Agdesc_t Agstrictdirected;
///< strict directed. A strict graph cannot have multi-edges or self-arcs.
CGRAPH_API extern Agdesc_t Agundirected;       ///< undirected
CGRAPH_API extern Agdesc_t Agstrictundirected; ///< strict undirected
/// @}

/** @defgroup cgraph_fast fast graphs
 *  @ingroup cgraph_misc
 *  @{
 */

/* this is expedient but a bit slimey because it "knows" that dict entries of
both nodes and edges are embedded in main graph objects but allocated separately
in subgraphs */
#define AGSNMAIN(sn) ((sn) == (&((sn)->node->mainsub)))
#define EDGEOF(sn, rep)                                                        \
  (AGSNMAIN(sn)                                                                \
       ? ((Agedge_t *)((unsigned char *)(rep)-offsetof(Agedge_t, seq_link)))   \
       : ((Dthold_t *)(rep))->obj)
/// @}

/// @addtogroup cgraph_app
/// @{

/// options for passing to `graphviz_acyclic`
typedef struct {
  FILE *outFile;
  bool doWrite;
  bool Verbose;
} graphviz_acyclic_options_t;

/// programmatic access to `acyclic`
///
/// See `man acyclic` for an explanation of the `acyclic` tool.
///
/// \param g Graph to operate on
/// \param opts Options to control acyclic algorithm
/// \param num_rev [inout] Running total of reversed edges
/// \return True if a cycle was found, indicating failure
CGRAPH_API bool graphviz_acyclic(Agraph_t *g,
                                 const graphviz_acyclic_options_t *opts,
                                 size_t *num_rev);

/// options for passing to `graphviz_tred`
typedef struct {
  bool Verbose;
  bool PrintRemovedEdges;
  FILE *out; ///< stream to write result(s) to
  FILE *err; ///< stream to print warnings to
} graphviz_tred_options_t;

/// @brief programmatic access to `tred` -
/// [transitive reduction](https://en.wikipedia.org/wiki/Transitive_reduction)
///
/// See `man tred` for an explanation of the `tred` tool.
///
/// \param g Graph to operate on
/// \param opts Options to control tred algorithm
CGRAPH_API void graphviz_tred(Agraph_t *g, const graphviz_tred_options_t *opts);

/// options for passing to `graphviz_unflatten`
typedef struct {
  bool Do_fans;
  int MaxMinlen;
  int ChainLimit;
} graphviz_unflatten_options_t;

/// programmatic access to `unflatten`
///
/// See `man unflatten` for an explanation of the `unflatten` tool.
///
/// \param g Graph to operate on
/// \param opts Options to control unflattening
CGRAPH_API void graphviz_unflatten(Agraph_t *g,
                                   const graphviz_unflatten_options_t *opts);

/** add to a graph any edges with both endpoints within that graph
 *
 * If `edgeset` is given as `NULL`, edges from the root graph of `g` will be
 * considered. In this case if `g` itself is the root graph, this call is a
 * no-op.
 *
 * If `g` is a connected component, the edges added will be all edges attached
 * to any node in `g`.
 *
 * \param g Graph to add edges to
 * \param edgeset Graph whose edges to consider
 * \return Number of edges added
 */
CGRAPH_API size_t graphviz_node_induce(Agraph_t *g, Agraph_t *edgeset);
/// @}

#ifdef __cplusplus
}
#endif
