/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: Details at https://graphviz.org
 *************************************************************************/

#include <cstdlib>
#include <cstring>
#include <gvc/gvc.h>
#include <string>
#include "gv_channel.h"

#define agfindattr(x, s) agattrsym(x, s)
#define agraphattr(g, n, s) agattr(g, AGRAPH, n, s)
#define agnodeattr(g, n, s) agattr(g, AGNODE, n, s)
#define agedgeattr(g, n, s) agattr(g, AGEDGE, n, s)

static char emptystring[] = {'\0'};

static GVC_t *gvc;

static void gv_init(void) {
  // list of builtins, enable demand loading
  gvc = gvContextPlugins(lt_preloaded_symbols, DEMAND_LOADING);
}

Agraph_t *graph(char *name) {
  if (!gvc)
    gv_init();
  return agopen(name, Agundirected, 0);
}

Agraph_t *digraph(char *name) {
  if (!gvc)
    gv_init();
  return agopen(name, Agdirected, 0);
}

Agraph_t *strictgraph(char *name) {
  if (!gvc)
    gv_init();
  return agopen(name, Agstrictundirected, 0);
}

Agraph_t *strictdigraph(char *name) {
  if (!gvc)
    gv_init();
  return agopen(name, Agstrictdirected, 0);
}

Agraph_t *readstring(char *string) {
  if (!gvc)
    gv_init();
  return agmemread(string);
}

Agraph_t *read(FILE *f) {
  if (!gvc)
    gv_init();
  return agread(f, nullptr);
}

Agraph_t *read(const char *filename) {
  FILE *f = fopen(filename, "r");
  if (!f)
    return nullptr;
  if (!gvc)
    gv_init();
  Agraph_t *g = agread(f, nullptr);
  fclose(f);
  return g;
}

//-------------------------------------------------
Agraph_t *graph(Agraph_t *g, char *name) {
  if (!gvc)
    gv_init();
  return agsubg(g, name, 1);
}

Agnode_t *node(Agraph_t *g, char *name) {
  if (!gvc)
    return nullptr;
  return agnode(g, name, 1);
}

Agedge_t *edge(Agraph_t *g, Agnode_t *t, Agnode_t *h) {
  if (!gvc || !t || !h || !g)
    return nullptr;
  // edges from/to the protonode are not permitted
  if (AGTYPE(t) == AGRAPH || AGTYPE(h) == AGRAPH)
    return nullptr;
  return agedge(g, t, h, nullptr, 1);
}

Agedge_t *edge(Agnode_t *t, Agnode_t *h) { return edge(agraphof(t), t, h); }

// induce tail if necessary
Agedge_t *edge(char *tname, Agnode_t *h) {
  return edge(node(agraphof(h), tname), h);
}

// induce head if necessary
Agedge_t *edge(Agnode_t *t, char *hname) {
  return edge(t, node(agraphof(t), hname));
}

// induce tail/head if necessary
Agedge_t *edge(Agraph_t *g, char *tname, char *hname) {
  return edge(g, node(g, tname), node(g, hname));
}

//-------------------------------------------------
static char *myagxget(void *obj, Agsym_t *a) {
  if (!obj || !a)
    return emptystring;
  char *val = agxget(obj, a);
  if (!val)
    return emptystring;
  if (strcmp(a->name, "label") == 0 && aghtmlstr(val)) {
    size_t len = strlen(val);
    auto hs = reinterpret_cast<char *>(malloc(len + 3));
    hs[0] = '<';
    strcpy(hs + 1, val);
    hs[len + 1] = '>';
    hs[len + 2] = '\0';
    return hs;
  }
  return val;
}
char *getv(Agraph_t *g, Agsym_t *a) { return myagxget(g, a); }
char *getv(Agraph_t *g, char *attr) {
  if (!g || !attr)
    return nullptr;
  Agsym_t *a = agfindattr(agroot(g), attr);
  return myagxget(g, a);
}
static void myagxset(void *obj, Agsym_t *a, char *val) {
  if (strcmp(a->name, "label") == 0 && val[0] == '<') {
    size_t len = strlen(val);
    if (val[len - 1] == '>') {
      std::string hs(val + 1, len - 2);
      val = agstrdup_html(agraphof(obj), hs.c_str());
    }
  }
  agxset(obj, a, val);
}
char *setv(Agraph_t *g, Agsym_t *a, char *val) {
  if (!g || !a || !val)
    return nullptr;
  myagxset(g, a, val);
  return val;
}
char *setv(Agraph_t *g, char *attr, char *val) {
  if (!g || !attr || !val)
    return nullptr;
  Agsym_t *a = agfindattr(agroot(g), attr);
  if (!a)
    a = agraphattr(g->root, attr, emptystring);
  myagxset(g, a, val);
  return val;
}
//-------------------------------------------------
char *getv(Agnode_t *n, Agsym_t *a) {
  if (!n || !a)
    return nullptr;
  if (AGTYPE(n) == AGRAPH) // protonode
    return nullptr;        // FIXME ??
  return myagxget(n, a);
}
char *getv(Agnode_t *n, char *attr) {
  if (!n || !attr)
    return nullptr;
  if (AGTYPE(n) == AGRAPH) // protonode
    return nullptr;        // FIXME ??
  Agraph_t *g = agroot(agraphof(n));
  Agsym_t *a = agattr(g, AGNODE, attr, nullptr);
  return myagxget(n, a);
}
char *setv(Agnode_t *n, Agsym_t *a, char *val) {
  if (!n || !a || !val)
    return nullptr;
  if (AGTYPE(n) == AGRAPH) // protonode
    return nullptr;        // FIXME ??
  myagxset(n, a, val);
  return val;
}
char *setv(Agnode_t *n, char *attr, char *val) {

  if (!n || !attr || !val)
    return nullptr;
  if (AGTYPE(n) == AGRAPH) { // protonode
    auto g = reinterpret_cast<Agraph_t *>(n);
    (void)agattr(g, AGNODE, attr,
                 val); // create default attribute in pseudo protonode
                       // FIXME? - deal with html in "label" attributes
    return val;
  }
  Agraph_t *g = agroot(agraphof(n));
  Agsym_t *a = agattr(g, AGNODE, attr, nullptr);
  if (!a)
    a = agnodeattr(g, attr, emptystring);
  myagxset(n, a, val);
  return val;
}
//-------------------------------------------------
char *getv(Agedge_t *e, Agsym_t *a) {
  if (!e || !a)
    return nullptr;
  if (AGTYPE(e) == AGRAPH) // protoedge
    return nullptr;        // FIXME ??
  return myagxget(e, a);
}
char *getv(Agedge_t *e, char *attr) {
  if (!e || !attr)
    return nullptr;
  if (AGTYPE(e) == AGRAPH) // protoedge
    return nullptr;        // FIXME ??
  Agraph_t *g = agraphof(agtail(e));
  Agsym_t *a = agattr(g, AGEDGE, attr, nullptr);
  return myagxget(e, a);
}
char *setv(Agedge_t *e, Agsym_t *a, char *val) {
  if (!e || !a || !val)
    return nullptr;
  if (AGTYPE(e) == AGRAPH) // protoedge
    return nullptr;        // FIXME ??
  myagxset(e, a, val);
  return val;
}
char *setv(Agedge_t *e, char *attr, char *val) {
  if (!e || !attr || !val)
    return nullptr;
  if (AGTYPE(e) == AGRAPH) { // protoedge
    auto g = reinterpret_cast<Agraph_t *>(e);
    (void)agattr(g, AGEDGE, attr,
                 val); // create default attribute in pseudo protoedge
                       // FIXME? - deal with html in "label" attributes
    return val;
  }
  Agraph_t *g = agroot(agraphof(agtail(e)));
  Agsym_t *a = agattr(g, AGEDGE, attr, nullptr);
  if (!a)
    a = agattr(g, AGEDGE, attr, emptystring);
  myagxset(e, a, val);
  return val;
}
//-------------------------------------------------
Agraph_t *findsubg(Agraph_t *g, char *name) {
  if (!g || !name)
    return nullptr;
  return agsubg(g, name, 0);
}

Agnode_t *findnode(Agraph_t *g, char *name) {
  if (!g || !name)
    return nullptr;
  return agnode(g, name, 0);
}

Agedge_t *findedge(Agnode_t *t, Agnode_t *h) {
  if (!t || !h)
    return nullptr;
  if (AGTYPE(t) == AGRAPH || AGTYPE(h) == AGRAPH)
    return nullptr;
  return agfindedge(agraphof(t), t, h);
}

Agsym_t *findattr(Agraph_t *g, char *name) {
  if (!g || !name)
    return nullptr;
  return agfindattr(g, name);
}

Agsym_t *findattr(Agnode_t *n, char *name) {
  if (!n || !name)
    return nullptr;
  return agfindattr(n, name);
}

Agsym_t *findattr(Agedge_t *e, char *name) {
  if (!e || !name)
    return nullptr;
  return agfindattr(e, name);
}

//-------------------------------------------------

Agnode_t *headof(Agedge_t *e) {
  if (!e)
    return nullptr;
  if (AGTYPE(e) == AGRAPH)
    return nullptr;
  return aghead(e);
}

Agnode_t *tailof(Agedge_t *e) {
  if (!e)
    return nullptr;
  if (AGTYPE(e) == AGRAPH)
    return nullptr;
  return agtail(e);
}

Agraph_t *graphof(Agraph_t *g) {
  if (!g || g == g->root)
    return nullptr;
  return agroot(g);
}

Agraph_t *graphof(Agedge_t *e) {
  if (!e)
    return nullptr;
  if (AGTYPE(e) == AGRAPH)
    return reinterpret_cast<Agraph_t *>(
        e); // graph of protoedge is itself recast
  return agraphof(agtail(e));
}

Agraph_t *graphof(Agnode_t *n) {
  if (!n)
    return nullptr;
  if (AGTYPE(n) == AGRAPH)
    return reinterpret_cast<Agraph_t *>(
        n); // graph of protonode is itself recast
  return agraphof(n);
}

Agraph_t *rootof(Agraph_t *g) {
  if (!g)
    return nullptr;
  return agroot(g);
}

//-------------------------------------------------
Agnode_t *protonode(Agraph_t *g) {
  return reinterpret_cast<Agnode_t *>(g); // gross abuse of the type system!
}

Agedge_t *protoedge(Agraph_t *g) {
  return reinterpret_cast<Agedge_t *>(g); // gross abuse of the type system!
}

//-------------------------------------------------
char *nameof(Agraph_t *g) {
  if (!g)
    return nullptr;
  return agnameof(g);
}
char *nameof(Agnode_t *n) {
  if (!n)
    return nullptr;
  if (AGTYPE(n) == AGRAPH)
    return nullptr;
  return agnameof(n);
}
char *nameof(Agsym_t *a) {
  if (!a)
    return nullptr;
  return a->name;
}

//-------------------------------------------------
bool ok(Agraph_t *g) { return g != nullptr; }
bool ok(Agnode_t *n) { return n != nullptr; }
bool ok(Agedge_t *e) { return e != nullptr; }
bool ok(Agsym_t *a) { return a != nullptr; }
//-------------------------------------------------
Agraph_t *firstsubg(Agraph_t *g) {
  if (!g)
    return nullptr;
  return agfstsubg(g);
}

Agraph_t *nextsubg(Agraph_t *g, Agraph_t *sg) {

  if (!g || !sg)
    return nullptr;
  return agnxtsubg(sg);
}

Agraph_t *firstsupg(Agraph_t *g) { return g->parent; }

Agraph_t *nextsupg(Agraph_t *, Agraph_t *) { return nullptr; }

Agedge_t *firstout(Agraph_t *g) {
  if (!g)
    return nullptr;
  for (Agnode_t *n = agfstnode(g); n; n = agnxtnode(g, n)) {
    Agedge_t *e = agfstout(g, n);
    if (e)
      return e;
  }
  return nullptr;
}

Agedge_t *nextout(Agraph_t *g, Agedge_t *e) {
  if (!g || !e)
    return nullptr;
  Agedge_t *ne = agnxtout(g, e);
  if (ne)
    return (ne);
  for (Agnode_t *n = agnxtnode(g, agtail(e)); n; n = agnxtnode(g, n)) {
    ne = agfstout(g, n);
    if (ne)
      return ne;
  }
  return nullptr;
}

Agedge_t *firstedge(Agraph_t *g) { return firstout(g); }

Agedge_t *nextedge(Agraph_t *g, Agedge_t *e) { return nextout(g, e); }

Agedge_t *firstout(Agnode_t *n) {
  if (!n)
    return nullptr;
  return agfstout(agraphof(n), n);
}

Agedge_t *nextout(Agnode_t *n, Agedge_t *e) {
  if (!n || !e)
    return nullptr;
  return agnxtout(agraphof(n), e);
}

Agnode_t *firsthead(Agnode_t *n) {
  if (!n)
    return nullptr;
  Agedge_t *e = agfstout(agraphof(n), n);
  if (!e)
    return nullptr;
  return aghead(e);
}

Agnode_t *nexthead(Agnode_t *n, Agnode_t *h) {
  if (!n || !h)
    return nullptr;
  Agraph_t *g = agraphof(n);
  Agedge_t *e = agfindedge(g, n, h);
  if (!e)
    return nullptr;
  do {
    e = agnxtout(g, AGMKOUT(e));
    if (!e)
      return nullptr;
  } while (aghead(e) == h);
  return aghead(e);
}

Agedge_t *firstedge(Agnode_t *n) {
  if (!n)
    return nullptr;
  return agfstedge(agraphof(n), n);
}

Agedge_t *nextedge(Agnode_t *n, Agedge_t *e) {
  if (!n || !e)
    return nullptr;
  return agnxtedge(agraphof(n), e, n);
}

Agedge_t *firstin(Agraph_t *g) {
  if (!g)
    return nullptr;
  Agnode_t *n = agfstnode(g);
  if (!n)
    return nullptr;
  return agfstin(g, n);
}

Agedge_t *nextin(Agraph_t *g, Agedge_t *e) {
  if (!g || !e)
    return nullptr;
  Agedge_t *ne = agnxtin(g, e);
  if (ne)
    return (ne);
  Agnode_t *n = agnxtnode(g, aghead(e));
  if (!n)
    return nullptr;
  return agfstin(g, n);
}

Agedge_t *firstin(Agnode_t *n) {
  if (!n)
    return nullptr;
  return agfstin(agraphof(n), n);
}

Agedge_t *nextin(Agnode_t *n, Agedge_t *e) {
  if (!n || !e)
    return nullptr;
  return agnxtin(agraphof(n), e);
}

Agnode_t *firsttail(Agnode_t *n) {
  if (!n)
    return nullptr;
  Agedge_t *e = agfstin(agraphof(n), n);
  if (!e)
    return nullptr;
  return agtail(e);
}

Agnode_t *nexttail(Agnode_t *n, Agnode_t *t) {
  if (!n || !t)
    return nullptr;
  Agraph_t *g = agraphof(n);
  Agedge_t *e = agfindedge(g, t, n);
  if (!e)
    return nullptr;
  do {
    e = agnxtin(g, AGMKIN(e));
    if (!e)
      return nullptr;
  } while (agtail(e) == t);
  return agtail(e);
}

Agnode_t *firstnode(Agraph_t *g) {
  if (!g)
    return nullptr;
  return agfstnode(g);
}

Agnode_t *nextnode(Agraph_t *g, Agnode_t *n) {
  if (!g || !n)
    return nullptr;
  return agnxtnode(g, n);
}

Agnode_t *firstnode(Agedge_t *e) {
  if (!e)
    return nullptr;
  return agtail(e);
}

Agnode_t *nextnode(Agedge_t *e, Agnode_t *n) {
  if (!e || n != agtail(e))
    return nullptr;
  return aghead(e);
}

Agsym_t *firstattr(Agraph_t *g) {
  if (!g)
    return nullptr;
  g = agroot(g);
  return agnxtattr(g, AGRAPH, nullptr);
}

Agsym_t *nextattr(Agraph_t *g, Agsym_t *a) {
  if (!g || !a)
    return nullptr;
  g = agroot(g);
  return agnxtattr(g, AGRAPH, a);
}

Agsym_t *firstattr(Agnode_t *n) {
  if (!n)
    return nullptr;
  Agraph_t *g = agraphof(n);
  return agnxtattr(g, AGNODE, nullptr);
}

Agsym_t *nextattr(Agnode_t *n, Agsym_t *a) {
  if (!n || !a)
    return nullptr;
  Agraph_t *g = agraphof(n);
  return agnxtattr(g, AGNODE, a);
}

Agsym_t *firstattr(Agedge_t *e) {
  if (!e)
    return nullptr;
  Agraph_t *g = agraphof(agtail(e));
  return agnxtattr(g, AGEDGE, nullptr);
}

Agsym_t *nextattr(Agedge_t *e, Agsym_t *a) {
  if (!e || !a)
    return nullptr;
  Agraph_t *g = agraphof(agtail(e));
  return agnxtattr(g, AGEDGE, a);
}

bool rm(Agraph_t *g) {
  if (!g)
    return false;
  // The rm function appears to have the semantics of agclose, so
  // we should just do that, and let cgraph take care of all the
  // details.
  agclose(g);
  return true;
}

bool rm(Agnode_t *n) {
  if (!n)
    return false;
  // removal of the protonode is not permitted
  if (strcmp(agnameof(n), "\001proto") == 0)
    return false;
  agdelete(agraphof(n), n);
  return true;
}

bool rm(Agedge_t *e) {
  if (!e)
    return false;
  // removal of the protoedge is not permitted
  if (strcmp(agnameof(aghead(e)), "\001proto") == 0 ||
      strcmp(agnameof(agtail(e)), "\001proto") == 0)
    return false;
  agdelete(agroot(agraphof(aghead(e))), e);
  return true;
}

bool layout(Agraph_t *g, const char *engine) {
  if (!g)
    return false;
  (void)gvFreeLayout(gvc, g); // ignore errors
  int err = gvLayout(gvc, g, engine);
  return err == 0;
}

// annotate the graph with layout information
bool render(Agraph_t *g) {
  if (!g)
    return false;
  attach_attrs(g);
  return true;
}

// render to stdout
bool render(Agraph_t *g, const char *format) {
  if (!g)
    return false;
  int err = gvRender(gvc, g, format, stdout);
  return err == 0;
}

// render to an open FILE
bool render(Agraph_t *g, const char *format, FILE *f) {
  if (!g)
    return false;
  int err = gvRender(gvc, g, format, f);
  return err == 0;
}

// render to an open channel
bool renderchannel(Agraph_t *g, const char *format, const char *channelname) {
  if (!g)
    return false;
  gv_channel_writer_init(gvc);
  int err = gvRender(gvc, g, format, (FILE *)channelname);
  gv_writer_reset(gvc); // Reset to default
  return err == 0;
}

// render to a filename
bool render(Agraph_t *g, const char *format, const char *filename) {
  if (!g)
    return false;
  int err = gvRenderFilename(gvc, g, format, filename);
  return err == 0;
}

typedef struct {
  char *data;
  int sz;  // buffer size
  int len; // length of array
} BA;

// render to string result, using binding-dependent gv_string_writer()
char *renderresult(Agraph_t *g, const char *format) {
  if (!g)
    return nullptr;
  if (!GD_alg(g))
    return nullptr;
  BA ba;
  ba.sz = BUFSIZ;
  // must be freed by wrapper code
  ba.data = reinterpret_cast<char *>(malloc(ba.sz * sizeof(char)));
  ba.len = 0;
  gv_string_writer_init(gvc);
  (void)gvRender(gvc, g, format, reinterpret_cast<FILE *>(&ba));
  gv_writer_reset(gvc); // Reset to default
  *reinterpret_cast<int *>(GD_alg(g)) = ba.len;
  return ba.data;
}

// render to string result, using binding-dependent gv_string_writer()
void renderresult(Agraph_t *g, const char *format, char *outdata) {
  if (!g)
    return;
  gv_string_writer_init(gvc);
  (void)gvRender(gvc, g, format, reinterpret_cast<FILE *>(outdata));
  gv_writer_reset(gvc); // Reset to default
}

// render to a malloc'ed data string, to be free'd by caller.
char *renderdata(Agraph_t *g, const char *format) {
  if (!g)
    return nullptr;
  char *data;
  unsigned int length;
  int err = gvRenderData(gvc, g, format, &data, &length);
  if (err)
    return nullptr;
  return data;
}

bool write(Agraph_t *g, FILE *f) {
  if (!g)
    return false;
  int err = agwrite(g, f);
  return err == 0;
}

bool write(Agraph_t *g, const char *filename) {
  if (!g)
    return false;
  FILE *f = fopen(filename, "w");
  if (!f)
    return false;
  int err = agwrite(g, f);
  fclose(f);
  return err == 0;
}

bool tred(Agraph_t *g) {
  if (!g)
    return false;
  int err = gvToolTred(g);
  return err == 0;
}
