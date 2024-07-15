/**
 * @file
 * @brief Graphviz context library
 * @ingroup gvc_api
 *
 * **libgvc** provides a context for applications wishing to manipulate
 * and render graphs. It provides command line parsing,
 * common rendering code, and a plugin mechanism for renderers.
 *
 * [man 3 gvc](https://graphviz.org/pdf/gvc.3.pdf)
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

#include <stdbool.h>

#include "types.h"
#include "gvplugin.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef GVDLL
#ifdef GVC_EXPORTS
#define GVC_API __declspec(dllexport)
#else
#define GVC_API __declspec(dllimport)
#endif
#endif

#ifndef GVC_API
#define GVC_API /* nothing */
#endif

/// @defgroup gvc_api Graphviz context library (GVC) API
/// @ingroup public_apis
/// @{
	
#define LAYOUT_DONE(g) (agbindrec(g, "Agraphinfo_t", 0, true) && GD_drawing(g))

/* misc */
/* FIXME - this needs eliminating or renaming */
GVC_API void gvToggle(int);

/* set up a graphviz context */
GVC_API GVC_t *gvNEWcontext(const lt_symlist_t *builtins, int demand_loading);

/*  set up a graphviz context - and init graph - retaining old API */
GVC_API GVC_t *gvContext(void);
/*  set up a graphviz context - and init graph - with builtins */
GVC_API GVC_t *gvContextPlugins(const lt_symlist_t *builtins, int demand_loading);

/* get information associated with a graphviz context */
GVC_API char **gvcInfo(GVC_t*);
GVC_API char *gvcVersion(GVC_t*);
GVC_API char *gvcBuildDate(GVC_t*);

/* parse command line args - minimally argv[0] sets layout engine */
GVC_API int gvParseArgs(GVC_t *gvc, int argc, char **argv);
GVC_API graph_t *gvNextInputGraph(GVC_t *gvc);
GVC_API graph_t *gvPluginsGraph(GVC_t *gvc);

/* Compute a layout using a specified engine */
GVC_API int gvLayout(GVC_t *gvc, graph_t *g, const char *engine);

/* Compute a layout using layout engine from command line args */
GVC_API int gvLayoutJobs(GVC_t *gvc, graph_t *g);

/* Check if a layout has been done */
GVC_API bool gvLayoutDone(graph_t *g);

/* Render layout into string attributes of the graph */
GVC_API void attach_attrs(graph_t *g);

/* Render layout in a specified format to an open FILE */
GVC_API int gvRender(GVC_t *gvc, graph_t *g, const char *format, FILE *out);

/* Render layout in a specified format to a file with the given name */
GVC_API int gvRenderFilename(GVC_t *gvc, graph_t *g, const char *format, const char *filename);

/* Render layout in a specified format to an GVC_APIal context */
GVC_API int gvRenderContext(GVC_t *gvc, graph_t *g, const char *format, void *context);

/* Render layout in a specified format to a malloc'ed string */
GVC_API int gvRenderData(GVC_t *gvc, graph_t *g, const char *format, char **result, unsigned int *length);

/* Free memory allocated and pointed to by *result in gvRenderData */
GVC_API void gvFreeRenderData (char* data);

/* Render layout according to -T and -o options found by gvParseArgs */
GVC_API int gvRenderJobs(GVC_t *gvc, graph_t *g);

/* Clean up layout data structures - layouts are not nestable (yet) */
GVC_API int gvFreeLayout(GVC_t *gvc, graph_t *g);

/* Clean up graphviz context */
GVC_API void gvFinalize(GVC_t *gvc);
GVC_API int gvFreeContext(GVC_t *gvc);

/* Return list of plugins of type kind.
 * kind would normally be "render" "layout" "textlayout" "device" "loadimage"
 * The size of the list is stored in sz.
 * The caller is responsible for freeing the storage. This involves
 * freeing each item, then the list.
 * Returns NULL on error, or if there are no plugins.
 * In the former case, sz is unchanged; in the latter, sz = 0.
 */
GVC_API char **gvPluginList(GVC_t *gvc, const char *kind, int *sz);

/** Add a library from your user application
 * @param gvc Graphviz context to add library to
 * @param lib library to add
 */
GVC_API void gvAddLibrary(GVC_t *gvc, gvplugin_library_t *lib);

/** Perform a Transitive Reduction on a graph
 * @param g  graph to be transformed.
 */
GVC_API int gvToolTred(graph_t *g);

/// @}

#undef GVC_API

#ifdef __cplusplus
}
#endif
