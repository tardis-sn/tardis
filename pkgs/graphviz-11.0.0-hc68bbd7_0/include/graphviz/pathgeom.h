/// @file
/// @ingroup public_apis
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

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef GVDLL
#ifdef PATHPLAN_EXPORTS
#define PATHGEOM_API __declspec(dllexport)
#else
#define PATHGEOM_API __declspec(dllimport)
#endif
#endif

#ifndef PATHGEOM_API
#define PATHGEOM_API /* nothing */
#endif

#ifdef HAVE_POINTF_S
    typedef struct pointf_s Ppoint_t;
    typedef struct pointf_s Pvector_t;
#else
    typedef struct Pxy_t {
	double x, y;
    } Pxy_t;

    typedef struct Pxy_t Ppoint_t;
    typedef struct Pxy_t Pvector_t;
#endif

    typedef struct Ppoly_t {
	Ppoint_t *ps;
	size_t pn;
    } Ppoly_t;

    typedef Ppoly_t Ppolyline_t;

    typedef struct Pedge_t {
	Ppoint_t a, b;
    } Pedge_t;

/* opaque state handle for visibility graph operations */
    typedef struct vconfig_s vconfig_t;

    PATHGEOM_API void freePath(Ppolyline_t* p);

#undef PATHGEOM_API

#ifdef __cplusplus
}
#endif
