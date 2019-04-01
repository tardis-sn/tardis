#ifndef MBCN_PROTOS_H
#define MBCN_PROTOS_H

#include "moab/MOABConfig.h"

#if defined(MOAB_FC_FUNC_)
#define MBCN_FC_WRAPPER MOAB_FC_FUNC_
#elif defined(MOAB_FC_FUNC)
#define MBCN_FC_WRAPPER MOAB_FC_FUNC
#else
#define MBCN_FC_WRAPPER(name,NAME) name
#endif

#define MBCN_GetBasis MBCN_FC_WRAPPER( mbcn_getbasis, MBCN_GETBASIS )
#define MBCN_SetBasis MBCN_FC_WRAPPER( mbcn_setbasis, MBCN_SETBASIS )
#define MBCN_EntityTypeName MBCN_FC_WRAPPER( mbcn_entitytypename, MBCN_ENTITYTYPENAME )
#define MBCN_EntityTypeFromName MBCN_FC_WRAPPER( mbcn_entitytypefromname, MBCN_ENTITYTYPEFROMNAME )
#define MBCN_Dimension MBCN_FC_WRAPPER( mbcn_dimension, MBCN_DIMENSION )
#define MBCN_VerticesPerEntity MBCN_FC_WRAPPER( mbcn_verticesperentity, MBCN_VERTICESPERENTITY )
#define MBCN_NumSubEntities MBCN_FC_WRAPPER( mbcn_numsubentities, MBCN_NUMSUBENTITIES )
#define MBCN_SubEntityType MBCN_FC_WRAPPER( mbcn_subentitytype, MBCN_SUBENTITYTYPE )
#define MBCN_SubEntityVertexIndices MBCN_FC_WRAPPER( mbcn_subentityvertexindices, MBCN_SUBENTITYVERTEXINDICES )
#define MBCN_AdjacentSubEntities MBCN_FC_WRAPPER( mbcn_adjacentsubentities, MBCN_ADJACENTSUBENTITIES )
#define MBCN_SideNumberInt MBCN_FC_WRAPPER( mbcn_sidenumberint, MBCN_SIDENUMBERINT )
#define MBCN_SideNumberUint MBCN_FC_WRAPPER( mbcn_sidenumberuint, MBCN_SIDENUMBERUINT )
#define MBCN_SideNumberLong MBCN_FC_WRAPPER( mbcn_sidenumberlong, MBCN_SIDENUMBERLONG )
#define MBCN_SideNumberUlong MBCN_FC_WRAPPER( mbcn_sidenumberulong, MBCN_SIDENUMBERULONG )
#define MBCN_SideNumberVoid MBCN_FC_WRAPPER( mbcn_sidenumbervoid, MBCN_SIDENUMBERVOID )
#define MBCN_SideNumber MBCN_FC_WRAPPER( mbcn_sidenumber, MBCN_SIDENUMBER )
#define MBCN_OppositeSide MBCN_FC_WRAPPER( mbcn_oppositeside, MBCN_OPPOSITESIDE )
#define MBCN_ConnectivityMatchInt MBCN_FC_WRAPPER( mbcn_connectivitymatchint, MBCN_CONNECTIVITYMATCHINT )
#define MBCN_ConnectivityMatchUint MBCN_FC_WRAPPER( mbcn_connectivitymatchuint, MBCN_CONNECTIVITYMATCHUINT )
#define MBCN_ConnectivityMatchLong MBCN_FC_WRAPPER( mbcn_connectivitymatchlong, MBCN_CONNECTIVITYMATCHLONG )
#define MBCN_ConnectivityMatchUlong MBCN_FC_WRAPPER( mbcn_connectivitymatchulong, MBCN_CONNECTIVITYMATCHULONG )
#define MBCN_ConnectivityMatchVoid MBCN_FC_WRAPPER( mbcn_connectivitymatchvoid, MBCN_CONNECTIVITYMATCHVOID )
#define MBCN_setPermutation MBCN_FC_WRAPPER( mbcn_setpermutation, MBCN_SETPERMUTATION )
#define MBCN_resetPermutation MBCN_FC_WRAPPER( mbcn_resetpermutation, MBCN_RESETPERMUTATION )
#define MBCN_permuteThisInt MBCN_FC_WRAPPER( mbcn_permutethisint, MBCN_PERMUTETHISINT )
#define MBCN_permuteThisUint MBCN_FC_WRAPPER( mbcn_permutethisuint, MBCN_PERMUTETHISUINT )
#define MBCN_permuteThisLong MBCN_FC_WRAPPER( mbcn_permutethislong, MBCN_PERMUTETHISLONG )
#define MBCN_permuteThisVoid MBCN_FC_WRAPPER( mbcn_permutethisvoid, MBCN_PERMUTETHISVOID )
#define MBCN_revPermuteThisInt MBCN_FC_WRAPPER( mbcn_revpermutethisint, MBCN_REVPERMUTETHISINT )
#define MBCN_revPermuteThisUint MBCN_FC_WRAPPER( mbcn_revpermutethisuint, MBCN_REVPERMUTETHISUINT )
#define MBCN_revPermuteThisLong MBCN_FC_WRAPPER( mbcn_revpermutethislong, MBCN_REVPERMUTETHISLONG )
#define MBCN_revPermuteThisVoid MBCN_FC_WRAPPER( mbcn_revpermutethisvoid, MBCN_REVPERMUTETHISVOID )
#define MBCN_HasMidEdgeNodes MBCN_FC_WRAPPER( mbcn_hasmidedgenodes, MBCN_HASMIDEDGENODES )
#define MBCN_HasMidFaceNodes MBCN_FC_WRAPPER( mbcn_hasmidfacenodes, MBCN_HASMIDFACENODES )
#define MBCN_HasMidRegionNodes MBCN_FC_WRAPPER( mbcn_hasmidregionnodes, MBCN_HASMIDREGIONNODES )
#define MBCN_HasMidNodes MBCN_FC_WRAPPER( mbcn_hasmidnodes, MBCN_HASMIDNODES )
#define MBCN_HONodeParent MBCN_FC_WRAPPER( mbcn_honodeparent, MBCN_HONODEPARENT )
#define MBCN_HONodeIndex MBCN_FC_WRAPPER( mbcn_honodeindex, MBCN_HONODEINDEX )

#endif
