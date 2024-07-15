# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [11.0.0] – 2024-04-28

### Added

- `gv2gml` gained a `-y` option to output the yWorks.com variant of GML instead
  of the default.
- A new command line option, `--filepath=…` has been added to perform the
  function previously served by the `$GV_FILE_PATH` environment variable, use of
  which was removed in Graphviz 6.0.1. Unlike the old `$GV_FILE_PATH` mechanism,
  `--filepath=…` takes effect regardless of the setting of the `$SERVER_NAME`
  environment variable. #2396

### Changed

- `gvpack`, in addition to recognizing a “cluster” name prefix as a mark of a
  cluster, now recognizes this case insensitively as well as recognizing the
  `cluster` attribute. This is more consistent with how the Graphviz libraries
  work.
- **Breaking**: `pkg-config` (.pc) files shipped with Graphviz now include
  `${prefix}/include` in the include path in addition to
  `${prefix}/include/graphviz`. Previously this missing path meant building
  Graphviz demo examples against an installation of Graphviz in a non-system
  path would not work. #2474
- The core PostScript output format (`-Tps`) warns if using an
  out-of-specification font name. To avoid this, use a more sophisticated output
  format like Cairo (`-Tps:cairo`) that does font name lookup and translation.
  #218
- **Breaking**: The libpack functions `putRects`, `packRects`, `putGraphs`,
  `packGraphs`, `packSubgraphs`, `pack_graph`, `shiftGraphs`, `ccomps`,
  `cccomps`, and `pccomps` now take the number of items they are operating on
  (`ng`) as a `size_t`.
- **Breaking**: The `bsearch_cmpf` and `qsort_cmpf` typedefs have been removed.
- `dot -c -v`, when constructing the config6 file, includes comments explaining
  any attempted actions that failed during plugin loading. #2456
- **Breaking**: The `Ndim` global is now a `unsigned short`.
- fdpgen no longer truncates graph names when inferring new names for connected
  component subgraphs.
- **Breaking**: The `nodequeue` type has been removed.
- **Breaking**: The field `Agraphinfo_t.n_nodes` has been removed. The function
  `agnnodes` is a more robust way of retrieving the number of nodes.
- The `-q` command line option will now suppress “no hard-coded metrics…”
  and other font lookup warnings. #2379
- **Breaking**: The `CMP` and `SGN` macros have been removed.
- The CMake build system no longer early-binds all enabled plugins into
  `dot`/`dot.exe`. This early binding was a change introduced in 10.0.1, but was
  not noted in this changelog. Traditionally, of the three Graphviz build
  systems (Autotools, CMake, MS Build), only changes to the Autotools build
  system were noted in this changelog under the assumption that packaging
  ecosystems making use of the other two build systems would need finer grained
  details and would be monitoring the Git commit history instead. This seems to
  not be the case, so in future side-effecting changes to any of the three build
  systems will be included here. #2527, #2528
- The precision of `sep`- and `esep`-based calculations has been improved.
- **Breaking**: Defines `AGRAPH`, `AGNODE`, `AGOUTEDGE`, `AGINEDGE`, and `AGEDGE` are
  replaced with `enum`.
- **Breaking**: The `obj_state_t.url_bsplinemap_poly_n` field is now a `size_t`
  and the `obj_state_t.url_bsplinemap_n` field is now a `size_t *`.
- **Breaking**: The `Ppoly_t.pn` (`Ppolyline_t.pn`) field is now a `size_t`.
- **Breaking**: The `Proutespline` function takes its `n_barriers` parameter as
  a `size_t`.
- **Breaking**: The `gvattr_t` type and the `GVJ_t.selected_obj_attributes` and
  `GVJ_t.selected_obj_type_name` fields have been removed.
- **Breaking**: The `gv_argvlist_t` type and functions that operate on it have
  been removed.
- Control characters in some error messages are escaped, preventing certain
  types of text injection that could cause user confusion.
- **Breaking**: `GVJ_t.numkeys` is a `size_t`.

### Fixed

- Indexing within `gvNextInputGraph` no longer incorrectly retains the index
  from prior use of the GVC context. When using Graphviz libraries
  programmatically, this could previously cause crashes or misbehavior. #2484
- Color schemes, typically controlled through the `colorscheme` attribute are
  now pushed and popped as they are applied and released. Previously processing
  multiple graphs wherein the first uses color schemes but later ones do not
  could result in color schemes being incorrectly retained and reapplied or
  use-after-free memory accesses.
- The GDI+ plugin, when asked to render a graphic metafile, no longer references
  uninitialized memory. This bug was introduced in Graphviz 2.24.0.
- A `free` of an invalid pointer in `edgepaint` was fixed. #2513
- `gvmap` no longer references uninitialized variables when trying to process
  triangles and encountering only 2 points.
- Using the `point` shape in combination with `peripheries=0` no longer causes
  out of bounds memory writes. This was a regression in Graphviz 7.0.0. #2497
- Unsafe use of a dangling pointer in `ccomps` has been removed. This was a
  regression in Graphviz 7.1.0.
- `gvcolor` no longer crashes when processing color names longer than 127
  characters.
- Interleaving calls to `colorxlate` and `gvrender_resolve_color` no longer
  confuse internal caching mechanisms. Callers should now get the correct color
  back.
- The `nop2` layout engine provided by the neato layout plugin is now equivalent
  to `neato -n2` as intended instead of mistakenly being equivalent to
  `nop`/`nop1`/`neato -n1`.
- An off-by-one error in rank installation was corrected. Previously, an unusual
  `rank=same` constraint could cause a crash when installing ranks. #1308
- `gxl2gv` no longer crashes or misbehaves when symlinked to a non-ASCII file
  name. This is a rare scenario that normal users should not encounter.
- `mm2gv` no longer crashes or misbehaves when reading malformed Matrix Market
  files with non-ASCII bytes in the header.
- A stack buffer overflow in `mm2gv` when processing malformed Matrix Market
  files has been fixed.
- The `newrank` attribute is treated as a boolean instead of any value
  (including `"false"`) being coerced into `"true"`. #2521
- Crashes and misbehavior no longer occur when the `sides` attribute contains
  non-ASCII characters.
- Graphviz binaries like `dot.exe` and `neato.exe` no longer crash or misbehave
  when symlinked to a non-ASCII file name on Windows. This is a rare scenario
  that normal users should not encounter.
- GVPR programs that use `tolower` or `toupper` on strings containing non-ASCII
  characters no longer crash. These functions do not lowercase/uppercase
  non-ASCII characters, so users probably still do not want to use non-ASCII
  strings in a GVPR program.
- Some `routesplines` miscalculations that led to lost edges and fatal errors
  have been avoided. #2368
- An inaccuracy involving an edge case when constructing lines within libpack
  has been corrected.
- A bug in the internal heap implementation used in the network simplex
  algorithm has been corrected. This would previously cause certain runs to
  infer incorrect ordering or subtrees. This was a regression in Graphviz
  2.40.0. #2391, #2529
- Compass points may be more accurately placed on the node boundary in some cases.
- A very small random adjustment in the calculation of the space available for
  edge routing around ellipse shaped nodes in fdp and neato layouts, has been
  removed.
- Incorrect edge splines for ellipse shaped nodes with ports using fdp or
  neato. #2168
- Incorrect edge splines for ellipse and polygon shaped nodes with ports and
  large penwidths using fdp or neato, causing the same symptoms as #2168.
- Incorrect edge splines for polygon shaped nodes with ports more than one
  periphery using fdp or neato, causing the same symptoms as #2168.
- Adjust the space available for edge routing based on penwidth when
  using fdp or neato and `splines=ortho`.

## [10.0.1] – 2024-02-11

### Added

- Releases now include packages for [Rocky Linux](https://rockylinux.org/) 8 and
  9.
- A new output format, `-Tsvg_inline`, has been added to generate a header-less
  SVG suitable for inlining into HTML. #2285
- The functionality of the `acyclic`, `tred` and `unflatten` command line tools
  are now exposed via the `graphviz_acyclic`, `graphviz_tred` and
  `graphviz_unflatten` API functions in libcgraph. #2194
- `graphviz_node_induce` is available as a new API function in cgraph.h.
- `tred` gained a `-o` command line option to redirect its output to a file.

### Changed

- The Criterion unit tests have been removed and migrated to Pytest. This is
  primarily relevant to downstream packagers of Graphviz. #2443
- **Breaking**: `Dtdisc_t.memoryf` and its associated macros has been removed.
- **Breaking**: The `Dt_t.type` field has been removed.
- **Breaking**: The `dtfound`, `DT_FOUND`, `dtleast`, and `dtmost` macros have
  been removed.
- The nrtmain.c test program has been removed from the portable tarball.
- The TCL Graphviz packages for inter-release versions/snapshots report
  themselves as `<next release>b<internal number>` instead of
  `<next release>~dev.<internal number>`. This fixes a problem wherein TCL would
  see `~dev` as being invalid characters to appear in a version. #2370
- Support for discovering Lua via `lua-config*` has been removed from the
  Autotools build system.
- Lua discovery in the Autotools build system should now respect the location of
  your Lua installation and not unconditionally attempt installation into
  `/usr`. #2152
- The GTK plugin is no longer built or distributed. This plugin relies on GTK 2
  and X11. If you use this plugin, please contact the maintainers to let
  them know it is worthwhile re-enabling this and forward porting it to GTK 3/4
  and Wayland. #1848
- In the Autotools build system, `LIBPOSTFIX=` can now be used to suppress `64`
  being appended to the library installation path.
- The `-m` command line option, whose functionality was disabled in Graphviz
  3.0.0, has been removed.
- Man page typography has been slightly improved.
- macOS release artifacts no longer include `vimdot`. This may be restored in
  future. #2423
- macOS release artifacts no longer include `smyrna`. This may be restored in
  future. #2422
- The PDF output format, `-Tpdf`, respects the environment variable
  `$SOURCE_DATE_EPOCH` for overriding `CreationDate` when built against Cairo
  ≥ 1.16.0. #2473
- The legacy C# viewer app is no longer distributed in the portable source
  tarball.
- Graphviz headers no longer define the `FALSE` and `TRUE` constants.
- The Autotools build system no longer supports Darwin 9 (Mac OSX Leopard).
- **Breaking**: `Agraph_t.link` has been split into `Agraph_t.id_link` and
  `Agraph_t.seq_link`. `Agraph_t.g_dict` has been split into `Agraph_t.g_id`
  and `Agraph_t.g_seq`.
- **Breaking**: `gvpropts.n_outgraphs` is now a `size_t`.
- The OCaml bindings have been removed. If you use these bindings, please contact
  the maintainers to notify them of the existence of users.
- **Breaking**: `polygon_t.sides` and `polygon_t.peripheries` are now `size_t`s.
- **Breaking**: liblab_gamut is no longer included in a Graphviz installation.
  This library had no accompanying header, so using it was not easy. If you are
  using this library, please contact the maintainers to notify them of the
  existence of users. #2489
- **Breaking**: `bezier.size` and `splines.size` are now `size_t`s.
- **Breaking**: the gv.i and gv.cpp SWIG inputs are no longer included in a
  Graphviz installation. #2491
- **Breaking**: the `gvrender_engine_t.beziercurve`,
  `gvrender_engine_t.library_shape`, `gvrender_engine_t.polygon`,  and
  `gvrender_engine_t.polyline` callbacks now take the number of points, `n`, as
  a `size_t`.
- **Breaking**: the `AVG` macro has been removed.
- **Breaking**: the `inside_t.s` union member gained members `lastn`, `radius`,
  `last_poly`, `last`, `outp`, `scalex`, `scaley`, `box_URx`, and `box_URy`.
  Zero initialize these when you construct instances of this type. #2498

### Fixed

- The paper size for Doxygen docs generation in the Autotools build system has
  been corrected to `a4`.
- References to `eventf` and `hashf` data structures in the libcdt man page
  have been removed. These data structures were removed in Graphviz 9.0.0.
- References to `DTOFFSET` in the libcdt man page have been removed. This macro
  was removed in Graphviz 2.40.0.
- A number of further updates to the libcdt man page have been made to reflect
  other changes that happened in Graphviz 9.0.0.
- Use of the non-portable `PATH_MAX` constant has been removed. This was a
  regression in Graphviz 7.0.1. In addition to fixing the regression, code has
  been adjusted to remove assumptions on the maximum path length and treat it as
  unbounded. #2452
- Compilation on NetBSD has been repaired. This was a regression in Graphviz
  9.0.0.
- Compilation on SunOS has been repaired. This appears to have been broken since
  the xlib plugin was added some time prior to Graphviz 2.38.0.
- gvpr programs that attempt to close out of range file descriptors no longer
  cause out of bounds memory accesses.
- When large edge weights are used that cause an integer overflow when summing
  them up, Graphviz now gracefully exits with an error message instead of
  crashing. #2450
- Support for the `%n` specifier in `scanf` in gvpr has been restored. This was
  a regression in Graphviz 9.0.0. #2454
- In the Autotools build system, `make dist` now processes cmd/gvedit correctly
  when Qt is not installed. Generating Qt “mocables” is postponed from configure
  time to build time. #2463
- The Autotools build system correctly detects Ruby headers, even when
  pkg-config support is unavailable. #2464
- Escaped characters in xdot fields no longer lead to the containing text being
  truncated. #2460
- When building against a libgd that is configured with `!gif && (jpeg || png)`,
  the GD plugin is once again compilable. This was a regression in Graphviz
  2.46.0.
- edgepaint spline intersection code would previously incorrectly use the second
  spline in one instance where it should have used the first. #1464
- In the Autotools build, libexpat discovery on macOS has been improved. #2477
- A bug that caused compound edges to sometimes be drawn in the wrong direction
  has been corrected. This was a regression in Graphviz 8.0.3. #2478
- When operating on multiple graphs, `unflatten` no longer retains chain node
  and size internal state across graphs.
- Repeated runs of a graph with subgraphs now produce a stable subgraph
  ordering. #2242
- The `dot` and `gml2gv` tools are now built with case-insensitive parsing
  by the CMake and MSBuild systems, as they always were by autotools, and
  in accordance with the graphviz specification. #2481
- Putting nodes in a subgraph no longer causes their layout order to be
  reversed. #1585
- Edges are no longer lost when using subgraphs and record shapes in
  combination. #1624
- A malformed config6 file that leads to plugin search failing no longer causes
  out-of-bounds memory reads. This now causes an error message and graceful
  failure. #2441
- Discovery of `php` in the Autotools build system has been improved.
- Text in the PIC output format is no longer forced to font size 1. This was a
  regression in Graphviz 8.0.2. Even with this fix, the PIC output format is
  limited in its utility. #2487
- When encountering a syntactically invalid HTML-like label, Graphviz.app no
  longer aborts. The abort was an intentional change in Graphviz 8.0.1 to avoid
  invalid memory reads in `dot`, but had the undesirable side effect of the
  graphical Graphviz.app exiting with no obvious cause. #2488
- Use of an uninitialized variable in `poly_inside` has been corrected. #2498
- Input containing UTF-8 data that is destined to appear as-is in the output
  (e.g. UTF-8 characters in a label when using the `-Tdot` output format) is
  once again processed correctly. On platforms with a signed `char` this could
  previously crash. This was a regression in Graphviz 2.49.0. #2502

## [9.0.0] - 2023-09-11

### Added

- On non-Windows platforms, new `-Tkitty` and `-Tkittyz` output formats are
  available that render to the Kitty terminal emulator’s graphvics protocol.
- HTML/CSS-style 3 letter hex colors are supported. Each R/G/B letter is
  duplicated to form a 6 letter hex color. E.g. `#09c` is equivalent to
  `#0099cc`. #2377

### Changed

- **Breaking**: The definition of `adjmatrix_t` is no longer exposed in public
  headers.
- **Breaking**: The upper limit for minimum edge length (`Agedgeinfo_t.minlen`)
  has been expanded from `USHRT_MAX` to `INT_MAX`. #2413
- **Breaking**: The libcdt macros `DTTREEMATCH`, `DTTREESEARCH`, `dtvnext`,
  `dtvcount`, `dtvhere`, and `dtcharhash` have been removed.
- **Breaking**: The libcgraph macros `AGHEADPOINTER`, `AGRIGHTPOINTER`,
  `AGLEFTPOINTER`, `FIRSTNREF`, `NEXTNREF`, `PREVNREF`, `LASTNREF`, `NODEOF`,
  `FIRSTOUTREF`, `LASTOUTREF`, `FIRSTINREF`, `NEXTEREF`, and `PREVEREF` have
  been removed.
- **Breaking**: The libcgraph types `Agnoderef_t` and `Agedgeref_t` have been
  removed.
- **Breaking**: The libcgraph function `agflatten` has been removed.
- **Breaking**: The `Agdesc_s.flatlock` field has been removed.
- **Breaking**: The `str` parameter from `gvPluginList` has been removed.
- **Breaking**: The definition of the `elist_append` and `alloc_elist` macros
  have been changed to use newer allocation functions. Users were/are not
  expected to call these macros.
- The functions `ageqedge`, `agtail`, `aghead`, `agopp`, `agmkout`, and `agmkin`
  have been reintroduced. These were previously removed in Graphviz 3.0.0. #2433
- **Breaking**: The first parameter `dt` to the `makef` and `freef` callbacks
  defined in cdt.h has been removed.
- Gvedit now identifies itself with organization name “Graphviz” and application
  name “gvedit” when reading and writing Qt-based settings. It would previously
  use organization name “Trolltech” and application name “MDI Example”. If you
  have existing settings under the old identification, Gvedit will attempt to
  migrate them to the new identification the first time it reads then writes
  settings. #2383
- **Breaking**: `gvprintf` is now tagged with
  `__attribute__((format(printf, …)))` when compiling with Clang or GCC. This
  enables the compiler to spot more misuses of this function. #2373
- **Breaking**: The `hashf` and `eventf` members of `Dtdisc_t` have been
  removed. Correspondingly, the `hshf` and `evf` parameters to the `DTDISC`
  macro have been removed. Also the `_DTHSH` macro has been removed.
- **Breaking**: The `Dtdata_t.minp` field has been removed.
- The print functionality of the macOS Graphviz app scales the graph to fit the
  output page size.
- **Breaking**: The libcdt containers `Dtbag`, `Dthash`, `Dtlist`, `Dtorder`,
  `Dtdeque`, and `Dtstack` have been removed.
- **Breaking**: The libcdt `dtappend` and `dtattach`  macros have been removed.
- Support for Lua 5.0 has been removed. Building the Graphviz Lua bindings now
  requires Lua ≥ 5.1.
- **Breaking**: The `Dt_t*` parameter to the callback for `dtwalk` has been
  removed.
- **Breaking**: The `POINTS_PER_PC` macro has been removed.
- **Breaking**: The `INITIAL_XDOT_CAPACITY` macro has been removed.
- **Breaking**: The `type` parameter to `dtdisc` has been removed.
- **Breaking**: The `h` parameter to `dtstrhash` has been removed.
- In addition to Guile 2.0 and Guile 2.2, Guile 3.0 is now supported by the
  Graphviz Guile bindings.
- **Breaking**: The concept of “memory allocator discipline” has been removed,
  along with the type `Agmemdisc_t` and fields `Agdisc_t.mem` and
  `Agdstate_t.mem`.
- **Breaking**: The `agcallbacks` function and `Agclos_t.callbacks_enabled` have
  been removed.
- **Breaking**: `pack_info.doSplines` is now a C99 `bool`. Correspondingly, the
  `doSplines` parameter to `shiftGraphs` is now a C99 `bool`.

### Fixed

- Processing large graphs that induce ranks containing more than 46340
  (`floor(√INT_MAX)`) nodes no longer results in an integer overflow during
  crossing matrix allocation. Ranks of up to `floor(√SIZE_MAX)` nodes are now
  supported.
- Double arrow head types like `invdot` and `onormalonormal` once again display
  correctly. This was a regression in Graphviz 8.0.1. #2406
- The `lvee` and `rvee` edge arrow shapes are slighty incorrect for
  penwidths &gt; 1. #2399
- Small gap between `lcurve` or `rcurve` arrow shaft and node. #2426
- Failure of arrowhead and arrowtail to respect penwidth #372 \
  Fixed also for the `normal` and `inv`
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html)
  when using the `l` or `r`
  [arrow shape modifiers](https://graphviz.org/doc/info/arrows.html#shape-modifiers). \
  Slightly improved for the `normal` and `inv`
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html)
  when not using any
  [arrow shape modifier](https://graphviz.org/doc/info/arrows.html#shape-modifiers). \
  Fixed also for the `crow` and `vee`
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html#primitive-shapes)
  and [record based nodes](https://graphviz.org/doc/info/shapes.html#record).
- Various incorrect calls to `gvprintf` have been corrected. On some platforms
  like x86-64, these problems were benign thanks to coincidences in the
  Application Binary Interface (ABI). On other platforms, these problems may
  have caused stack corruption and crashes. #2373
- The font dictionary is now initialized even if a configuration file is not
  found. Previously this situation could lead to a crash when using Graphviz
  programmatically. This problem was present as least as far back as Graphviz
  2.38.0. #1520
- **Breaking**: The `vt100` output format has been renamed to `vt`. This fixes a
  problem where it was not discoverable on macOS. #2429
- Escape sequences like `\"` are recognized in strings and double escaping
  (producing `\\"`) is avoided. #2397
- The Autotools build system no longer uses headers and libraries from the
  `--prefix` path given on the command line. This previously caused
  cross-compilation to incorrectly pick up host headers and libraries. #2442

## [8.1.0] – 2023-07-06

### Added

- On non-Windows platforms, new `-Tvt100` and `-Tvt100-24bit` output formats are
  available that do rudimentary rendering to a terminal that supports ANSI
  escape sequences.
- Some notes about the interaction with wide-oriented streams were added to the
  cgraph man page.

### Changed

- When memory allocation failures cause Graphviz to exit, information about the
  failing allocation is included in the error message.

### Fixed

- Failure of arrowhead and arrowtail to respect penwidth #372 \
  Fixed also for the `curve` and `icurve`
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html#primitive-shapes).
- Iteration calculations based on `nslimit` and/or `nslimit1` attributes are
  clamped to `[0, INT_MAX]`. That is, calculations that result in a negative
  number of iterations are rounded up to `0` and those that result in a number
  of iterations that exceeds `INT_MAX` are rounded down to `INT_MAX`. Iteration
  numbers outside this range do not have useful behavior, but could be caused
  unintentionally by users.
- Setting `xdotversion` attributes that parse as greater than 17 (`xdotversion`
  is parsed as a series of digits, ignoring all non-digits) no longer causes an
  out of bounds read when rendering to xdot. #2390
- Icon size in the macOS Graphviz.app has been fixed so icons are no longer
  invisible.
- Compiling the portable source tarball on a machine with SWIG ≥ 4.1.0 no
  longer fails due to missing PHP files. #2388
- Syntax for a loop in `gvmap.sh` has been corrected. This was a regression in
  Graphviz 2.50.0. #2404

## [8.0.5] – 2023-04-30

### Changed

- Support for versions of Pango prior to 1.22.0 has been removed.
- On Windows, the Pango plugin now uses the newer `pango_layout_get_baseline`
  API.
- `dot` no longer installs a signal handler for `SIGINT`. This means typing
  Ctrl+C while `dot` is running will no longer attempt a partial render and exit
  with 0 status. Ctrl+C will have the standard behavior, typically aborting
  `dot`.

### Fixed

- A minor inaccuracy for some cases of clipping an edge to a polygon node
  boundary has been fixed.
- A minor inaccuracy in node width and height for some cases of rendering
  polygon nodes has been fixed.
- A minor inaccuracy for some cases of calculating text height from `fontsize`
  in the GD plugin has been fixed.
- A minor vertical misalignment of text in the GD plugin has been fixed.
- Instead of using the actual font height as given by the font metrics, an
  approximation based on font size was used in the Pango plugin.
- A minor inaccuracy for some cases of calculating text width and height in the
  Pango plugin has been fixed.
- A minor vertical misalignment of text in the Pango plugin has been fixed.
- Ensure `HAVE_PANGOCAIRO` is set when using CMake and the library is available.
- A minor inaccuracy in node width and height for some cases of defining
  polygon-based nodes has been fixed.
- A minor inaccuracy for some cases of calculating margin for record-based nodes
  has been fixed.
- A minor inaccuracy in node width and height for some cases of defining
  record-based nodes has been fixed.
- On all known supported platforms except 32-bit MinGW, graphs involving small
  distance constraints no longer cause a crash during majorization. #1554

## [8.0.3] – 2023-04-16

### Added

- A pre-processor script for resolving external image references
  (`image="http…"`) is now available. This enables graphs to reference images
  from intranet or internet locations. #1664

### Changed

- The accuracy of box overlapping checks during routing has been improved.

### Fixed

- makeCompoundEdge: Assertion `bez->sflag` failed. #1879
- Graphviz.app’s export functionality has been restored. This was a regression
  in Graphviz 5.0.0. #2375

## [8.0.2] – 2023-04-10

### Changed

- The Autotools build system can now detect a MacPorts-installed libANN. #1854
- Support for versions of Cairo prior to 1.10 has been removed.
- Graphs that generate rectangles with areas in the range [2³², 2⁶⁴ - 1]
  are now supported. Previously areas greater than 2³² - 1 would be
  rejected. #2371

### Fixed

- Head and tail of `digraph` edges with `dir = both` were inverted if
  `splines = ortho` was used. The bug was only exposed on straight edges.
  Edges with at least one corner were unaffected. #144
- `_Gdtclft_Init` link errors when builting libtcldot_builtin using the
  Autotools build system have been resolved. #2365
- Incorrect string construction in the `pov` output formatter has been fixed.

## [8.0.1] – 2023-03-27

### Added

- When specifying a color in HSV format, it is now possible to give an
  additional fourth component for the alpha channel (opacity). #510

### Changed

- Graphviz will now exit when encountering a syntactically invalid HTML label
  instead of attempting to recover and continue. #1311
- **Breaking**: the `url_map_n` field in the `obj_state_t` struct is now a
  `size_t`.
- The limit of 5 unique `samehead` or `sametail` values per node has been
  removed. The maximum number of unique `samehead` or `sametail` values is now
  limited only by available memory. #452
- **Breaking**: The `size` field of the `elist` struct is now a `size_t`.
- **Breaking**: The `size` field of the `nlist` struct is now a `size_t`.
- **Breaking**: The `n_nodes` field of the `Agraphinfo_t` struct is now a
  `size_t`.
- **Breaking**: The `nspans` field of `textlabel_t.u.txt` is now a `size_t`.
- **Breaking**: The `sflag` and `eflag` fields of the `bezier` struct are now
  `uint32_t`s.
- **Breaking**: The `nvertices` field of the `stroke_t` struct is now a
  `size_t`.
- “no hard-coded metrics” warnings are now only printed once per font.
- The Autotools build system now discovers Python paths using `sysconfig`
  instead of `distutils.sysconfig`, making it compatible with Python 3.12. This
  alters the installation path of the Python Graphviz module to something more
  correct. #2332

### Fixed

- The `pic` output renderer uses PIC or troff comments where appropriate, fixing
  a problem that resulted in comments being misinterpreted by `groff` and being
  visible in the final output. #2341
- `gv2gxl` and `gxl2gv` correctly detect their mode (gv→gxl or gxl→gv) on
  Windows when called via an absolute or relative path. #2344
- Using C pre-processor line directives (`^\s*#(line )?\d+`) claiming a line
  number greater than `INT_MAX` no longer causes an integer overflow. #1318
- fdp cluster→cluster edges were correct but now drawn incorrectly. This was a
  regression in Graphviz 7.0.0. #2345
- Failure of arrowhead and arrowtail to respect penwidth #372 \
  Fixed also for the `cylinder`
  [node shape](https://graphviz.org/doc/info/shapes.html#polygon).
- Second periphery of a cylinder shaped node is not correct. #2297
- Graphs with more than 127 layers no longer cause out of bound writes. #2355
- htmltable.c assertions are no longer triggered by using HTML table cells too
  small to fit their content. #1622
- `dot2gxl -d` no longer crashes when receiving a closing `node` tag following a
  closing `graph` tag. #2094
- A buffer overflow in Smyrna when loading GVPR scripts has been corrected.
- A buffer overflow when loading a plugin with a long type string has been
  corrected.
- Graphs that involve more than 2000 stroke points during tapering calculations
  no longer cause out of bounds writes.
- Using `arrowsize=0.0` no longer triggers an assertion failure or crash during
  miter calculations. This was a regression in Graphviz 7.0.0. #2342
- When using the `beautify=true` attribute, beautification no longer confuses
  leaves and dimensions. This previously could have led to skipping calculations
  or infinite loops.
- When using the `beautify=true` attribute, the first and last nodes around a
  circular layout are no longer placed on top of each other. #2283
- Applying `concentrate=true` to duplicate edges no longer results in errors due
  to non-normal edges being found. #2087
- `splines=ortho` and `concentrate=true` when used in combination no longer
  causes crashes during spline construction. #2361
- Externally referenced SVG files with their opening `<svg` tag on the same line
  as their XML declaration are no longer ignored. #2352

### Removed

- The VML output renderer has been removed. This format has been superseded by
  SVG. #799
- Legacy man page references to `dotty` have been removed. `dotty` was removed
  in Graphviz 4.0.0.
- **Breaking**: The definition of the `elist_fastapp` macro has been removed.
- Versions of Librsvg prior to 2.36.0 are no longer supported.
- Versions of GDK prior to 2.0.0 are no longer supported.
- Versions of Glib prior to 2.36.0 are no longer supported.
- **Breaking**: The `Agnodeinfo_t.inleaf` field and its `ND_inleaf` accessor
  have been removed.
- **Breaking**: The `Agnodeinfo_t.outleaf` field and its `ND_outleaf` and
  `GD_outleaf` accessors have been removed.
- **Breaking**: The `Agraphinfo_t.has_sourcerank` field and its
  `GD_has_sourcerank` accessor has been removed.
- **Breaking**: The `Agraphinfo_t.has_sinkrank` field and its
  `GD_has_sinkrank` accessor has been removed.
- Support for the legacy Microsoft Visio VDX format has been removed.
- **Breaking**: The `arrow_at_start` and `arrow_at_end` parameters from the
  `gvrender_engine_t.beziercurve` callback have been removed.
- **Breaking**: The `GVRENDER_DOES_ARROWS` constant has been removed.
- The extra cmpnd.c code is no longer shipped in the Graphviz distribution
  tarball.
- **Breaking**: The `STROKE_CLOSED`, `STROKE_FILLED`, `STROKE_PENDOWN`, and
  `STROKE_VERTICES_ALLOCATED` constants have been removed.
- **Breaking**: The `stoke_t.flags` field has been removed.

## [7.1.0] – 2023-01-21

### Added

- The command line option `--help` has been added as an alias for `-?`. #1618
- The command line option `--version` has been added as an alias for `-V`. #1618

### Fixed

- The Autotools build system no longer errors when attempting libANN discovery
  during cross-compilation. This was a regression in Graphviz 7.0.6. #2335
- Graphs with more than 46341 (⌈√INT_MAX⌉) nodes no longer crash `twopi`. #1999
- Compatibility with `/bin/sh` has been restored in the Autotools build system.
  This was a regression in Graphviz 7.0.6. This restores the ability to compile
  on NetBSD which was fixed in 7.0.4 but regressed in 7.0.6. #2340
- `ccomps` no longer crashes when failing to open files.

## [7.0.6] – 2023-01-06

### Changed

- The Autotools build system no longer looks for `python` binaries. The Python
  interpreter is unconditionally assumed to be `python3`. The configure option
  `--enable-python` is now an alias for `--enable-python3`.
- The Autotools and CMake build systems, when building `gvedit`, will now look
  for and use Qt6 in preference over Qt5. #2233
- Reserved stack size on Windows for the `dot.exe` binary has been increased
  from ~3.8MB to 32MB. #1710
- Reserved stack size on macOS for the `dot` binary when built with CMake has
  been increased from 8MB to 32MB. #1710
- The Autotools build system can now find libANN on Debian-based operating
  systems, enabling compilation of libmingle and `mingle`. #1835
- The `webdot` web server interface to Graphviz has been removed. For a modern
  replacement, see for example https://github.com/magjac/d3-graphviz. #1131

### Fixed

- The modgraph.php example no longer includes gv.php, which is no longer
  generated by SWIG 4.1.0. #2322

## [7.0.5] – 2022-12-23

### Fixed

- Using `style` attributes in excess of 128 bytes and/or 63 individual styles no
  longer results in out-of-bounds memory accesses. #2325

## [7.0.4] – 2022-12-03

### Fixed

- The `alt` attributes are once again set in the cmap output. This was a
  regression in Graphviz 7.0.2, that intentionally removed these but did not
  account for the W3C specification making these attributes required when the
  `href` attribute is set. #265, #2319
- Building Graphviz from source using the Autotools build system in now possible
  on NetBSD. #2317
- The ortho library now allocates trapezoid structures on-demand, removing the
  “Trapezoid-table overflow” error that previously occurred when its upfront
  estimation was exceeded. #56, #1880

## [7.0.3] – 2022-11-26

### Changed

* Support for the Intel C Compiler in the Autotools build system has been
  removed. #2298
* Fallback typedefs for `ssize_t` have been removed from the CMake build system.

### Fixed

- The CMake build system no longer builds auxiliary tools beyond `gvpack` with
  demand loading disabled.
- `gvpack` built with the CMake build system can now find plugins correctly at
  run time. #1838

## [7.0.2] – 2022-11-18

### Added

- The `cluster`, `dot_builtins`, `dot2gxl`, `gv2gxl`, `gvedit`, and `prune`
  utilities are now included in the CMake build system. #1753, #1836

### Changed

- `gvedit` now uses a relative path from its own executable to discover its
  attributes file, `../share/graphviz/gvedit/attrs.txt`. This should make it
  more relocatable and make its behavior more consistent across operating
  systems.
- `alt` tags are no longer set in the cmap output. #265

### Fixed

- `gxl2gv`, when dealing with `name` attributes, may be less likely to crash. We
  say “may be less likely” because a bug remains that makes a crash still
  the most likely outcome. #2300
- Gradient URL references in SVG output once again align with their targets
  (linear or radial gradients) when `id` attributes are in use. This was
  arguably a regression in Graphviz 6.0.1. #2307
- The CMake build system’s interactions with a Zlib installed in a non-system
  location has been improved.
- Do not try to install `gv.php` if using SWIG-4.1.0. Graphviz 7.0.1 changes
  listed SWIG 4.1.0 as supported, but there was an oversight that is fixed in
  7.0.2. Complete #2277, #2303
- Several compilation errors when building Smyrna on macOS have been fixed. This
  was a regression in Graphviz 7.0.1.
- A crash when using neato layout with large inferred per-node edge counts was
  fixed. #42

## [7.0.1] – 2022-11-09

### Added

- SWIG 4.1.0 is now supported in the Autotools build system. #2277, #2303

### Changed

- When built with zlib support, Graphviz will unconditionally use
  `deflateBound`. The only user-visible effect should be slightly decreased
  memory usage when using a zlib-compressed output format.
- The test suite only detects Graphviz companion programs adjacent to the first
  `dot` found in `$PATH` #2201

### Fixed

- Failure of arrowhead and arrowtail to respect penwidth #372 \
  Fixed also for the `diamond` and `tee`
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html#primitive-shapes).
- The CMake build system no longer uses the final install location as the
  staging area for example graphs and templates during `cpack`. This bug was
  introduced in Graphviz 4.0.0. #2232
- The CMake build system uses corrected paths to DevIL and FreeType headers and
  libraries when discovered.
- The CMake build system under MinGW no longer attempts to install third party
  Expat and Getopt libraries.

## [7.0.0] – 2022-10-22

### Changed

- **Breaking**: An `outline_width` field has been added to the `Agnodeinfo_t`
  struct.
- **Breaking**: An `outline_height` field has been added to the `Agnodeinfo_t`
  struct.
- When using the CMake build system, the minimum requirement has been increased
  from CMake 3.9 to CMake 3.13.0.
- When compiling Graphviz with the Intel C Compiler, the Autotools build system
  no longer suppresses `-ffast-math`. Users relying on this are encouraged to
  investigate what in their build environment is appending a flag their C
  compiler does not support.
- The `-ffast-math` compiler flag is no longer enabled by the Autotools build
  system.
- Reserved stack size on Windows for the `dot.exe` binary has been increased
  from the default 1MB to ~3.8MB.

### Fixed

- Issues with GCC 8, `-O3` and `-ffast-math` #1412
- Problem building Graphviz-5.0.1: undefined symbols `__*_finite` #2296
- Failure of arrowhead and arrowtail to respect penwidth #372 \
  Fixed for all
  [polygon-based node shapes](https://graphviz.org/doc/info/shapes.html#polygon)
  (except for the `cylinder` node shape) and for the
  [edge arrow shapes](https://graphviz.org/doc/info/arrows.html)
  `normal`, `inv`, `box` and `dot`

### Removed

- Support for the MOSEK commercial solver has been removed. The `MOSEK`
  build-time macro no longer does anything.

## [6.0.2] - 2022-10-11

### Fixed

- Using `aspect` with a custom value in the `packmode` attribute is no longer
  misparsed.
- Smyrna bounding box computation has been corrected. There was a regression in
  4.0.0 that resulted in a degenerate bounding box being computed for any input
  graph. See #2279 for details.
- Smyrna warnings about the `on_attrRB0_clicked` and `on_attrSearchBtn_clicked`
  handlers being not found have been fixed and the buttons to which they are
  wired should be functional again. This was a regression in 2.50.0 See #2279
  for details.
- Smyrna warnings about the `mAttributesSlot` handler being not found have been
  fixed and the button to which it is wired should be functional again. This was
  a regression in 2.49.2 See #2279 for details.
- Graphviz no longer fails to load private Ghostscript symbols ("Could not load
  `libgvplugin_gs.so.6`) #2280
- trailing commas issue with fdp layout #2282
- Missing `-S` in `sccmap` man page usage summary.
- In `sccmap`, a `-v` option following a `-S` option now re-enables strongly
  connected component output that the man page implied.

## [6.0.1] – 2022-09-11

### Changed

- **Breaking**: libxdot fields for the size and number of operations, the
  statistics counts, and polygon line points are now `size_t` values instead of
  `int` values
- Accuracy of the bounding boxes printed by the `showboxes` feature have been
  improved.

### Fixed

- Id attribute is not used in linearGradient. #2258
- Graphviz 5.0.1 undocumented change of automatically generated output filename
  with -O flag (missing dot separator). This was a regression in 5.0.1. #2270
- Assert fail in `aaglex` for multiple calls to `agmemread`. This was a
  regression in 5.0.1. #2272

### Removed

- The `$GV_FILE_PATH` sandboxing mechanism has been removed #2257

## [5.0.1] – 2022-08-20

### Fixed

- -Tx11: Assertion `xev.xbutton.button >= 1 && xev.xbutton.button <= 5 && "Xlib
  returned invalid button event"` failed #2256
- missing Perl includes patch #2262
- smyrna: incorrect tokenization in frmobjectui.c:set_attr_object_type #2259
- [Dot] Record shape+UTF+" | "=Eats spaces. #925
- Memory leak in osage
- Segmentation fault when running test example neatopack.c #1800
- Memory leak in dot when using clusters
- Memory leak in patchwork when using clusters
- Subgraph layout and rendering
- Generated SVG files no longer use `transparent` paint or color as this keyword
  does not exist in SVG 1.1, and instead use `none` or a zero-opacity color.
- Unnecessary space in 'points' attribute for 'polyline' in SVG output
- `graphml2gv` no longer calls itself `cvtgxl` in error messages

### Added
 - GVContext::version() to lib/gvc++
 - GVContext::buildDate() to lib/gvc++

## [5.0.0] – 2022-07-07

### Changed

- `use_sanitizers` option has been removed from the CMake build system.

### Fixed

- **Breaking**: The 4.0.0 change replacing the `Agiodisc_t` struct member
  `putstr` by `printf` has been reverted
- graphviz-4.0.0: build error: cmd/tools/gvcolor.c:159: undefined reference to
  `fmax` #2246
- Failed assertion in `chkSgraph` for twopi layout and ortho splines. #14
- Failed assertion in `chkSgraph` for dot layout and ortho splines. #1408
- Failed assertion in `chkSgraph` for circo layout and ortho splines. #1990
- Segmentation Fault with splines="ortho". #1658
- Transparent Label appear in SVG output #146
- Binary tcl modules should compile with -module #1285
- b15.gv crashes dot #827
- heap overflow in function startElementHandler in gxl2gv.c #2093
- Crash on assertion #121
- `xdotversion` attribute is no longer misparsed. This was a regression in
  Graphviz 2.47.2. #358

## [4.0.0] – 2022-05-29

### Changed

- **Breaking**: The `mark` field of the `Agnodeinfo_t` struct is now a
  `size_t` instead of a `char`.
- **Breaking**: The unused `shape_t` struct has been removed from the public
  header `types.h`
- **Breaking**: The `Agiodisc_t` struct member `putstr` that was previously an
  `fputs` analog is replaced by `printf` that is required to behave similar to
  `fprintf`.
- the `mingle`, `diffimg`, `gvmap`, and `edgepaint` binaries are now included in
  the CMake build system
- the `gvmap.sh` and `vimdot` scripts are now installed by the CMake build
  system on operating systems other than Windows
- a brief note about the (previously undocumented) behavior of Graphviz when
  sent `SIGUSR1` is now mentioned in the man page
- build system support for `dotty`, `lefty`, and `lneato` has been removed
- the CMake build system now includes the DevIL, GDK, GhostScript, GTK, LASi,
  Poppler, Quartz, Rsvg, Visio, WebP, and Xlib plugins
- `awk` is no longer a build-time dependency #2118

### Fixed

- `agcanon`, `agcanonStr`, and `agwrite` now return error values on memory
  allocation failures instead of crashing or corrupting data
- `gvpr` programs can now pass dynamically allocated arguments to user-defined
  functions without corrupting their content. Some cases of this were a
  regression in Graphviz 2.46.0. Other cases have existed since the first
  release of `gvpr`. #2185
- spurious "no hard-coded metrics" warnings on labels with empty lines #2179
- fixed corruption of user shape characteristics during EPSF initialization
- output formats canon, dot, and xdot are not completely faithful to input #2184
- gvpr index function produces wrong results #2211. This was a regression in
  Graphviz 2.47.0.
- Error on more than 128 cluster subgraphs #2080
- `dot2gxl` no longer crashes on input `<node id="">` #2092
- remove itos #2229
- `sfdp` no longer crashes on certain graphs with cycles. #2225
- `gml2gv` does not handle integer penwidth correctly #1871

### Removed

- the glitz plugin has been removed. The plugin was never complete and
  distributions no longer ship glitz.

## [3.0.0] – 2022-02-26

### Changed

- **Breaking**: Using Graphviz as a library on Windows now requires the `GVDLL`
  symbol to be set to ensure correct linking.
- **Breaking**: Graphviz headers no longer define the `boolean` type. A
  replacement is C99 `bool` in the C standard library’s stdbool.h.
- **Breaking**: The `insidefn` member of the `shape_functions` struct must now
  be a pointer to a function returning a C99 `bool` instead of a
  Graphviz-specific `boolean`.
- **Breaking**: The `swapEnds` and `splineMerge` members of the `splineInfo`
  struct must now be pointers to functions returning a C99 `bool`s instead of
  Graphviz-specific `boolean`s. Similarly the `ignoreSwap` and `isOrtho` members
  of this struct must now be C99 `bool`s instead of a Graphviz-specific
  `boolean`s.
- **Breaking**: The `defined`, `constrained`, `clip`, and `dyna` fields of the
  `port` struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `set` and `html` fields of the `textlabel_t` struct are now
  C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `usershape` field of the `shape_desc` struct is now a C99
  `bool` instead of a Graphviz-specific `boolean`.
- **Breaking**: The `candidate` and `valid` fields of the `rank_t` struct are
  now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `filled`, `landscape`, and `centered` fields of the
  `layout_t` struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `has_images`, `has_flat_edges`, `has_sourcerank`,
  `has_sinkrank`, `expanded`, and `exact_ranksep` fields of the `Agraphinfo_t`
  struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `clustnode` and `has_port` fields of the `Agnodeinfo_t`
  struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `conc_opp_flag` field of the `Agedgeinfo_t` struct is now a
  C99 `bool` instead of a Graphviz-specific `boolean`.
- **Breaking**: The `must_inline` and `nocache` fields of the `usershape_t`
  struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `device_sets_dpi`, `external_context`, `fit_mode`,
  `needs_refresh`, `click`, `has_grown`, and `has_been_rendered` fields of the
  `GVJ_t` struct are now C99 `bool`s instead of Graphviz-specific `boolean`s.
- **Breaking**: The `loadimage` member of the `gvloadimage_engine_t` struct must
  now accept a C99 `bool` parameter instead of a former Graphviz-specific
  `boolean` parameter.
- **Breaking**: The `textlayout` member of the `gvtextlayout_engine_t` struct
  must now return a C99 `bool` instead of a Graphviz-specific `boolean`.
- **Breaking**: The `config` and `auto_outfile_names` members of the
  `GVC_common_s` struct are now C99 `bool`s instead of Graphviz-specific
  `boolean`s.
- **Breaking**: The `fixed` member of the `pack_info` struct is now an array of
  C99 `bool`s instead of an array of Graphviz-specific `boolean`s. Similarly,
  `pack_graph` now takes a `bool` array instead of a `boolean` array.
- **Breaking**: `pccomps` now takes a C99 `bool` instead of a `boolean` output
  parameter.
- **Breaking**: `gvusershape_file_access` now returns a C99 `bool` instead of a
  Graphviz-specific `boolean`.
- **Breaking**: 1-bit fields of the `obj_state_s` struct are now unsigned
  instead of signed.
- **Breaking**: Graphviz headers no longer define the constant `MAXSHORT`. A
  drop-in replacement is `SHRT_MAX` in the C standard library’s limits.h.
- **Breaking**: Graphviz headers no lnger define `NIL` macros. A drop-in
  replacement is `NULL` in the C standard library’s stddef.h.
- **Breaking**: Graphviz headers no longer define the `NOT` macro. A drop-in
  replacement is the C/C++ operator `!`.
- **Breaking**: Graphviz headers no longer (re-)define the C constants `INT_MIN`
  and `INT_MAX`. Replacements can be found in the C standard library’s limits.h.
- **Breaking**: Graphviz headers no longer define the constant `_DUMMY_ELEM`.
- **Breaking**: The `-m` memory test option to Graphviz command-line programs
  has been deprecated. Tools such as
  [Leak Sanitizer](https://clang.llvm.org/docs/LeakSanitizer.html) are a more
  effective modern way of diagnosing memory safety issues.
- **Breaking**: Graphviz headers no longer define the constant `MAXFLOAT`. A
  replacement is `FLT_MAX` in the C standard library’s float.h.
- The Ming plugin that produced Shockwave files has been removed. This format
  was EOLed by Adobe in April 2019. #2160
- CentOS 7 packages now include libmingle and the `mingle` program.
- The tclpkg Makefile no longer suppresses `-fstack-clash-protection` nor
  other compiler options containing `-x`
- Lefty is no longer enabled in the portable source tarball.
- on Linux, the CMake build system uses the standard `GNUInstallDirs` to locate
  target installation paths

### Fixed

- **Breaking**: GVPR now typedefs `ssize_t` as `SSIZE_T` on Windows instead of
  `int` #1804
- **Breaking**: `vgpanecmd` in the TCL tclpathplan library no longer accepts
  abbreviated commands (e.g. `r` for `rotate`) and commands must be given in
  full #1961
- fix detection of unavailable output format
- SVG layout doesn't always respect requested size #1855
- mismatched format string in `mingle`
- Building from scratch with Visual Studio fails #2175
- Plugins are not configured on MinGW #2176
- gvpr on MinGW does not support absolute paths #1780
- PNG format not available in CMake builds with MinGW
- tclpkg Makefile corrupts CFLAGS #2177
- lneato -? sometimes fails with STATUS_STACK_BUFFER_OVERRUN on Windows #1934
- expr misinterprets `<<` and `>>` #2103
- stdout and stderr are not flushed at exit on MinGW #2178
- Gvedit on macOS now understands the `-?` help command line argument
- CMAKE_LIBRARY_PATH is not honored #1973
- assert failure with `nslimit1=0` #1902
- `gvpr` usage output has been restored. This was a regression in Graphviz
  2.46.0.
- C++ API not usable after install #2196

## [2.50.0] – 2021-12-04

### Added

- hard-coded lookup tables for fallback font metrics for more fonts and font
  variants
- a new `gvputs_nonascii` API function has been implemented for GVC I/O with C
  escaping

### Changed

- Check for existence of `dl_iterate_phdr(3)` and if it is available, prefer
  using it instead of iterating over `/proc/self/maps` for computing `libdir`.
- A limit on GVC config files of 100000 bytes has been removed.
- MD5 checksums of release artifacts are no longer provided. SHA256 checksums
  are still provided and these should be used instead.
- when cross-compiling, the `dot -c` is no longer run during installation
- `$CMAKE_INCLUDE_PATH` is no longer manually configured in the CMake build
  system

### Fixed

- remove Bashism from `gvmap.sh` #2151
- Lefty artifacts are no longer installed when Lefty is disabled #2153
- Smyrna artifacts are no longer installed when Smyrna is disabled
- calling convention mismatches in delaunay.c’s GTS code
- impure assertion in `jacobi`
- undefined behavior in libgvc’s reading of little endian numbers
- boldness of `agnxtsubg` in cgraph man page
- parameter `name` in `gvusershape_find` prototype corrected to a const pointer,
  to match the implementation
- xdot JSON output is not valid JSON #1958
- fix uninitialized read of `pid` in `_sfpopen` on Windows
- claimed minimum CMake version supported has been corrected to 3.9

## [2.49.3] – 2021-10-22

### Fixed

- gvpr "split", "tokens", and "index" functions produce incorrect results #2138.
  This was a regression that occurred between 2.47.1 and 2.47.2.

## [2.49.2] – 2021-10-16

### Changed

- Lefty is disabled by default in the Autotools build system. To re-enable it,
  pass `--enable-lefty` to `./configure`. In a future release, Lefty will be
  removed.
- remove PHP5 support in SWIG bindings

### Fixed

- Msys experimental packages are included in release artifacts #2130
- CMake build system incorrectly aliases gv2gml to gml2gv #2131
- Gv2gml Doesn't escape quotes in attributes #1276
- GVPR incorrectly understands color schemes #1956

## [2.49.1] – 2021-09-22

### Changed

- the CMake build system installs gzipped man pages if `gzip` is available #1883
- CMake projects using Graphviz as a subproject (`add_subdirectory`) can now
  link against and use `gvc`.

### Fixed

- various problems in the generation of Javascript bindings
- 2.48.0: test suite is failing #2112
- Ensure correct file-level dependency for generated file in cmake generated
  projects #2119
- compile failures with a C++20-compatible toolchain #2122
- compile errors on macOS when using Bison 3.8 #2127
- Make Graphviz buildable as a cmake subproject/subdirectory #1477
- Header not found in Cmake project #2109

## [2.49.0] – 2021-08-28

### Added
- a very basic C++ API for a subset of the functions in lib/cgraph and
  lib/gvc, allowing a graph to be rendered from DOT source to a
  specified format. The new API is available through two new
  libraries: lib/cgraph++ and lib/gvc++. It is experimental, meaning
  that it might have breaking changes also in upcoming patch or minor
  releases (towards #2001)
- CMake builds now support an `with_expat` option that allows the support for
  using HTML-like labels through the optional expat library to be explicitly
  enabled (default) or disabled
- CMake builds now support an with_zlib option that allows the support for
  raster image compression through the optional zlib library to be explicitly
  enabled (default) or disabled

### Changed

- the CMake build system now enables `-Wextra` when building C++
- some Cgraph functions that take `char*` arguments that they do not modify have
  been updated to take `const char*` arguments #634
- incorrectly using the `layout` attribute on anything other than a graph now
  results in a warning about this being invalid #2078
- `edgepaint` accepts more standard `--` prefixed command line arguments and
  rejects invalid options #1971
- improved detection of Lefty dependencies in the Autotools build system
- libexpr rejects printing the time (`%t`) if no format is provided
- `-DDATE=…` option in the CMake build system has been removed
- the Autotools build system no longer writes the DATE file and the portable
  source tarball no longer includes this

### Fixed

- The attached dot file causes a segfault when processed #2095
- fix typos and update URLs in `edgepaint` usage text and man page
- Fix clang's undefined behavior warning in dotLayout
- gvpr doesn't build on macOS but MKDEFS_EXECUTABLE points to wrong
  directory #2101
- the generated gdefs.h header is no longer installed
- `ccomps` out-of-memory message no longer incorrectly refers to `gc`
- do not abort when `calloc(0, x)` or `calloc(x, 0)` in `gcalloc` return `NULL`
- failed Exshort_t type discrimination #1799
- dot manpage is in wrong directory on Windows #1936
- CMake builds fail when when the ltdl library is not available even if the
  `enable_ltdl` option is `ON`
- CMake builds fail when when the optional `zlib` library is not available
- fix graph rotation in quartz plugin

## [2.48.0] - 2021-07-17

### Added
- a new C++ test infrastructure based on ctest and Catch2 towards #2002
- support for test coverage analysis with
  [lcov](http://ltp.sourceforge.net/coverage/lcov.php) and
  [test coverage visualization in GitLab MRs](https://docs.gitlab.com/ee/user/project/merge_requests/test_coverage_visualization.html)

### Changed

- updated Graphviz bug report URL in the Autotools build system
- Fix `WIN32` path of `gvc.def` specified in `libgvc_la_LDFLAGS`
- the CMake build system now not only checks for Bison, but also ensures the
  found version is recent enough #1916

### Fixed

- ortho's eqEndSeg incorrectly modifies its arguments #2047
- Autotools enables -Wtrampolines and -Wlogical-op for Clang #2066
- node_distinct_coloring failure due to out-of-memory now reported correctly
  instead of referring to a failure to open lab_gamut
- Fix a typo `GD_LIBS` to `GDLIB_LIBS` in `tclpkg/tcldot/Makefile.am` !2022
- Autotools build system sets libgd variables now instead of incorrectly setting
  GTK variables
- HTML strings used as labels are distinguishable in GXL output by
  `kind="HTML-like string"`
- a Bashism removed from the Autotools build system
- when Criterion is available, the `command_line` test binary is no longer built
  and installed by default, but rather during `make check`
- round-tripping a file through ``gv2gxl`` and then ``gxl2gv`` no longer causes
  HTML-like labels to become non-HTML like labels #517
- fix ODR violation by including the ortho object files in the gvc lib also for
  CMake and MSbuild #2096

## [2.47.3] - 2021-06-19

### Changed

- marginally more accurate computations in Smyrna sphere projection
- Perl is no longer required to build Graphviz #2067
- nop more reliably returns success and failure exit statuses
- implicit 63 character limit on plugin names is removed in GVC
- the packaging work flow on CentOS 7 now selects the Python 3 bindings, instead
  of Python 2 bindings
- remove Python 2 bindings #1992
- improved thread-safety in Graphviz bindings Makefile

### Fixed

- no longer pass libcommon to the linker twice in mm2gv when building with CMake
- Quartz plugin is now compiled with explicit `--tag=CC` to libtool #2065
- out of bounds read when attempting to load a plugin whose name is ≥63
  characters
- out of bounds read when examining a registered plugin whose name is ≥63
  characters
- do not `fclose(NULL)` in gvmap
- Assertion error when using `dot` in `ortho.c` in `convertSPtoRoute` in
  graphviz 2.47.2 #2082. This was a regression introduced in 2.47.2.

## [2.47.2] - 2021-05-26

### Added

- CMake option `use_sanitizers` which enables building with address and
  undefined behavior sanitizer

### Changed

- $PATH is no longer assumed to be "/bin:/usr/bin:/usr/local/bin" if unset
- test suite no longer assumes `python3` exists #2049
- CMake build system locates Python 3 before calling it
- diff and grep are no longer required to build Graphviz on Windows

### Fixed

- Uninitialized variable read in delaunay_tri
- potentially mismatched format string in tclpkg
- `gvToolTred` is now exported from gvc.dll on Windows mirroring the behavior on
  Unix platforms.

## [2.47.1] - 2021-04-17

### Changed

- More detailed error messages when opening input file fails

### Fixed

- Windows build thinks xdg-open can be used to open a web browser #1954
- lab_gamut_data misses a value #1974
- xdot man page does not document some functions #1957
- Superfluous empty `@param` in documentation #1977
- PIC renderer does not work and probably never has #131
- dot conversion to dia format #689
- memory leak of reference-counted HTML strings
- Align rank from bottom in dot graph #1339
- Fix for TBbalance attribute code #1980
- HTML parser error with single closing square bracket in table row #1893
- reference counted strings put the HTML bit in the middle of the reference
  count #1984
- &amp;amp; escape disappearing #797
- miscalculation of minimum rank on large graphs
- AddressSanitizer: strcpy-param-overlap in gvconfig_libdir when
  running dot -c #1994
- fix reuse of va_list in pov rendering

## [2.47.0] - 2021-03-15

### Changed
- The edges in JSON output are ordered now !1728
- remove regex usage #1919
- RxSpencer is no longer a dependency on Windows
- gvmap.sh is compatible with POSIX shells in addition to ksh
- sed is no longer a build dependency on Windows
- SHA256 checksum generation? #1955

### Fixed
- Fix gvpr -? to actually print usage and exit non-zero
- gvpr is not built by CMake #1878
- typos in gpcanvas.c #1927
- memory leak in libmingle
- private inheritance in IncVPSC #1874
- broken sorting in nearest_neighbor_graph_ann.cpp #1938
- memory leak in ANN bridge
- gvpr on Windows does not support absolute paths #1780
- buffer overflow in unflatten
- agxbputc macro does not bracket its arguments #1814

## [2.46.1] - 2021-02-13

### Added
- Support for building against Guile 2.2
- Portable source is now also offered as a .tar.xz

### Changed
- CentOS/RHEL 6 is no longer supported
- Vestiges of Qt4 support have been removed
- C++11 support is now required of the C++ compiler used to build Graphviz
- C99 support is now required of the C compiler used to build Graphviz
- Question about userout() function in agerror.c #1924
- The minimum version of Python required to run the test suite is 3.6

### Fixed
- memory leak in label construction
- gvedit compilation errors out, but works if manually compiled with qt5 #1862
- incorrect HTML BR attribute parsing code #1913
- broken overflow checks in RectArea #1906
- various memory leaks !1699
- Fix bad free in lefty !1709
- typo in pathcross #1926
- Out-of-bounds write caused by incorrect error handling of malloc in genUserdata #1928
- Offer .tar.xz files too #454
- Header file graphviz_version.h has no include guards #1929
- regression: newlines embedded in quoted labels / node names are not preserved in 2.46.0 #1931
- Properly fill graphviz_version.h !1706

## [2.46.0] - 2021-01-18

### Added
- Cgraph's agxbuf API gained a new function agxbdisown(), for dissociating
  backing memory from the managed buffer
- Build system support for the Elbrus 2000 CPU, thanks to Michael Shigorin

### Changed
- Cgraph's agheap() API has been removed
- Autotools build system support for eFence has been removed
- Building Graphviz with ICC defaults to -O2 instead of -O0
- Build system work arounds for GCC 3 have been removed
- Incomplete support for running the test suite under CMake has been removed
- Portable source tarballs now use the “ustar” POSIX format
- Minimum version of Flex required to build Graphviz is now 2.5.2
- Minimum version of Bison required to build Graphviz is now 3.0
- Minimum version of CMake required to build Graphviz using CMake is now 3.1

### Fixed
- gvpr: line numbers in gvpr errors/warnings are incorrect #1594
- URL typo in patchwork man page
- Escaped backslashes are not correctly handled when producing xdot with dot #165
- heap-over-flow(off-by-null) in lib/common/shapes.c #1700
- Windows MSBuild executables have the wrong version #1745
- Cast Overflow at pango_textlayout #1314
- x11 back end segfaults if display is unavailable #1776
- typo in cmd/gvpr/lib/clustg #1781
- Segfault in dot #1783
- Incorrect 'Arrow type "s" unknown' error #1444
- segfault on reading 0x10 #1724
- Null-dereference READ (144736912) #1676
- "Warning! PATH too long installer unable to modify PATH!" using CMake Windows installer and PATH length > 1024 #1770
- gvedit -? gives "option - unrecognized - ignored" instead of showing usage #1813
- lefty is not built for Windows (fixed for MSBuild builds only) #1818
- a failure to detect OpenGL glGenTextures() errors has been corrected
- sfio does compile time benchmarknig #1422
- iffe "lib" check always succeeds when compiler optimises #1521
- syntax error near text who is not present #1411
- Explicitly links with libstdc++; should allow libc++ if appropriate #163
- A macOS file that was erroneously excluded from portable source tarballs has
  been restored
- Add option -? for usage to diffimg
- Add option -? for usage to dotty
- Add option -? for usage to lneato
- Add option -? for usage to vimdot
- Fix smyrna -? to actually print usage instead of error
- Fix edgepaint -? to actually print usage instead of error
- Remove '"' from usage text in non-Windows version of dotty
- Correct misspelled 'smyrna' in usage
- Fix edgepaint -o option
- Correct shebang of gvmap.sh to use ksh
- Fix gvmap.sh -? option to exit with zero exit status
- Graphviz doesn't build on MacOS with the latest libc++ #1785
- make fails if ps2pdf is not installed (using autotools) #1763
- multiple graphs to file output causes a segfault #1845
- lefty PTY functionality relies on file descriptor implementation details #1823
- buffer overflow in fdpgen
- Crashes by VRML output when current directory is not writable #793
- Segmentation fault when newrank=true #1221
- sfdp craches #236
- fdp segmentation fault with GK=0 #1290
- fdp crash #1865
- Graphviz always crash with this simple dot file #167
- Seg fault in dot #1771
- gml2gv doesn't handle some attributes correctly #1869
- Add missing circo, fdp, neato, osage, patchwork, sfdp & twopi tools to Windows builds (copies of dot)
- Add gv2gml tool to CMake (copy of gml2gv on Windows, symlink to gml2gv otherwise)
- Regression: fdp generates internal names in the output #1876
- Regression: fdp assertion error on cluster in edge #1877
- Regression in id / &lt;title&gt; in svg for twopi #1907

## [2.44.1] - 2020-06-29

### Added
- applied RH patches (from graphviz-2.42.2-8.fc32.src.rpm)
  - graphviz-2.42.2-coverity-scan-fixes.patch
  - graphviz-2.42.2-dotty-menu-fix.patch
  - graphviz-2.42.2-ocaml-allow-const-cast.patch
- some allocation failures that could previously allow memory corruption now exit
- lab_gamut.3.pdf is no longer included in release archives

### Changed
- Windows binaries available at https://www2.graphviz.org/Packages/ instead of
  https://ci.appveyor.com/project/ellson/graphviz-pl238
- Retarget Windows builds to Visual Studio 2019 and toolset v142

### Fixed
- Released Ubuntu packages does not contain language bindings for Python3 #1737
- Neato's hier mode is broken since v2.44.0 #1726
- Segmentation fault (core dumped) #1436

## [2.44.0] - 2020-04-08

### Added
- New SGD mode in neato (thanks [Jonathan Zheng](https://gitlab.com/jxz12/graphviz/-/tree/sgd))
- Add pkg-config files !1322
- tred: add feature to output removed edges to stderr upon request !1326
- Fix issue #1671: Workaround: avoid creating a virtual edge loop. !1328
- Add riscv64 to host_cpu configure.ac !1329
- lib/cgraph: include empty malloc.h from subdir include !1332
- lib/gvpr: compile mkdefs with $(HOSTCC) rather than $(CC) !1333
- lib/vpsc: rename bcopy->b_copy !1334

### Fixed
- MSB4018 The NativeCodeAnalysis task failed unexpectedly. #1481

## [2.42.4] - 2020-04-05

### Added
- Include all test files in distro !1341
- host_cpu add mips64 platform !1325
- Correct description of 'port' syntax in manual !1324

### Fixed
- svg output displays TITLE of %3 if graph had no name #1376
- XML errors in generated SVG when URL attribute contains ampersand (&) #1687
- Test files missing from source distributions #1647
- SVG error for "g.transform.scale " in graphviz version 2.43 #1605

## [2.42.3] and earlier

```

October 9, 2019
    - Release 2.42.3
	- Merge 1316, 1317, 1319, 1320
	- Patches #1591, #1596
	- Add Fedora 32 builds
	- CI/CD fixes
	- Documentation (Warning about HTML label usage)

September 17, 2019
    - Release 2.42.2 - ( Never fully released due to CI/CD hardware issues )
    - Fix deployment issues.  Builds can now be found under:
	             http://www2.graphviz.org/Packages/
July 17, 2019
    - Release 2.42.1
    - Fix deployment issues.  Builds can now be found under:
	             http://www2.graphviz.org/Packages/
July 4, 2019
    - Release 2.42.0
    - Fixes quite a few bugs
September 22, 2017
    - Move master repo to GitLab: https://gitlab.com/graphviz/graphviz
December 21, 2016
	- Remove usage of ast_common.h
December 20, 2016
    - Release 2.40.0
        - network-simplex fixes and optimization (Stephen North)
	- built-in tred tool now available in the various swig generated
	language bindings (John Ellson)
	- number rounding added to SVG renderer (same as PS and TK rounding)
	to aid regression testing. (John Ellson)
	- additional regressson test framework, used in Travis CI builds. (Erwin Janssen)
	- PHP7 support (requires swig-3.0.11 or later). (John Ellson)
	- Allow user to specify clustering algorithm in gvmap. (Emden Gansner)
	- Add Sierpinski graph generator to gvgen. (Emden Gansner)
	- Extensive code cleanup (Erwin Janssen)
	- Removal of libgd source - use vanilla libgd from separate install
	- Windows builds (Erwin Janssen)
	- Appveyor CI for automated Windows build testing (Erwin Janssen)
	- Travis CI for Fedora/Centos builds (Erwin Janssen)
	- Added JSON output format, -Tjson  (Emden Gansner)
	- New curved arrowhead, cylinder node shape.
	- Resolves bugs: 2599, 1172
June 18, 2016
	- Experimenting with Travis CI
February 13, 2016
	- Add cylinder shape for databases.
	- Free installed plugins
	- Update makefile for dot so that the using libpanco_C in the static build include PANGOFT2
        as well as PANGOCAIRO_LIBS (needed for some versions of Ubuntu)
February 1, 2016
	- Add json output format
April 26, 2015
	- output class value in svg files
September 9, 2014
	- Add plain shape for use with HTML-like labels.
August 12, 2014
	- Add icurve arrowhead.
July 28, 2014
	- Revert to old, translate to origin semantics in neato, etc. Add flag notranslate if that is
          what the user desires.
April 13, 2014
	- Release 2.38.0
	- Resolves bugs: 2409, 2413, 2417, 2420, 2422, 2423, 2425
March 27, 2014
	- Enable packing for dot
	- Allow scaling to work for all non-dot layouts
March 9, 2014
	- Add overline text characteristic.
March 4, 2014
	- Fix bugs in gvpr and gv.cpp so edges can be created in subgraphs.
	- Add edgepaint program for coloring edges to make them easier to tell apart.
	- Modify neato to avoid unnecessary translations of output. This allows positions
	given on input to remain the same on output.
	- Fix swig java package to work and support gv.renderresult.
	- Fix test for the absence of layout (old test relied on statically allocated Agraphinfo_t).
	- HTML-like tables and cells can now specify which borders should be drawn.
	- The fixedsize attribute now takes the value "shape" which allows labels much larger than the
	node shape.
January 11, 2014
	- Release 2.36.0
	- Resolves bugs: 2372, 2384, 2388, 2391, 2392, 2383, 2395, 2401, 2406
	- Various MacOS Fixes from Pixleglow.
	- Remove old libgraph sources from distributions.
	- Move master git repo to github.com
September 15, 2013
	- Add <S> element for strike-through to HTML-like labels.
September 6, 2013
	- Release 2.34.0
	- New version of xdot supporting inline text characteristics such as <b> and
	version-specific output based on xdotversion
	- Resolves bugs: 2325, 2326, 2333, 2334, 2337, 2338, 2340, 2343,
		2345, 2346, 2349, 2350, 2351, 2352, 2353, 2354, 2357, 2359
	- Resolves Redhat bug: BZ#847458
August 21, 2013
	- Added mingle command and library for edge bundling
August 1, 2013
	- Release 2.32.0
	- New version of xdot format, annotating gradient color schemes
	- Support for reading pdf images using poppler
	- Lefty/dotty/lneato now accept anonymous graphs
July 2, 2013
	- Add star node shape
	- Add two-tone (non-gradient) fill
February 14, 2013
	- Release 2.30.1
	- various build fixes
January 13, 2013
	- Release 2.30.0
	- Replaced libgraph with libcgraph; use of libgraph is now deprecated
	- New ranking algorithm that allows multiple subgraph constraints
November 27, 2012
	- Add graphml2gv and gv2gml to Windows package.
September 25, 2012
	- Support edges using curved arcs.
August 16, 2012
	- Added new shapes used by the synthetic biology community.
July 12, 2012
	- For HTML-like labels, provide rounded cells, and dashed or dotted borders.
	- Add lcurve and rcurve arrowheads.
	- Add prototype pie chart and striped fills.
	- Support insets in treemaps to make containment clear
June 7, 2012
	- Add random rooted tree generation to gvgen
February 29, 2012
	- Allow GVPRPATH to specify that the default path be prepended or appended to it.
February 27, 2012
	- Support arbitrary lists of layers; allow the user to specify arbitrary layers for output.
February 24, 2012
	- A collection of gvpr scripts, which were part of the source package, are now
	installed in <prefix>/share/graphviz/gvpr, and the that path is used as part of th
	default built-in path for gvpr.
February 15, 2012
	- Update libexpr to reflect Glenn Fowler's changes including scopes for variables.
February 9, 2012
	- Add next graph variable to gvpr
February 8, 2012
	- Modify dot and fdp so that a cluster's margin attribute will affect the space
	  between the bounding box and nodes
January 26, 2012
	- Modify the dijkstra tool to use only directed edges
	- Output numbers without quotes if quotes are not needed on input
	- Support gradient fill
January 23, 2012
	- Provide support for webp images
January 17, 2012
	- Fix tapered edges to use the dir attribute and arrowhead
September 21, 2011
	- Add imagepath attribute
	- Add help functionality to Graphviz.app
August 24, 2011
	- Add <B>,<I>,<U> to html strings via svg
August 16, 2011
	- Add tapered edges
August 3, 2011
	- Add support for external labels
July 14, 2011
	- Add initial implementation of graphml2gv
July 8, 2011
	- Add basic horizontal and vertical rules to html tables
May 6, 2011
	- Release 2.28.0
	- incremented library api version:
	    libcdt, libgraph, libcgraph, libgvpr, libgvc
	- Add gvmap, cluster and gvmap.sh
	- Deprecate dotty; prefer gvedit
	- Add patchwork supporting squarified tree maps
	- Add ordering as a node attribute
	- Fix problems with font resolution
	- Fix problems with text placement
	- Fix twopi to set root attribute
	- Make available layouts and formats available via the API
	- Fix error message system so that an application can capture the messages
	- New Qt-based version of gvedit
	- New attributes and features for sfdp
	- gvgen now allows the user to specify graph name and node name template
	- Make overlap=false denote overlap=prism
	- More efficient xdot library
	- HTML-like labels provide ID
	- Fixed bugs: 1480 1980 2044 2087 2088 2089 2091 2093 2094
	2095 2101 2102 2103 2104 2112 2113 2118 2128 2129 2139 2149
	2157 2113 2159 2160 2161 2163
March 31, 2011
	- Add many new gvpr scripts to release package
	- Add scale attribute to twopi
October 14, 2010
	- Add <B>,<I>,<U> to html strings via cairo
February 15, 2010
	- migrated to 2005 version of cdt
January 26, 2010
	- Release 2.26.3
	- libcgraph.so   version bumped from 4 to 5 due to API changes
	- Allow ranksep to specify multiple radial differences in twopi
	- Allow the user to specify all pairwise distances in neato with
	- Fixed bugs: 1280 1409 1567 1583 1624 1631 1655 1708 1709
	1727 1784 1792 1798 1800 1813 1814 1830 1831 1833 1836 1839
model=mds
December 10, 2009
	- Release 2.26.0
	- Core features:
		- added: "smyrna" - a new opengl-based viewer for large graphs
		- added: rudimentary "gml2gv", "gv2gml" converters
		- extended support for various image formats in node images
		- removed vestiges of codegens, now all putput formats supported
		  through plugins.  Dropped some output formats for which
		  plugins have not been developed: -Tdia, -Tmif
		- gvpr converted to a library; additional array handling and
		  text processing functions added; language extended to allow
		  multiple BEG_G/N/E blocks.
		- allow background images specified via xdot
	- Ports added/dropped from nightly builds:
	  (The dropped ports could probably be re-added if there was demand.)
		- added MacOSX SnowLeopard  (multiarch:  i386/x86_64/ppc)
		- added Fedora 12 (i386, x86_64)
		- added Fedora 13 (Rawhide) (i386, x86_64)
		- dropped Fedora 7 (i386, x86_64)
		- dropped Fedora 8 (i386, x86_64)
		- dropped RHEL 3 (i386, x86_64, ia64)
		- dropped Ubuntu 8 (i386)
	- Fixed bugs: 1683 1713 1718 1720 1738 1747 1759 1770 1776 1786
	  1799 1816 1827

June 16, 2009
	- Release 2.24.0
	- Core:
		- add new layout engine for large graphs: sfdp
		- add new layout engine for nested graphs: osage
        - pack library extended to handle array-based packing modes
          using array bounds, aspect ratio, user-controlled sorting, etc.
	- Fixed bugs: 1515 1590 1598 1601 1605 1607 1609 1610 1611 1614
	1615 1617 1625 1628 1634 1635 1640 1641 1642 1646 1649 1651 1652

March 13, 2009
	- Release 2.22.2
		- fix for buffer overflow (present in 2.22.0 and 2.22.1)
	- Fixed bugs:
		1602

March 9, 2009
	- Release 2.22.1
		- build fixes for Visual Studio and for FreeBSD
March 3, 2009
	- Release 2.22.0
	- Core:
		- libgvc api changed, version bumped.  Affects third party
		  applications using libgvc.
		- plugin api changed, version bumped.  Affects third party
		  plugins for graphviz.
		- 90% conversion to cgraph has been done, but not enabled yet,
		  (and yes, its true what they say about the last 10% )
		- drop libagraph from distribution  (use libcgraph)
		- layout code completely converted to floating point.
		- new "dot -P" option for generating a graph of available
		  plugins.
		- registered MIME type:  text/vnd.graphviz for .gv files
		- rename files from .dot to .gv to avoid conflict with
		  Word templates.  .dot still supported, but deprecated.
		- new command: mm2gv   (matrix-market graph file conversion)
		- rename commands:	dot2gxl -> gv2gxl
					gxl2dot -> gxl2gv
	- Plugins:
		- new rsvg plugin for support of node shapes in SVG format
		- new gs plugin for support of node shapes in PS format
		- new lasi plugin for support of UTF-8 characters in PS output
		  (the above thee plugins are Linux only, at the moment)
		- new quartz plugin (MacOSx only)
		- new gdiplus plugin (Windows only)
		- new -Tvml support in core plugin (thanks Steve Roush)
		- new -Ttk support in core plugin (also used by Tcldot and
		  gv_tcl language bindings.)
		- disabled old style codegens completely
	- Linux:
		- new Ubuntu8 builds
		- new Fedora 10 and 11 builds
	- MacOSx:
		- Universal binary for Leopard: i386, x86_64, ppc, ppc64
		- Should not conflict with parallel install of MacPorts
		  version of graphviz
		- Improved GUI
	- Windows:
		- VisualC project files now available, in addition to the GNU
		  Makefiles that are used the mingw builds.
	- Language Bindings:
		- fixed problem with writing dot, xdot, plain, canon to
		  memory or to Tcl_Channels
		- renamed man pages to Debian style:  gv.3tcl, gv.3perl, etc
	- Fixed bugs: 827 1365 1366 1367 1368 1374 1375 1376 1378 1380 1382
	1383 1385 1386 1388 1390 1391 1392 1394 1395 1397 1398 1399 1405
	1407 1410 1412 1414 1415 1416 1421 1424 1425 1427 1429 1431 1433
	1435 1436 1437 1438 1440 1441 1444 1446 1451 1452 1453 1456 1457
	1459 1460 1461 1462 1463 1464 1465 1466 1470 1474 1475 1476 1477
	1478 1484 1485 1489 1490 1492 1493 1495 1496 1499 1500 1501 1502
	1503 1505 1509 1513 1521 1523 1525 1530 1531 1532 1533 1535 1536
	1539 1540 1542 1543 1546 1547 1551 1553 1554 1561 1565 1566 1568
	1569 1570 1571 1573 1577 1578 1579 1580 1581 1582 1584 1586

June 25, 2008
	- Release 2.20.2
	- Fix bug in HTML-like labels
June 23, 2008
	- Release 2.20.1
	- Fix bug in ccomps related to conversion to cgraph
June 20, 2008
	- Release 2.20.0
	- Preparing for Release 2.20
	- Fixed bugs: 1315, 1317, 1324, 1336, 1343, 1364
	- Add new "folder" shape for nodes.
	- Migration of gvpr tools to libcgraph.
	- New output format -Teps  (encapsulated postscript)
	- Various NetBSD and SuSE fixes incorporated
	- ./configure now provides a summary
	- RPM specfile updates for fedora-10 (no more string comparisons)
	- Add MacOS support (Glen Low)
March 10, 2008
	- Release 2.18
	- Fixed bugs: 1249 1255 1256 1268 1276 1289 1295 1300
		Fedora BZ#247376,
	- in -Tps use a new number formatter that suppresses trailing 0.
	- support tcl/tk-8.5
	- support gcc-4.3
	- support for node usershapes/images in svg format (thanks Alex Poylisher)
	- install: perl, php, python, ruby, tcl, bindings in language-specified directories
	- add arrowhead scaling with edge penwidth
	- add "folder" node shape (thanks Pander)
	- many windows and mac fixes (thanks Glen)
	- add "smyna" large graph view (thanks Arif) (not yet included in binary distros)
December 12, 2007
	- Release 2.16.1
	- Fixed bugs: 1228 1234 1238 1239 1245
	- Improvements to PHP binding
	- Improvements to OCAML binding
	- Make regression tests run from the build tree, rather than require installation
	- Repair freetype detection on RedHat-7 (Yes, people still use it!!)
	- Fix zoom-at-mouse-location in -Txlib and -Tgtk
	- Fix some dotty regressions
November 9, 2007
	- Release 2.16
	- Fixed bugs: 456 473 1021 1153 1154 1155 1159 1160 1162 1165 1166
	1168 1169 1170 1172 1173 1174 1175 1177 1178 1179 1181 1182 1183
	1185 1187 1189 1192 1193 1195 1196 1199 1204 1207 1210 1215 1216
	1217 1218 1219 1220 1223
	- new regression test suite
	- new cgraph library (will eventually replace graph and agraph)
	- add "image" and "imagescale" for simpler support for images in nodes
	- add "tab" "box3d" and "component" shapes.  - Diomidis Spinellis
	- replace arith.h in distro
	- add functions to access version info to avoid need for gvcint.h
	- Fix problem with irregular character spacing at 96dpi in pango/cairo output formats.
	- Add gdk_pixbuf plugin providing: .bmp .ico .jpg .png .tif
	- Add DevIL plugin providing: .bmp .jpg .png .tif .tga
	- Extend GD plugin to provide a backend to cairo for: .gif .jpg .png .gd .gd2 .wbmp  <- gifs are now antialiased
	- Rework plugin framework to separate device from renderer, and to autoload load dependendent plugins
	- show defaults in output from: ./configure --help
	- add more info to dot -v  and dot -v2 debug outputs
	- various issues with CR/LF in windows, but not in binary outputs.
August 15, 2007
	- release 2.14.1
	- Fixed bugs: 1163, 1167
	- Windows build fixes
	- Add xdot parsing library to source distros
	- graphviz.spec fixes for rpm distros from Gareth Armstrong
	- moved language binding man pages to mann (gv_php.n, gv_ocaml.n, etc.)
	- New access functions for version info in GVC_t - permits gvcint.h to
	be private.
August 2, 2007
	- release 2.14
	- Fixed (or otherwise closed) bugs:
		74 130 162 184 190 197 219 223 281 295 311 316
		324 352 364 385 393 404 420 447 455 474 489 507
		530 532 537 543 551 564 571 574 577 583 587 588
		590 592 595 599 638 647 650 660 675 667 668 669
		676 684 685 686 721 725 734 740 746 747 748 749
		752 755 756 765 778 780 781 782 785 794 803 814
		822 828 836 840 847 852 862 866 868 893 928 944
		948 950 955 961 976 985 992 1024 1057 1064 1065
		1066 1069 1072 1074 1079 1085 1086 1089 1091 1092
		1093 1094 1096 1107 1111 1123 1124 1130 1138 1145
		1151 1152 1156
	- Fixed Redhat bugs: 218191 237497
	- Fixed Debian bugs: 321128 422862 422873
	- Fixed Gentoo bugs: 173676
	- Using system version of libgd if gd-2.0.34 or later. (Fedora 7 and 8 distros)
	        internal copy of gd updated to gd-2.0.35.
	- Updated GVGUI viewer for Windows
	- Windows build process now uses GNU autoconf and UWIN
	- Added support for selection of edge routing types:
		line, polyline, orthogonal, spline
	- Added -Tvml support
December 5, 2006
	- release 2.12
	- Bug fix release for 2.10
	- The gd plugin for font handlers was not being used at all if the build
	did not use fontconfig, e.g., on Windows. In addition, the code had
	dropped the name mapping to Windows font names.
	- PostScript output had an extraneous '%' character on the first line,
	which would cause printing to fail.
	- Text handling, during both sizing and layout, incorrectly handled
	empty lines such as label="\nabc".
	- HTML-like tables had been changed to use too much vertical space,
	to possibly use the wrong font in calculating the height of a line,
	and to use the wrong offset when moving the baseline from one line to
	the next.
November 27, 2006
	- release 2.10
	- dot - New pango+cairo renderer plugin (was in separate graphviz-cairo tree).
	  -- -Tpng now uses cairo   (-Tpng:gd for old gd based renderer)
	  -- -Tpdf now available
	  -- -Tps:cairo now available (-Tps is a direct ps renderer not based on cairo)
	  -- -Tsvg:cairo now available (-Tsvg is a direct svg renderer not based on cairo)
	  -- -Txlib now available -- "dot -Tx11 foo.dot"  watches foo.dot with inotify and updates
	  -- -Tgtk now available -- eventually to provide a graph editing capability - not fully working
	  -- -Tswf "Flash" now available using the ming library. Currently has incomplete font support and not yet in Fedora rpms because ming not yet available as rpm.
	- remove hard gd dependencies from dot.  gd renderers now provided
	  as optional plugin.   Deprecated, but required for -Tjpg, -Tgif and -Tvrml.
	- gvpr - Add kindOf function, plus functions to set and get default values
	- dot - Implement esep attribute to allow graph to specify room
	around nodes for spline routing.
	- neato - add vpsc library and DIGCOLA
	- neato - add IPSEPCOLA additions from Tim Dwyer
	- move: -Tps, -Tfig, -Tsvg, -Timap/ismap/cmap/cmapx, -Tdot/xdot,
	from codegens to a "core" plugin.
	- dot - new usershape plugin mechanism potentially supporting
	  a wider range of input shape format -> output format combinations.
	display on changes
	- Fixes for builds on Mac OS/X
	- dot - new -O switch to automatically generate output file
	names based on the input filename and the -T value.
	 e.g.  "dot -Tpng -O *.dot"
	Also works for case of multiple graphs in a single input file.
	- add support for "Brewer" color nameset
	- move reusable .so libraries to $PREFIX/lib per frequent request
	from Debian community.   Plugin .so's remain in $PREFIX/lib/graphviz.
	- Fix bugs 882 884 886 896 902 905 906 911 918 919 933 936 938 940
	   948 955 958 967 979 987 993 1005 1006 1011 1012 1013 1014 1016
	   1018 1025 1030 1034 1035 1039 1040 debian#37300

February 3, 2006
	- release 2.8
	- (POTENTIAL INCOMPATIBILITY) The default input scaling, in the
	absence of a "-s" switch, has been changed from inches to points.
	The new behavior of "neato" is equivalent to "neato -s72".
	The old behavior can be restored with "neato -s1".
	The purpose of this change is to avoid a Frequently-Made-Mistake
	when using "neato -n" to process a previously generated layout.
	Previously it was necessary to use "neato -n -s72", but with this
	change the default matches dot's output and the "-s72" is not required.
	- Added pseudo layout engines: "dot -Knop" and dot -Knop1" equivalent
	to "neato -n"
	- Added pseodo layout engine: "dot -Knop2" equivalent to "neato -n2"
	- Add support for color namespaces; add Brewer color data
	- Add support for simulated duplex edges using parallel edges:
	head arrow takes first color, tail arrow takes second color.
	- source code management moved back to CVS until GIT matures a bit more
	- distribute separe rpms for binares of language bindings :
	- Add a small pad region around graph renderings to allow for finite
	penwidths at the drawing edges
	- Add protonode(g) and E=protoedge(g) functions to simplify
	language bindings.
	- Add special purpose code to deal with html labels from language
	bindings.
	- Various portability fixes for: HPUX, Mac OS/X, Cygwin, Windows.
	- Fix bugs 784 786 787 788 789 790 791 793 795 796 798 799
	    800 801 804 806 811 812 817 820 821 823 824 825 830
	    837 839 841 842 843 848 850 851 854 855 856 857 858
	    859 861 863 866 867 869 872 874 876 877

August 28, 2005
	- release 2.6
	- experimentally moved source code management from CVS to GIT
	- added iterator functions to script bindings
	- more C-API tuning
	- add "-c" switch to dot to explicitly generate plugin "config" file
		instead of generating it as a side-effect of "dot -V"
	- better support for binary relocation.
	- plugin versioning and version checking
	- clean up of header files
	- provide statically linked "dot_static" (not incl. in rpms)
	- additional "event" support for GUIs (e.g. "DotEdit" graphviz-cairo)
	- add some information about plugins to "dot -v" output.
	- lefty/dotty fixes
	- fix bugs 746 750 752 753 754 756 761 763 764 765 768
		771 772 773 774 775 776 777 778
	- not a bug 757 760 770
July 20, 2005
	- release 2.4
	- major code restructuring
	- new plugin architecture (e.g. see separate package: graphviz-cairo )
	- new script-language bindings using swig (perl, ruby, python, tcl, java ... )
	- C-API now in libgvc (no more dotneato.[ch] or dotneato-config.sh]
	- pkgconfig now used for reusable libraries
	- lefty upgrade
	- fix bugs 156 255 492 631 641 647 659 662 665 670 690 691
			701 702 703 705 730 731 732 741 743
April 7, 2005
	- release 2.2.1
	- correct license headers to CPL in .cpp files
	- undo indentation cleanup to dynagraph .h files
	- fix bugs: 183 247 419 615 616 625 626 627 643
		646 651 658 661 664 674
	- fix buffer overrun in Gvfilepath construction
January 19, 2005
	- release 2.2
	- fix bugs: 86 345 517 579 580 597 600 601 604
	- use the original cpl1.0.txt as the license master, instead of CPL.html        - fix for bug generating in memory bitmaps that was affecting webdot
	- fixes for windows builds
	- documentation updates
December 11, 2004
	- release 2.0
	- new CPL license
	- re indent all sources
December 11, 2004
	- release 1.18
	dotneato
	- fix bugs: 451 536 545 547 548 559 561 565 572
	- increase max size  of HTML tables.
	- spline cluster edges in fdp
	- center userimages in nodes
	- support user images in HTML table cells
	- syntax extension for node:port:compass as well as node:compass
	- FreeBSD fixes
	- sync with gd-2.0.32
	- attempt to catch some out-of-memory conditions with very large graphs
	- support background and node-fill partial transparency when truecolor=true

September 14, 2004
	- release 1.16
	dotneato
	- fix bugs: 275 523 526 527 529 534
August 30, 2004
	- release 1.14
    dotneato
	- the official gd now has support support for GIFs again - the
		internal gd is now closely sync'ed with the official version
		and will eventually be removed in favor of using a
		separate installation of the official version.
	- gd has new support for FontConfig (thanks to Dag Lem)
		NB. the fontname attribute in graphs is now a font pattern
		as understood by fontconfig (e.g. fontname="Times-Italic"),
		unless it contains a '/' in which case it is interpreted as
		a font path as before.
	- gd provides support for html4 entities in decimal, hex or named, e.g "&lt;"
	- "dot -v" debugging output now reports fontname -> fontpath resolutions

	- PostScript generated by -Tps now uses "xshow" operator for strings
		for better matching of bitmap and PostScript outputs.

	- ability to use an external gd-2.0.29 version of libgd (EXPERIMENTAL)

	- new feature: parallel edges by using a ":" separated list of edge colors
	- new feature: rankdir=BT and rankdir=RL  (thanks to Dag Lem)

	- new layout engine: fdp - force directed placement (EXPERIMENTAL)
		a neato-like undirected layout engine that produces
		clustered symmetric layouts.
		Supports edges between clusters and nodes.

	- updated neato engine: now using stress majorization as the default,
		which avoids the potential for cycling
	- model=subset in neato provides a third distance function, where
		two nodes sharing many nodes will be place farther apart
	- shape=none now equivalent to shape=plaintext
	- fix label justification with \l and \r
	- first cut at <FONT> support added to html labels
	- various color transparency fixes
	- various fixes for UTF8 and Latin[12] character encodings.
	- various cluster fixes.
	- improved hyperlink support in -Tsvg
	- support tooltips on clusters in client-side imagemaps

    gvpr
	- add support for scanf and friends

    general
	- greater use of shared libraries.
	- pkg-config files provided for shared libraries (EXPERIMENTAL)
	- "./configure --disable-shared --enable-static" works if needed
	- C++ wrappers on all header files (thanks to Victor Wodecki)
	- various configuration and portablity fixes
	- provide pdf version of man pages
	- Windows package provides graphviz libraries and header files
	- Closed bugs: 195 198 234 321 330 399 401 406 410 411
		412 413 415 416 417 423 424 427 430 431 433 434 435
		438 441 442 444 445 449 450 452 454 457 458 462 463
		464 467 468 469 471 475 480 482 485 495 496 498 499
		500 501 504 508 511 512 514

March 5, 2004
    - added glyphwidths.ps support utility

March 1, 2004
    - release 1.12
    - general
	- rename bcc -> bcomps to avoid name conflict with "Bruce's C Compiler"
		on Redhat distributions.
	- all build without X11 (fix problem in lefty tree)
	- remove from distribution:
		dag, fdp, geo, grid, incr, shape, tcldgr, tcldgl
    - dotneato
	- fix "brown-bag" problem resulting in PNG and JPEG errors on RH8 and RH9.
February 23, 2004
    - release 1.11
    - general
	- fix windows builds
	- add tool "bcc" to distribution
    - dotneato
	- add -Gviewport="X,Y,Z,x,y"  where XY are the dimensions of a viewport
	  in device coordinates (pixels), Z is a zooming factor, x,y is the
	  location of the center of the viewport in graph coordinates.
	  Supported in bitmap and imagemap outputs only.
	- fix memory leak in gd/gdft.c
	- clean up calculation of whitespace around labels
    - dotty, lefty
	- fix for bug #400
December 23, 2003
	- added dijkstra (single source distance) filter
September 10, 2003
    - general
	- removed CVS directories from .tar.gz distributions
	- add "config" directory to contain some of the autoconf clutter
	- only remove flex products with "make maintainer-clean" to
	  avoid trying to regenerate them after "make distclean"
	  basically this is to avoid the broken flex on Debian.
	- suppress complaints from ./configure about config.rpath
	- doc/build.html updated with notes about Windows builds
	- build fixes for Forte 6sp2 compiler on Sun -xarch=v9a (64bit)
	- build fixes for OpenBSD
	- improved configure testing for Tcl/Tk
	- various bug fixes, internal restructuring, etc
    - dotneato
	- fix problem with extra escape chars in .fig output
	- support for "setlinewidth" in -Tfig
	- improved splines in -Tfig
	- add manpage for dotneato-config
	- neato: add defaultdist graph attribute to set distance
	  between components
	- first cut at html table formatter add. not ready for use yet
	  as the syntax is going to change some more.
    - tools
	- renamed "colorize" to "gvcolor" to avoid conflict on Debian
	- renamed "gpr" to "gvpr" to avoid conflict on Debian
	- add fflush() to acyclic, ccomps, gvcolor, tred, dot2gxl
	  to try to fix truncated output when used in php or perl cgi scripts
July 9, 2003
	- rerelease 1.10 with ast_common.h fix in -devel rpms
July 3, 2003
	- declare this version 1.10
	- general
	    - "mkdir obj;cd obj;../configure;make"   now works (bug #293)
	    - "make prefix=xxx"   now works (bug #274)
	    - "--with-wish=xxx"   now works (bug #270)
	    - remove generated file: ast_common.h from source distributions
	    - make GIF support configurable
	    - added .cvsignore throughout source tree to reduce CVS noise
	    - FAQ updates
	    - documentation updates for gpr
	    - improve portability of dotneato-config, but requires libtool now
	    - improvements to error processing for library users
	-gd
	    - sync with gd-2.0.15
	    - optimize line drawing code
	- dot, neato, twopi
	    - fix bugs 240 270 274 293 298 303
	    - support "peripheries=0" without crashing
	    - add support for "dia" output format (-Tdia)
	    - espf fixes (use of showpage)
	    - svg fixes (coordinates and viewBox)
	    - ismap/imap, fixes (quoting of label strings)
	    - fix to "point" shape
	    - improve (m|c|re)alloc usage
	    - improve handling of very-small fonts in bitmap outputs.
	    - various fixes for multiple -T -o feature
	    - add support for splines to records and ports (neato)
	    - various improvements to libpack
	    - dot_init_graph and neato_init_graph external for library users
	    - cluster improvements (neato)
	    - fix support for truecolor
	    - normalize splines so that they now always go from tail to head
	    - add some simple help text for any unrecognized option
		(e.g. -?  -h  --help)
	- tools
	    - extend gpr language to allow access to command-line arguments
	    - add sqrt() function to gpr
	    - add new tool - gvpack
	- tcldot
	    - use .dll extension if on windows
	    - doted demo
		- use tcl's file requestor instead of homebrew
		- add zooming controlled by mousewheel
		- support additional export formats

January 31, 2003
	- declare this version 1.9
		(3-level version numbering has been dropped now
		that we have nightly snapshot builds with their
		own extended numbering.)
	- general
	    - config.h is no longer installed.  config.h is generated by
		./configure for the current build only.  It may not be
		applicable for derivative builds.
	    - improve ICONV configure tests
	    - lots of janitor-work to clean up warning messages from -Wall
	    - use @OBJEXT@ in Makefile.am so that .obj is used under cygwin
	    - fixes for Solaris builds
	    - use libpng-config if available
	    - reduce long build times due to touching ast_common.h too often
	    - improve dependency tracking.  "make -j8" now works with distcc
	    - autogen.sh fixes to work on RH7.3, RH8.0, and Solaris.
	    - eliminate use of suffix rules which confused some makes.
	    - DOT language allows '+' for concatenation of quoted strings
	- dot, neato, twopi
	    - fix bugs 209 210 214 216 217 222 224 225 229
			230 233 236 237
	    - update gd into alignment with gd-2.0.9
	    - change to make libagraph output compatible with libgraph input
	    - add shapes: septagon, pentagon, a_ediamond, rect, rectangle
	    - introduce "ND_...", "ED_...", "GD_...", node/edge/graph-data
		accessor macros in partial preparation for use of
		libagraph in dot.
	    - add libdotneato.so, dotneato.h, dotneato-config
		to aid use of dot libraries by user apps based
	        on installed graphviz-devel rpm and without access
		to graphviz sources.
	    - new xdot output format providing detailed drawing instructions
	    - new -y command line flag, inverts y coordinates
	    - support multiple -T when -o given, as in:
			cat xxx.dot | dot -Tpng -Tcmap -o xxx
		which produces xxx.png and xxx.cmap from a single
		layout computation.   Intended for use in CGI programs.
	- agraph
	    - correct callback ordering for deletions
	- tools
	    - add gxl2dot and dot2gxl for GXL language conversions
	    - gvui now provides *map output
	- tcldot, tcldgr, tcldgl
	    - improve tcl8.4 support
	    - extend search path for tcl.h to include /usr/local/include/tcl8.4/
		in support of BSD install conventions.
	- dynagraph
	    - many fixes
	    - change to not build dynagraph by default (use --with-dynagraph)
	- docs
	    - dotguide updates
September 27, 2002
		- declare this version 1.8.10
	- general
	    - various configure.in fixes and simplifications
	    - change configure to now build dynagraph by default
	    	"--without-dynagraph" is supported
	    - fix graphviz.spec.in to partition packages properly
	    	graphviz no longer depends on graphviz-tcl.
	    -  Makefile.old cleanups
	    - configure.old now set version number automatically from
	      configure.in
	- dot, neato, twopi
	    - Initial support for image node shapes + URL fetch.
	    - Made number of dimensions a runtime variable in neato.
	    - Bug fix in vrmlgen for degenerate splines.
	    - Bug fix - ordering=in should now work
	    - Bug fix - layers no numbered from 0 to match PS requirements
	    - Bug fix - don't draw arrows on invisible edges
	    - Bug fix - when pack=true and ratio is set
	    - Bug fix - agraph/scan.l to work with latest flex beta

August 2, 2002
		- declare this version 1.8.9
	- general
	    - split rpm into:
	        graphviz, graphviz-tcl, graphviz-graphs, graphviz-devel
	    - gcc3 warning cleanup
	    - Install lincdt, libgraph, libagraph, libgd, libpathplan, libexp,
		and libpack so that they can be used by other programs.
		Headers and man3 in graphviz-devel
	- dynagraph, graphsearch
 	    - New tools based on libagraph and written in C++
	- dot, neato, twopi
	    - Add node and edge tooltips for use with -Tcmap
	    	\N,\E,\H,\T substitutions also work in tooltips.
	    - Add alt="label_string" to -Tcmap
	    - Add edge-label and port mappings to -Tps and -Tps2 so
	        that edges can be hyperlinked in PDF documents.
	    - Add support for \E (edge name), \H (head-node name),
	        \T (tail-node name) substitutions in edge labels and edge URLs
	    - Add support for stylesheet="file.css" for use in -Tsvg
	    - Fix -Tpic to work with recent gpic (Bruce Lilly)
	    - Fix alignment of imagemaps to images.
	    - Fix "transparent" color support in -Tsvg
	    - Fix support for graph [URL="default.html"] in -Tsvg and -Tcmap.
	    - Fix '&' escaping in URLs in -Tsvg
	    - Fix infinite loop in dot layout algorithm
	    - Fix text rotations again (hopefully freetype is stable now.)
	    - Cluster layout improvements
	    - Clean up warning messages from pathplan
	    - Consolidation of mapping code from imapgen.c and ismapgen.c into mapgen.c
	- gpr
	    - Added additional mode to extract components based sharing an
	        edge or a cluster
	    - Fix test for getopt
	- tcl-based tools
	    - Disable tcl-based tool building if tcl/tk not available
	        with stubs support.
	- documentation updates: FAQ, dotguide, dot.1
July 5, 2002
	    - declare 1.8.7 a "brown bag" release
		 and declare this version 1.8.8
	- remove wrong assert in gdgen.c
	- fix graph centering in bitmap outputs
	- provide enough margins
	- fix line widths after scaling
		(test with directed/proc3d.dot)
	- fix text rotations (requires libfreetype.so.6.3.1)
		(test with directed/NaN.dot)
July 5, 2002
	    - declare this version 1.8.7
	- Fix missing "]" in ihi demo.
July 2, 2002
	- Add URL mappings for clusters: svg,svgz,ps,ismap,imap,cmap.
	- Fix to avoid white edges in bitmap outputs when bgcolor is set.
	- Improve sizing and position of strings in bitmap outputs
	  when using builtin fonts (when font file not found).
	- Fix \N substitution in edge URLs in imap and cmap outputs.
	- Add -Tcmap for client-side imagemaps.
	- Generate warnings instead of access violation for EPSF file problems.
	- Various spline fixes in neato.
	- Fixes to pack.c
	- Add feature to ccomps to allow extraction of individual component
	  by number or node.
	- Cdt make to use iffe provided in the tools directory.
	- Various Makefile.old fixes.
	- Use HAVE_LIBZ to remove GD2 format if libz not available.
	  Now bare-bones programs can be built without any add-on libraries.
	- Modified dot grammar to allow simple name attributes in attribute
	  lists.  Thus, [splines] is equivalent to [splines=true]. Adopted
	  the same convention for command line attributes -G, -E and -N.
	  In addition, such command line attributes now override any
	  competing initial attribute statements.
	- HP-UX 11.11 build fixes for struct dioattr.
	- Fix for bug #158 "Nodes disappear with ports"
	- Various Windows-specific #ifdefs
	- Fix edge coordinates in -Tplain.

May 24, 2002
	    - declare this version 1.8.6
May 19, 2002
	- Fixed segfault from use of bgcolor in clusters.
May 15, 2002
	- Changed install location of architecture-independent demo
	  scripts and graphs to <prefix>/share/graphviz/ to conform to FHS.
	- Avoid multiple linking of libfreetype (and others) which caused
	  problems on SunOS-2.8.
May 6, 2002
	- Factored out some duplicated arrow code from dotgen/splines.c
	  and neatorgen/splines.c into common/arrows.c.
	- Added new arrow types:  halfopen, box, obox, crow.
	- Touched up the arrow designs so that they look better at default size.
	- Modified/extended graphs/directed/newarrows.dot to show new arrows.
May 3, 2002
        - Added some UML arrow types from Diomidis Spinellis <dds@aueb.gr>
	  empty, invempty, open, diamond, odiamond.
May 2, 2002
	- Added new pack option to neato. This causes each connected component
	  to be laid out separately, and then the resulting graphs are packed
	  together in a single layout.
	- Amended neato to accept new tee arrowhead.
April 19, 2002
	- Coords of rectangles changed to left/top right/bottom in -Timap.
	- Generate COPYING from LICENSE.html during ./authogen.sh,
	  remove COPYING from CVS.
April 16, 2002
	- Minor license file patches.
	- Corrected one of those reversed flat edge bugs again.

April 11, 2002
	     - declared this version 1.8.5
	- various portability fixes
	- various SVG fixes and optimizations
April 5, 2002:
	     - declared this version 1.8.4
	- SVG renderer:
		- make graph|node|edge ids unique, particularly for multiedges
		- put graph|node|edge names in <title>...</title>
		- use some property inheritance to reduce size of output
		- fix compile errors when no zlib
		- updated DTD reference
	- GD renderer:
		- Minimal Type1 font support:
			- look in /usr/lib/X11/fonts/Type1/
			- look for .pfa or .pfb font files based on fontname
		- run gdgen.c through dos2unix - problems with gcc on SuSE
	- fix Mac-OSX build problems:
		- improve strto[u]ll configure tests
		- add -fno-common for extern problem
		- function renamed to avoid conflicts (vis -> visibility)
		- add configure tests for search.h, malloc.h, getopt.h, errno.h
		- improve configure tests for FILE struct features
		- add configure tests for lrand48
	- add new demo graphs:
		- graphs/undirected/Heawood.dot
		- graphs/undirected/Petersen.dot
	- neato:
		- fix for -x implementation in neato (Bug 77)
		- fix spline problem (Bug 87)
		- fix some divide-by-zero problems
	- twopi:
		- fix Bug 117
		- update man pages for unconnected graphs capability
	- added arrowhead or arrowtail = tee
March 22, 2002:
	- add dotneato/pack code to twopi
	- add contrib/prune to gnu build and install
March 20, 2002:
	    - declared this version 1.8.3
	- fixed parse error for lines starting with '#' in .dot files
	- fixed a recently introduced bug that caused failure of:
		digraph G {  {rank = same;  A -> B; B -> A } }
	- updated DOCTYPE header in SVG outputs
	- added dotneato/common/xbuf.[ch] for dynamic string handling
	  to avoid sprintf buffer overruns.
	- twopigen - handle special case of graphs with < 3 nodes.
	- neato - handle point shapes
	- added fontcolor support to svg
March 14, 2002:
	- Fixed bug 109
	- Removed duplicate definitions for str[n]casecmp
	- Added missing declarations needed for Windows
	- Cleaned up warning messages from set but unused variables
	- Removed use of DOS preprocessor variable; uniformly replaced by MSWIN32
March 8, 2002:
	- declared this version 1.8.2
    - Mainly to fix a missed static buffer problem which trips up the
      Windows community
March 1, 2002:
	- declared this version 1.8.1
    - Bug fixes reported from user testing of 1.8.0, especially problem
      with SVG output
February 25, 2002:
	- updated dotguide.tex and moved to LaTeX article format
	- added webdot.cgi perl script, enhanced to accept the same
	    argument format as John's tcl version (so it can also
	    serve neato and twopi graph layouts).

February 7, 2002: graphviz-1.8.0 pre
	- declared this version 1.8.0

February 5, 2002: graphviz-1.7.17-0
    - various 64bit portability fixes
    - various bug fixes
January 2, 2002: graphviz-1.7.16-0
    - dotneato
	- fix bugs in -Tps output due to pen/fill color changes
	- various -Tfig.c fixes
	- various portability fixes
December 28, 2001: graphviz-1.7.15-0
    -dotneato
        - introduce damping factor into neato's solver
        - clean up pencolor v fillcolor code so that filled polygons are drawn
		just once if the renderer is capable (e.g. svg, fig)
        - complete -Tfig support (xfig format)
December 11, 2001: graphviz-1.7.14-0
    -dotneato
	- add -Tsvgz (compressed SVG) support
December 11, 2001: graphviz-1.7.13-0
    - dotneato
        - fontwidth fixes
	- remove some potential buffer overruns
	- escape '&' in SVG, unless it is already part of a UTF entity sequence
	- recognize Times_New_Roman and Courier_New as default font names.
	- improve -liconv support in configure
	- clean up some compiler warnings
    - dynagraph
	- change "round" to "ROUND" to avoid conflict with system headers on linux
December 03, 2001: graphviz-1.7.12-0
    - dotneato
        - add -Tplain-ext which includes port identifiers edge records
	- escape '>' with '&gt;' in edge ids and edge URLs in -Tsvg.
	- spline fixes
	- mincross fixes
	- improved text alignment in nodes - particularly in bitmap outputs.
	- fixed text scaling problems for 8-bit characters (e.g. umlauts)
	- add graph lexer and postscript support for extended characters
    - lefty
        - fix for X11 displays
    - pathplan
        - added workaround for gcc-0.96 bug when "-O2 -mcpu=686 -ffast-math"
October 22, 2001: graphviz-1.7.11-0
    - dotneato
	- svg - fix landscape "y" direction
	      - fix text rotation (works in batik, not yet in sodipodi or amaya)
	      - fix linewidth
	      - fix xmnls:xlink reference
    - doc
	- Dot.ref - updated
    - graphs/directed
        - newarrows.dot expanded
	- honda-tokoro.dot added
October 21, 2001: graphviz-1.7.10-0
    - lefty & dotty
	- realign code with EK's master tree.
	  includes fix for dirty trails when dragging nodes in dotty.
    - dotneato
	- svg - kludge escape of "<" & ">" characters in labels.
    - general
	- generate doxygen documentation on http://www.graphviz.org/
August 20, 2001: graphviz-1.7.9-0
    - general
	- first release from relocated cvs server
    - dotneato
        - fix for abort from spline code
        - fix for crash from gd tiling code
August 15, 2001: graphviz-1.7.8-0
    - general
        - Update gd to gd-2.0.1 with extensions
    - dotneato
        - more spline fixes
        - add suport for "#rgb" color specification
        - add twopi layout engine (circular layouts)
July 13, 2001: graphviz-1.7.7-0
    - Synchronization release prior to relocating CVS server.
    - general
    	- some Makefile fixes for OpenBSD
	- some FAQ updates
    - dotneato
        - self-edge fixes
        - spline fixes
    - libgraph
        - parser fixes
July 1, 2001: graphviz-1.7.6-3
    - general
	- portability fixes (including 14 charater file names !)
	- memory leak fixes
	- "make test" targets in graphs/directed, graphs/undirected
    - configure
	- add support for building without X11, Tk, Tcl
	- add hooks for dmalloc and ElectricFence debugging
    - dotneato
	- spline fixes
	- cluster fixes
	- fix label centering
	- fix support for graph margins in bitmapped outputs
	- correction to PostScript preamble
	- SVG generator improvement - now works with Amaya and SodiPodi
    - tcldot
	- now uses Tcl Channels properly for input
	- fixes for linewidth support
	- command extensions
	    - listattributes now accepts list
	    - queryattributes now accepts list
	    - setattributes now accepts list
	    - queryattributevalues - new command
		- generates list of pairs compatible with setattributes
    - dotty
	- passthrough keyboard events
    - doted
	- fix resizing problems
	- add PNG and SVG output formats

April 27, 2001: graphviz-1.7.6

    NEW FEATURES

    Added a collection of graph processing tools:

    acyclic : a filter that takes a directed graph as input
    and outputs a copy of the graph with sufficient edges
    reversed to make the graph acyclic.

    ccomps : decomposes graphs into their connected components,
    printing the components to standard output.

    colorize : is a filter that sets node colors from initial
    seed values. Colors flow along edges from tail to head.

    gc : a graph analogue to wc in that it prints to standard
    output the number of nodes, edges, connected components or
    clusters contained in the input files.

    gpr : a graph stream editor inspired by awk. It copies
    input graphs to its output, possibly transforming their
    structure and attributes, creating new graphs, or
    printing arbitrary information.

    nop : reads a stream of graphs and prints each in
    pretty-printed (canonical) format on stdout.

    sccmap : decomposes digraphs into strongly connected components
    and an auxiliary map of the relationship between components.

    tred : computes the transitive reduction of directed graphs,
    and prints the resulting graphs to standard output. This
    removes edges implied by transitivity.

    unflatten : is a preprocessor to dot that is used to improve
    the aspect ratio of graphs having many leaves or disconnected
    nodes. The usual layout for such a graph is generally very
    wide or tall. unflatten inserts invisible edges or adjusts
    the minlen on edges to improve layout compaction.


    FIXES

    Add FAQ

    Change PNG default background color from transparent to white
    because of the difficulty some viewers have with transparency.

    Add support for [color=transparent]

    Fix broken support for specific capitalized fontnames
    (Times Helvetica Arial Courier)

    Fix broken support for DOTFONTPATH

    Some bitmap font scaling fixes - we're still not happy with
    bitmap font scaling as some labels still exceed the area
    allocated by the layout engines.

    Some -Timap fixes for mouse sensitive graphs on web pages

    Some cluster layout fixes

    Fix for [rankdir=LR] problems when using neato layout engine

    Some neato layout fixes

    Updates to unix.dot

    Various OS and distro fixes


December 23, 2000: graphviz-1.7.5

   - update to gd-1.8.4 and freetype2
   - add support for font paths


December 15, 2000: graphviz-1.7.4
    -various cluster fixes
    -separate support for node fillcolor from pencolor (see dot.1)
    -add support for dotted and dashed lines to bitmap renderers (PNG, GIF etc)
    -add support for varying linewidth to bitmap renderers
    -remove libtcldot dependence on lingdtclft (already statically included)
    -various fixes to build processes, GNU and non-GNU


graphviz-1.7.3 .....

May 3, 2000: removed webdot into its own CVS module and rpm package

April 16, 2000: Use check for "gdImagePng" to make sure that we have
   recent version of libgd.  <ellson@graphviz.org>

April 14, 2000: Add Tcldgl and dge demo <ellson@graphviz.org>

April 14, 2000: Add dynagraph libraries <north@research.att.com>

April 14, 2000: Flatten directory hierarchy of sources <ellson@graphviz.org>

April 14, 2000: Fix X11 library detection for lefty:
	src/configure.in, src/lefty/Makefile.in
   <ellson@graphviz.org>

April 14, 2000: Fix pic support:
	src/dotneato/picgen.c,
	src/dotneato/emit.c,
	webdot/tcl/webdot.tcl
   <Bruce Lilly>

April 7, 2000: Upgrade webdot installation process:
	webdot/Makefile, webdot/README
    <ellson@graphviz.org>

March 13, 2000: Support for virtual hosts in webdot/webdot.tcl, add
   "puts $skt "Host: $server"     Michael Tillberg <mt@proteome.com>

March 13, 2000: Fix to src/graph/parser.y line 149
   "if ((e->head == t->node) && !(Agraph_type & AGDIGRAPH)) {"
   Stephen North  <north@research.att.com>

March 13, 2000: Use AM_PROG_LIBTOOL instead of AC_PROG_LIBTOOL
   in configure.in.  John Ellson <ellson@graphviz.org>
```

[11.0.0]: https://gitlab.com/graphviz/graphviz/compare/10.0.1...11.0.0
[10.0.1]: https://gitlab.com/graphviz/graphviz/compare/9.0.0...10.0.1
[9.0.0]: https://gitlab.com/graphviz/graphviz/compare/8.1.0...9.0.0
[8.1.0]: https://gitlab.com/graphviz/graphviz/compare/8.0.5...8.1.0
[8.0.5]: https://gitlab.com/graphviz/graphviz/compare/8.0.3...8.0.5
[8.0.3]: https://gitlab.com/graphviz/graphviz/compare/8.0.2...8.0.3
[8.0.2]: https://gitlab.com/graphviz/graphviz/compare/8.0.1...8.0.2
[8.0.1]: https://gitlab.com/graphviz/graphviz/compare/7.1.0...8.0.1
[7.1.0]: https://gitlab.com/graphviz/graphviz/compare/7.0.6...7.1.0
[7.0.6]: https://gitlab.com/graphviz/graphviz/compare/7.0.5...7.0.6
[7.0.5]: https://gitlab.com/graphviz/graphviz/compare/7.0.4...7.0.5
[7.0.4]: https://gitlab.com/graphviz/graphviz/compare/7.0.3...7.0.4
[7.0.3]: https://gitlab.com/graphviz/graphviz/compare/7.0.2...7.0.3
[7.0.2]: https://gitlab.com/graphviz/graphviz/compare/7.0.1...7.0.2
[7.0.1]: https://gitlab.com/graphviz/graphviz/compare/7.0.0...7.0.1
[7.0.0]: https://gitlab.com/graphviz/graphviz/compare/6.0.2...7.0.0
[6.0.2]: https://gitlab.com/graphviz/graphviz/compare/6.0.1...6.0.2
[6.0.1]: https://gitlab.com/graphviz/graphviz/compare/5.0.1...6.0.1
[5.0.1]: https://gitlab.com/graphviz/graphviz/compare/5.0.0...5.0.1
[5.0.0]: https://gitlab.com/graphviz/graphviz/compare/4.0.0...5.0.0
[4.0.0]: https://gitlab.com/graphviz/graphviz/compare/3.0.0...4.0.0
[3.0.0]: https://gitlab.com/graphviz/graphviz/compare/2.50.0...3.0.0
[2.50.0]: https://gitlab.com/graphviz/graphviz/compare/2.49.3...2.50.0
[2.49.3]: https://gitlab.com/graphviz/graphviz/compare/2.49.2...2.49.3
[2.49.2]: https://gitlab.com/graphviz/graphviz/compare/2.49.1...2.49.2
[2.49.1]: https://gitlab.com/graphviz/graphviz/compare/2.49.0...2.49.1
[2.49.0]: https://gitlab.com/graphviz/graphviz/compare/2.48.0...2.49.0
[2.48.0]: https://gitlab.com/graphviz/graphviz/compare/2.47.3...2.48.0
[2.47.3]: https://gitlab.com/graphviz/graphviz/compare/2.47.2...2.47.3
[2.47.2]: https://gitlab.com/graphviz/graphviz/compare/2.47.1...2.47.2
[2.47.1]: https://gitlab.com/graphviz/graphviz/compare/2.47.0...2.47.1
[2.47.0]: https://gitlab.com/graphviz/graphviz/compare/2.46.1...2.47.0
[2.46.1]: https://gitlab.com/graphviz/graphviz/compare/2.46.0...2.46.1
[2.46.0]: https://gitlab.com/graphviz/graphviz/compare/2.44.1...2.46.0
[2.44.1]: https://gitlab.com/graphviz/graphviz/compare/2.44.0...2.44.1
[2.44.0]: https://gitlab.com/graphviz/graphviz/compare/2.42.4...2.44.0
[2.42.4]: https://gitlab.com/graphviz/graphviz/compare/2.42.3...2.42.4
[2.42.3]: https://gitlab.com/graphviz/graphviz/compare/2.42.2...2.42.3
