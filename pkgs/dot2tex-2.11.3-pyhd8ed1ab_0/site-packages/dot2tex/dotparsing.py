# -*- coding: utf-8 -*-
"""Graphviz dot language parser.

This parser is derived from the the parser distributed with the pydot module.

Original authors are
Michael Krause <michael AT krause-software.de>
Ero Carrera <ero AT dkbza.org>
"""

__version__ = '2.11.3'
__author__ = ['Michael Krause', 'Ero Carrera', 'Kjell Magne Fauske']
__license__ = 'MIT'

import re
import itertools
import os
import logging
import string

import pyparsing
from pyparsing import __version__ as pyparsing_version
from pyparsing import (Literal, CaselessLiteral, Word, OneOrMore, Forward, Group, Optional, Combine, restOfLine,
                       cStyleComment, nums, alphanums,
                       ParseException, CharsNotIn, Suppress, Regex, removeQuotes)

from collections import OrderedDict

dot_keywords = ['graph', 'subgraph', 'digraph', 'node', 'edge', 'strict']

id_re_alpha_nums = re.compile('^[_a-zA-Z][a-zA-Z0-9_]*$')
id_re_num = re.compile('^-?(\.[0-9]+|[0-9]+(\.[0-9]*)?)$')
id_re_with_port = re.compile('^.*:([^"]+|[^"]*\"[^"]*\"[^"]*)$')
id_re_dbl_quoted = re.compile('^\".*\"$', re.S)
id_re_html = re.compile('^<<.*>>$', re.S)

log = logging.getLogger("dot2tex")


def needs_quotes(s):
    """Checks whether a string is a dot language ID.

    It will check whether the string is solely composed
    by the characters allowed in an ID or not.
    If the string is one of the reserved keywords it will
    need quotes too.
    """

    if s in dot_keywords:
        return True

    chars = [ord(c) for c in s if ord(c) > 0x7f or ord(c) == 0]
    if chars:
        return True

    res = id_re_alpha_nums.match(s)
    if not res:
        res = id_re_num.search(s)
        if not res:
            # res = id_re_dbl_quoted.match(s)
            if not res:
                res = id_re_html.match(s)
            pass

            ##                if not res:
            ##                    res = id_re_with_port.match(s)

    if not res:
        return True

    return False


def quote_if_necessary(s):
    if not isinstance(s, str):
        return s
    tmp = s
    if needs_quotes(tmp):
        tmp = '"%s"' % s  # .replace('"','\\"')
    tmp = tmp.replace('<<', '<')
    tmp = tmp.replace('>>', '>')
    return tmp


def flatten(lst):
    for elem in lst:
        if type(elem) in (tuple, list):
            for i in flatten(elem):
                yield i
        else:
            yield elem


# Code snippet from Python Cookbook, 2nd Edition by David Ascher, Alex Martelli
# and Anna Ravenscroft; O'Reilly 2005
def windows(iterable, length=2, overlap=0, padding=True):
    it = iter(iterable)
    results = list(itertools.islice(it, length))
    while len(results) == length:
        yield results
        results = results[length - overlap:]
        results.extend(itertools.islice(it, length - overlap))
    if padding and results:
        results.extend(itertools.repeat(None, length - len(results)))
        yield results


def nsplit(seq, n=2):
    """Split a sequence into pieces of length n

    If the lengt of the sequence isn't a multiple of n, the rest is discareded.
    Note that nsplit will strings into individual characters.

    Examples:
    >>> nsplit('aabbcc')
    [('a', 'a'), ('b', 'b'), ('c', 'c')]
    >>> nsplit('aabbcc',n=3)
    [('a', 'a', 'b'), ('b', 'c', 'c')]

    # Note that cc is discarded
    >>> nsplit('aabbcc',n=4)
    [('a', 'a', 'b', 'b')]
    """
    return [xy for xy in zip(*[iter(seq)] * n)]


# The following function is from the pydot project
def __find_executables(path):
    """Used by find_graphviz

    path - single directory as a string

    If any of the executables are found, it will return a dictionary
    containing the program names as keys and their paths as values.

    Otherwise returns None
    """

    success = False
    progs = {'dot': '', 'twopi': '', 'neato': '', 'circo': '', 'fdp': ''}

    was_quoted = False
    path = path.strip()
    if path.startswith('"') and path.endswith('"'):
        path = path[1:-1]
        was_quoted = True

    if os.path.isdir(path):
        for prg in progs:
            if progs[prg]:
                continue

            if os.path.exists(os.path.join(path, prg)):
                if was_quoted:
                    progs[prg] = '"' + os.path.join(path, prg) + '"'
                else:
                    progs[prg] = os.path.join(path, prg)

                success = True

            elif os.path.exists(os.path.join(path, prg + '.exe')):
                if was_quoted:
                    progs[prg] = '"' + os.path.join(path, prg + '.exe') + '"'
                else:
                    progs[prg] = os.path.join(path, prg + '.exe')

                success = True

    if success:
        return progs

    else:
        return None


# The following function is from the pydot project
# The multi-platform version of this 'find_graphviz' function was
# contributed by Peter Cock
#
def find_graphviz():
    """Locate Graphviz's executables in the system.

    Tries three methods:

    First: Windows Registry (Windows only)
    This requires Mark Hammond's pywin32 is installed.

    Secondly: Search the path
    It will look for 'dot', 'twopi' and 'neato' in all the directories
    specified in the PATH environment variable.

    Thirdly: Default install location (Windows only)
    It will look for 'dot', 'twopi' and 'neato' in the default install
    location under the "Program Files" directory.

    It will return a dictionary containing the program names as keys
    and their paths as values.

    If this fails, it returns None.
    """

    # Method 1 (Windows only)
    #
    if os.sys.platform == 'win32':
        try:
            import win32api, win32con

            # Get the GraphViz install path from the registry
            #
            hkey = win32api.RegOpenKeyEx(win32con.HKEY_LOCAL_MACHINE,
                                         "SOFTWARE\AT&T Research Labs\Graphviz", 0, win32con.KEY_QUERY_VALUE)

            path = win32api.RegQueryValueEx(hkey, "InstallPath")[0]
            win32api.RegCloseKey(hkey)

            # Now append the "bin" subdirectory:
            #
            path = os.path.join(path, "bin")
            progs = __find_executables(path)
            if progs is not None:
                # print("Used Windows registry")
                return progs

        except ImportError:
            # Print a messaged suggesting they install these?
            #
            log.debug('The win32api is not installed')
            pass
        except:
            log.debug('Failed to access the registry key')

    # Method 2 (Linux, Windows etc)
    #
    if 'PATH' in os.environ:
        for path in os.environ['PATH'].split(os.pathsep):
            progs = __find_executables(path)
            if progs is not None:
                return progs

    # Method 3 (Windows only)
    #
    if os.sys.platform == 'win32':
        # Try and work out the equivalent of "C:\Program Files" on this
        # machine (might be on drive D:, or in a different language)
        #
        if 'PROGRAMFILES' in os.environ:
            # Note, we could also use the win32api to get this
            # information, but win32api may not be installed.

            path = os.path.join(os.environ['PROGRAMFILES'], 'ATT', 'GraphViz', 'bin')

        else:
            # Just in case, try the default...
            path = r"C:\Program Files\att\Graphviz\bin"

        progs = __find_executables(path)

        if progs is not None:
            # print("Used default install location")
            return progs

    for path in (
            '/usr/bin', '/usr/local/bin',
            '/opt/local/bin',
            '/opt/bin', '/sw/bin', '/usr/share',
            '/Applications/Graphviz.app/Contents/MacOS/'):
        progs = __find_executables(path)
        if progs is not None:
            # print("Used path")
            return progs

    # Failed to find GraphViz
    #
    return None


ADD_NODE = 'add_node'
ADD_EDGE = 'add_edge'
ADD_GRAPH_TO_NODE_EDGE = 'add_graph_to_node_edge'
ADD_NODE_TO_GRAPH_EDGE = 'add_node_to_graph_edge'
ADD_GRAPH_TO_GRAPH_EDGE = 'add_graph_to_graph_edge'
ADD_SUBGRAPH = 'add_subgraph'
SET_DEF_NODE_ATTR = 'set_def_node_attr'
SET_DEF_EDGE_ATTR = 'set_def_edge_attr'
SET_DEF_GRAPH_ATTR = 'set_def_graph_attr'
SET_GRAPH_ATTR = 'set_graph_attr'


class DotDataParser(object):
    """Container class for parsing Graphviz dot data"""

    def __init__(self):
        pass
        self.dotparser = self.define_dot_parser()

    # parse actions
    def _proc_node_id(self, toks):
        if len(toks) > 1:
            return toks[0], toks[1]
        else:
            return toks

    def _proc_attr_list(self, toks):
        return dict(nsplit(toks, 2))

    def _proc_attr_list_combine(self, toks):
        if toks:
            first_dict = toks[0]
            for d in toks:
                first_dict.update(d)

            return first_dict
        return toks

    def _proc_attr_assignment(self, toks):
        return SET_GRAPH_ATTR, dict(nsplit(toks, 2))

    def _proc_node_stmt(self, toks):
        """Return (ADD_NODE, node_name, options)"""
        if len(toks) == 2:
            return tuple([ADD_NODE] + list(toks))
        else:
            return tuple([ADD_NODE] + list(toks) + [{}])

    def _proc_edge_stmt(self, toks):
        """Return (ADD_EDGE, src, dest, options)"""
        edgelist = []
        opts = toks[-1]
        if not isinstance(opts, dict):
            opts = {}
        for src, op, dest in windows(toks, length=3, overlap=1, padding=False):
            # is src or dest a subgraph?
            srcgraph = destgraph = False
            if len(src) > 1 and src[0] == ADD_SUBGRAPH:
                edgelist.append(src)
                srcgraph = True
            if len(dest) > 1 and dest[0] == ADD_SUBGRAPH:
                edgelist.append(dest)
                destgraph = True
            if srcgraph or destgraph:
                if srcgraph and destgraph:
                    edgelist.append((ADD_GRAPH_TO_GRAPH_EDGE, src[1], dest[1], opts))
                elif srcgraph:
                    edgelist.append((ADD_GRAPH_TO_NODE_EDGE, src[1], dest, opts))
                else:
                    edgelist.append((ADD_NODE_TO_GRAPH_EDGE, src, dest[1], opts))
            else:
                # ordinary edge
                edgelist.append((ADD_EDGE, src, dest, opts))

        return edgelist

    def _proc_default_attr_stmt(self, toks):
        """Return (ADD_DEFAULT_NODE_ATTR,options"""
        if len(toks) == 1:
            gtype = toks
            attr = {}
        else:
            gtype, attr = toks
        if gtype == 'node':
            return SET_DEF_NODE_ATTR, attr
        elif gtype == 'edge':
            return SET_DEF_EDGE_ATTR, attr
        elif gtype == 'graph':
            return SET_DEF_GRAPH_ATTR, attr
        else:
            return 'unknown', toks

    def _proc_subgraph_stmt(self, toks):
        """Returns (ADD_SUBGRAPH, name, elements)"""
        return 'add_subgraph', toks[1], toks[2].asList()

    def _main_graph_stmt(self, toks):
        return toks[0], toks[1], toks[2], toks[3].asList()

    # The dot grammar is based on the dot parser from the pydot project.
    def define_dot_parser(self):
        """Define dot grammar

        Based on the grammar http://www.graphviz.org/doc/info/lang.html
        """
        # punctuation
        colon = Literal(":")
        lbrace = Suppress("{")
        rbrace = Suppress("}")
        lbrack = Suppress("[")
        rbrack = Suppress("]")
        lparen = Literal("(")
        rparen = Literal(")")
        equals = Suppress("=")
        comma = Literal(",")
        dot = Literal(".")
        slash = Literal("/")
        bslash = Literal("\\")
        star = Literal("*")
        semi = Suppress(";")
        at = Literal("@")
        minus = Literal("-")
        pluss = Suppress("+")

        # keywords
        strict_ = CaselessLiteral("strict")
        graph_ = CaselessLiteral("graph")
        digraph_ = CaselessLiteral("digraph")
        subgraph_ = CaselessLiteral("subgraph")
        node_ = CaselessLiteral("node")
        edge_ = CaselessLiteral("edge")

        punctuation_ = "".join([c for c in string.punctuation if c not in '_']) + string.whitespace
        # token definitions

        identifier = Word(alphanums + "_").setName("identifier")

        # double_quoted_string = QuotedString('"', multiline=True,escChar='\\',
        #    unquoteResults=True) # dblQuotedString
        double_quoted_string = Regex(r'\"(?:\\\"|\\\\|[^"])*\"', re.MULTILINE)
        double_quoted_string.setParseAction(removeQuotes)
        quoted_string = Combine(double_quoted_string +
                                Optional(OneOrMore(pluss + double_quoted_string)), adjacent=False)
        alphastring_ = OneOrMore(CharsNotIn(punctuation_))

        def parse_html(s, loc, toks):
            return '<<%s>>' % ''.join(toks[0])

        opener = '<'
        closer = '>'
        try:
            html_text = pyparsing.nestedExpr(opener, closer,
                                             ((CharsNotIn(
                                                     opener + closer).setParseAction(lambda t: t[0]))
                                              )).setParseAction(parse_html)
        except:
            log.debug('nestedExpr not available.')
            log.warning('Old version of pyparsing detected. Version 1.4.8 or '
                        'later is recommended. Parsing of html labels may not '
                        'work properly.')
            html_text = Combine(Literal("<<") + OneOrMore(CharsNotIn(",]")))

        float_number = Combine(Optional(minus) +
                               OneOrMore(Word(nums + "."))).setName("float_number")

        ID = (alphastring_ | html_text | float_number |
              quoted_string |  # .setParseAction(strip_quotes) |
              identifier).setName("ID")



        righthand_id = (float_number | ID).setName("righthand_id")

        port_angle = (at + ID).setName("port_angle")

        port_location = ((OneOrMore(Group(colon + ID)) |
                          Group(colon + lparen + ID + comma + ID + rparen))).setName("port_location")

        port = Combine((Group(port_location + Optional(port_angle)) |
                        Group(port_angle + Optional(port_location)))).setName("port")

        node_id = (ID + Optional(port))
        a_list = OneOrMore(ID + Optional(equals + righthand_id) +
                           Optional(comma.suppress())).setName("a_list")

        attr_list = OneOrMore(lbrack + Optional(a_list) +
                              rbrack).setName("attr_list").setResultsName('attrlist')

        attr_stmt = ((graph_ | node_ | edge_) + attr_list).setName("attr_stmt")

        edgeop = (Literal("--") | Literal("->")).setName("edgeop")

        stmt_list = Forward()
        graph_stmt = (lbrace + Optional(stmt_list) +
                      rbrace + Optional(semi)).setName("graph_stmt")

        edge_point = Forward()

        edgeRHS = OneOrMore(edgeop + edge_point)
        edge_stmt = edge_point + edgeRHS + Optional(attr_list)

        subgraph = (Optional(subgraph_, '') + Optional(ID, '') + Group(graph_stmt)).setName("subgraph").setResultsName(
                'ssubgraph')

        edge_point <<= (subgraph | graph_stmt | node_id)

        node_stmt = (node_id + Optional(attr_list) + Optional(semi)).setName("node_stmt")

        assignment = (ID + equals + righthand_id).setName("assignment")
        stmt = (assignment | edge_stmt | attr_stmt | subgraph | graph_stmt | node_stmt).setName("stmt")
        stmt_list <<= OneOrMore(stmt + Optional(semi))

        graphparser = ((Optional(strict_, 'notstrict') + ((graph_ | digraph_)) +
                        Optional(ID, '') + lbrace + Group(Optional(stmt_list)) + rbrace).setResultsName("graph"))

        singleLineComment = Group("//" + restOfLine) | Group("#" + restOfLine)

        # actions
        graphparser.ignore(singleLineComment)
        graphparser.ignore(cStyleComment)
        node_id.setParseAction(self._proc_node_id)
        assignment.setParseAction(self._proc_attr_assignment)
        a_list.setParseAction(self._proc_attr_list)
        edge_stmt.setParseAction(self._proc_edge_stmt)
        node_stmt.setParseAction(self._proc_node_stmt)
        attr_stmt.setParseAction(self._proc_default_attr_stmt)
        attr_list.setParseAction(self._proc_attr_list_combine)
        subgraph.setParseAction(self._proc_subgraph_stmt)
        # graph_stmt.setParseAction(self._proc_graph_stmt)
        graphparser.setParseAction(self._main_graph_stmt)
        return graphparser

    def build_graph(self, graph, tokens):
        subgraph = None
        for element in tokens:
            cmd = element[0]
            if cmd == ADD_NODE:
                cmd, nodename, opts = element
                node = graph.add_node(nodename, **opts)
                graph.allitems.append(node)

            elif cmd == ADD_EDGE:
                cmd, src, dest, opts = element
                srcport = destport = ""
                if isinstance(src, tuple):
                    srcport = src[1]
                    src = src[0]
                if isinstance(dest, tuple):
                    destport = dest[1]
                    dest = dest[0]
                edge = graph.add_edge(src, dest, srcport, destport, **opts)
                graph.allitems.append(edge)
            elif cmd in [ADD_GRAPH_TO_NODE_EDGE, ADD_GRAPH_TO_GRAPH_EDGE, ADD_NODE_TO_GRAPH_EDGE]:
                cmd, src, dest, opts = element
                srcport = destport = ""
                if isinstance(src, tuple):
                    srcport = src[1]

                if isinstance(dest, tuple):
                    destport = dest[1]
                if not (cmd == ADD_NODE_TO_GRAPH_EDGE):
                    if cmd == ADD_GRAPH_TO_NODE_EDGE:
                        src = subgraph
                    else:
                        src = prev_subgraph
                        dest = subgraph
                else:
                    dest = subgraph

                edges = graph.add_special_edge(src, dest, srcport, destport, **opts)
                graph.allitems.extend(edges)

            elif cmd == SET_GRAPH_ATTR:
                graph.set_attr(**element[1])

            elif cmd == SET_DEF_NODE_ATTR:
                graph.add_default_node_attr(**element[1])
                defattr = DotDefaultAttr('node', **element[1])
                graph.allitems.append(defattr)
            elif cmd == SET_DEF_EDGE_ATTR:
                graph.add_default_edge_attr(**element[1])
                defattr = DotDefaultAttr('edge', **element[1])
                graph.allitems.append(defattr)
            elif cmd == SET_DEF_GRAPH_ATTR:
                graph.add_default_graph_attr(**element[1])
                defattr = DotDefaultAttr('graph', **element[1])
                graph.allitems.append(defattr)
                graph.attr.update(**element[1])
            elif cmd == ADD_SUBGRAPH:
                cmd, name, elements = element
                # print("Adding subgraph")
                if subgraph:
                    prev_subgraph = subgraph
                subgraph = graph.add_subgraph(name)
                subgraph = self.build_graph(subgraph, elements)
                graph.allitems.append(subgraph)

        return graph

    def build_top_graph(self, tokens):
        """Build a DotGraph instance from parsed data"""
        # get basic graph information
        strict = tokens[0] == 'strict'
        graphtype = tokens[1]
        directed = graphtype == 'digraph'
        graphname = tokens[2]
        # lets build the graph
        graph = DotGraph(graphname, strict, directed)
        self.graph = self.build_graph(graph, tokens[3])

    def parse_dot_data(self, data):
        """Parse dot data and return a DotGraph instance"""
        try:
            try:
                self.dotparser.parseWithTabs()
            except:
                log.warning('Old version of pyparsing. Parser may not work correctly')
            if os.sys.version_info[0] >= 3 and isinstance(data, bytes):
                data = data.decode()
            ndata = data.replace('\\\n', '')
            # lines = data.splitlines()
            # lines = [l.rstrip('\\') for l in lines]
            tokens = self.dotparser.parseString(ndata)
            self.build_top_graph(tokens[0])
            return self.graph

        except ParseException as err:
            # print(err.line)
            # print(" "*(err.column-1) + "^")
            # print(err)
            # return None
            raise

    def parse_dot_data_debug(self, data):
        """Parse dot data"""
        try:
            try:
                self.dotparser.parseWithTabs()
            except:
                log.warning('Old version of pyparsing. Parser may not work correctly')

            tokens = self.dotparser.parseString(data)
            self.build_top_graph(tokens[0])

            return tokens[0]

        except ParseException as err:
            print(err.line)
            print(" " * (err.column - 1) + "^")
            print(err)
            return None


class DotDefaultAttr(object):
    def __init__(self, element_type, **kwds):
        self.element_type = element_type
        self.attr = kwds

    def __str__(self):
        attrstr = ",".join(["%s=%s" % \
                            (quote_if_necessary(key), quote_if_necessary(val)) \
                            for key, val in self.attr.items()])
        if attrstr:
            attrstr = "[%s]" % attrstr
            return "%s%s;\n" % (self.element_type, attrstr)
        else:
            return ""


class DotParsingException(Exception):
    """Base class for dotparsing exceptions."""


class DotNode(object):
    """Class representing a DOT node"""

    def __init__(self, name, **kwds):
        """Create a Node instance

        Input:
            name - name of node. Have to be a string
            **kwds node attributes

        """
        self.name = name
        self.attr = {}
        self.parent = None
        self.attr.update(kwds)

    def __str__(self):
        attrstr = ",".join(["%s=%s" % \
                            (quote_if_necessary(key), quote_if_necessary(val)) \
                            for key, val in self.attr.items()])
        if attrstr:
            attrstr = "[%s]" % attrstr
        return "%s%s;\n" % (quote_if_necessary(self.name), attrstr)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        try:
            return self.name == other
        except:
            return False

    def __ne__(self, other):
        try:
            return self.name != other
        except:
            return False

    def __getattr__(self, name):
        try:
            return self.attr[name]
        except KeyError:
            raise AttributeError


class DotGraph(object):
    """Class representing a DOT graph"""

    def __init__(self, name='G', strict=True, directed=False, **kwds):
        self._nodes = OrderedDict()
        self._allnodes = {}
        self._alledges = {}
        self._allgraphs = []
        self._edges = {}
        self.strict = strict
        self.directed = directed
        self.subgraphs = []
        self.name = name
        self.padding = "    "
        self.seq = 0
        self.allitems = []

        self.attr = {}
        self.strict = strict
        self.level = 0
        self.parent = None
        self.root = self
        self.adj = {}

        self.default_node_attr = {}
        self.default_edge_attr = {}
        self.default_graph_attr = {}

        self.attr.update(kwds)
        self._allgraphs.append(self)
        pass

    def __len__(self):
        return len(self._nodes) + sum(len(s) for s in self.subgraphs)

    def __getattr__(self, name):
        try:
            return self.attr[name]
        except KeyError:
            raise AttributeError

    def get_name(self):
        if self.name.strip():
            return quote_if_necessary(self.name)
        else:
            return ""

    def add_node(self, node, **kwds):
        if not isinstance(node, DotNode):
            node = DotNode(str(node), **kwds)
        n = node.name

        ##        if n not in self.adj:
        ##            self.adj[n]={}
        if n in self._allnodes:
            self._allnodes[n].attr.update(kwds)
        else:
            node.attr.update(self.default_node_attr)
            node.attr.update(kwds)
            self._allnodes[n] = node
        if n not in self._nodes:
            self._nodes[n] = node

            # Todo: Adding a node to a subgraph should insert it in parent graphs
        ##        parent = self.parent
        ##        if parent:
        ##            print("Parent %s " % parent.name)
        ##        while parent:
        ##            print("Parent %s " % parent.name)
        ##            if n not in parent._nodes:
        ##                parent._nodes[n] = node
        ##                parent = parent.parent
        ##            else:
        ##                parent = None

        return node

    def add_edge(self, src, dst, srcport="", dstport="", **kwds):
        u = self.add_node(src)
        v = self.add_node(dst)
        edge = DotEdge(u, v, self.directed, srcport, dstport, **self.default_edge_attr)
        edge.attr.update(kwds)

        ##        if not self.strict:
        ##            self.adj[u][v]=self.adj[u].get(v,[])+ [edge]
        ##            if not self.directed:
        ##                self.adj[v][u]=self.adj[v].get(u,[])+ [edge]
        ##
        ##        else:
        ##            self.adj[u][v]=edge
        ##            if not self.directed:
        ##                self.adj[v][u]=edge

        edgekey = (u.name, v.name)

        if edgekey in self._alledges:
            edgs = self._alledges[edgekey]
            if not self.strict:
                # edge.parent = edge_parent
                if edgekey in self._edges:
                    self._edges[edgekey].append(edge)
                edgs.append(edge)

                ##            else:
                ##                edgs[0].attributes.update(edge.attributes)
                ##                return edgs[0]
        else:
            ##            edge.parent = edge_parent
            # edge.attr.update(self.default_edge_attr)
            self._alledges[edgekey] = [edge]
            self._edges[edgekey] = [edge]
        return edge

    def add_special_edge(self, src, dst, srcport="", dstport="", **kwds):
        src_is_graph = isinstance(src, DotSubGraph)
        dst_is_graph = isinstance(dst, DotSubGraph)
        edges = []
        if src_is_graph:
            src_nodes = src.get_all_nodes()
        else:
            src_nodes = [src]
        if dst_is_graph:
            dst_nodes = dst.get_all_nodes()
        else:
            dst_nodes = [dst]

        for src_node in src_nodes:
            for dst_node in dst_nodes:
                edge = self.add_edge(src_node, dst_node, srcport, dstport, **kwds)
                edges.append(edge)

        return edges

    def add_default_node_attr(self, **kwds):
        self.default_node_attr.update(kwds)

    def add_default_edge_attr(self, **kwds):
        self.default_edge_attr.update(kwds)

    def add_default_graph_attr(self, **kwds):
        self.default_graph_attr.update(kwds)

    ##            #nodecls = self._allnodes[name]
    ##            #nodeparent = nodecls.parent
    ##
    ##            parent = self.parent
    ##            if parent and (not (nodeparent == self)):
    ##                while parent:
    ##
    ##                    if name in parent._nodes:
    ##
    ##                        del parent._nodes[name]
    ##                        # changing a node parent may trigger a change in
    ##                        # edge parent
    ##                        nodecls.parent = self
    ##                        parent = None
    ##                        self._nodes[name] = nodecls
    ##                    else:
    ##                        parent = parent.parent
    ##
    ##
    ##
    ##
    ##        else:
    ##            nodecls.parent = self
    ##            self._allnodes[name]=nodecls
    ##            self._nodes[name] = nodecls

    def delete_node(self, node):
        if isinstance(node, DotNode):
            name = node.name
        else:
            name = node
        try:
            del self._nodes[name]
            del self._allnodes[name]
        except:
            raise DotParsingException("Node %s does not exists" % name)

    def get_node(self, nodename):
        """Return node with name=nodename

        Returns None if node does not exists.
        """
        return self._allnodes.get(nodename, None)

    ##    def add_edge(self, src, dst,attributes={},**kwds):
    ##        src = self.add_node(src)
    ##        dst = self.add_node(dst)
    ##        edge = DotEdge(src,dst,self,attributes=attributes,**kwds)
    ##        edgekey = (src.name,dst.name)
    ##        # Need to set correct edge parent
    ##        # The edge should belong to the node with the lowest level number
    ##        if src.parent.level <= dst.parent.level:
    ##            edge_parent = src.parent
    ##        else:
    ##            edge_parent = dst.parent
    ##
    ##        if edgekey in self._alledges:
    ##            edgs = self._alledges[edgekey]
    ##            if not self.strict:
    ##                edge.parent = edge_parent
    ##                edgs.append(edge)
    ##
    ##            else:
    ##                edgs[0].attributes.update(edge.attributes)
    ##                return edgs[0]
    ##        else:
    ##            edge.parent = edge_parent
    ##
    ##            self._alledges[edgekey] = [edge]
    ##            self._edges[edgekey] = [edge]
    ##
    ##
    ##
    ##
    ##        return edge

    def add_subgraph(self, subgraph, **kwds):
        if isinstance(subgraph, DotSubGraph):
            subgraphcls = subgraph
        else:
            subgraphcls = DotSubGraph(subgraph, self.strict, self.directed, **kwds)
        subgraphcls._allnodes = self._allnodes
        subgraphcls._alledges = self._alledges
        subgraphcls._allgraphs = self._allgraphs
        subgraphcls.parent = self
        subgraphcls.root = self.root
        subgraphcls.level = self.level + 1
        subgraphcls.add_default_node_attr(**self.default_node_attr)
        subgraphcls.add_default_edge_attr(**self.default_edge_attr)
        subgraphcls.add_default_graph_attr(**self.attr)
        subgraphcls.attr.update(self.default_graph_attr)
        subgraphcls.padding += self.padding
        self.subgraphs.append(subgraphcls)
        self._allgraphs.append(subgraphcls)
        return subgraphcls

    def get_subgraphs(self):
        return self.subgraphs

    def get_edges(self):
        return self._edges

    def get_all_nodes(self):
        nodes = []
        for subgraph in self.get_subgraphs():
            nodes.extend(subgraph.get_all_nodes())
        nodes.extend(self._nodes)
        return nodes

    def set_attr(self, **kwds):
        """Set graph attributes"""
        self.attr.update(kwds)
        # self.set_default_graph_attr(kwds)

    nodes = property(lambda self: self._nodes.values())
    allnodes = property(lambda self: self._allnodes.values())
    allgraphs = property(lambda self: self._allgraphs.__iter__())
    alledges = property(lambda self: flatten(self._alledges.values()))
    edges = property(get_edges)

    def __str__(self):
        s = ""
        padding = self.padding
        if len(self.allitems) > 0:
            grstr = "".join(["%s%s" % (padding, n) for n in map(str, flatten(self.allitems))])
            attrstr = ",".join(["%s=%s" % \
                                (quote_if_necessary(key), quote_if_necessary(val)) \
                                for key, val in self.attr.items()])
            if attrstr:
                attrstr = "%sgraph [%s];" % (padding, attrstr)
            if not isinstance(self, DotSubGraph):
                s = ""
                if self.strict:
                    s += 'strict '
                if self.directed:
                    s += "digraph"
                else:
                    s += "graph"
                return "%s %s{\n%s\n%s\n}" % (s, self.get_name(), grstr, attrstr)
            else:
                return "%s %s{\n%s\n%s\n%s}" % ('subgraph', self.get_name(), grstr, attrstr, padding)

        subgraphstr = "\n".join(["%s%s" % (padding, n) for n in map(str, self.subgraphs)])

        nodestr = "".join(["%s%s" % (padding, n)
                           for n in map(str, self._nodes.values())])
        edgestr = "".join(["%s%s" % (padding, n)
                           for n in map(str, flatten(self.edges.values()))])

        attrstr = ",".join(["%s=%s" %
                            (quote_if_necessary(key), quote_if_necessary(val))
                            for key, val in self.attr.items()])
        if attrstr:
            attrstr = "%sgraph [%s];" % (padding, attrstr)
        if not isinstance(self, DotSubGraph):
            s = ""
            if self.strict:
                s += 'strict '
            if self.directed:
                s += "digraph"
            else:
                s += "graph"
            return "%s %s{\n%s\n%s\n%s\n%s\n}" % (s, self.get_name(), subgraphstr, attrstr, nodestr, edgestr)
        else:
            return "%s %s{\n%s\n%s\n%s\n%s\n%s}" % (
                'subgraph', self.get_name(), subgraphstr, attrstr, nodestr, edgestr, padding)


class DotEdge(object):
    """Class representing a DOT edge"""

    def __init__(self, src, dst, directed=False, src_port="", dst_port="", **kwds):
        self.src = src
        self.dst = dst
        self.src_port = src_port
        self.dst_port = dst_port
        # self.parent = parent_graph
        self.attr = {}
        if directed:
            self.conn = "->"
        else:
            self.conn = "--"

        self.attr.update(kwds)

    def __str__(self):
        attrstr = ",".join(["%s=%s" % \
                            (quote_if_necessary(key), quote_if_necessary(val)) \
                            for key, val in self.attr.items()])
        if attrstr:
            attrstr = "[%s]" % attrstr
        return "%s%s %s %s%s %s;\n" % (quote_if_necessary(self.src.name), \
                                       self.src_port, self.conn, \
                                       quote_if_necessary(self.dst.name), self.dst_port, attrstr)

    def get_source(self):
        return self.src.name

    def get_destination(self):
        return self.dst.name

    def __getattr__(self, name):
        try:
            return self.attr[name]
        except KeyError:
            raise AttributeError


class DotSubGraph(DotGraph):
    """Class representing a DOT subgraph"""

    def __init__(self, name='subgG', strict=True, directed=False, **kwds):
        DotGraph.__init__(self, name, strict, directed, **kwds)


testgraph = r"""
/* Test that the various id types are parsed correctly */
digraph G {
    TPR [label=TехПрог];
    LFP [label=ЛиФП, pos="420,686", width="0.86"];
    //"aa\\" -> b [label="12"];
}
"""

if __name__ == '__main__':
    import pprint

    print("Creating parser")
    gp = DotDataParser()
    tok = gp.parse_dot_data_debug(testgraph)
    # dg = parse_dot_data(testgraph)
    pprint.pprint(tok)
    print(gp.graph)
