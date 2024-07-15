# -*- coding: utf-8 -*-
"""
    WideFormat layout engine
    ------------------------

    In this layout for each nesting level the table is extended by
    one or more columns.

    :copyright: Copyright 2017-2021, Leo Noordergraaf
    :licence: GPL v3, see LICENCE for details.
"""

from sys import version_info
from copy import deepcopy
from pathlib import Path
from docutils import statemachine
from docutils import nodes
from docutils.nodes import fully_normalize_name as normalize_name

if version_info[0] == 2:
    str_unicode = unicode
else:
    str_unicode = str

NOESC = ':noesc:'  # prefix marker to indicate string must not be escaped.

class WideFormat(object):
    KV_SIMPLE = [
        'multipleOf', 'maximum', 'exclusiveMaximum', 'minimum',
        'exclusiveMinimum', 'maxLength', 'minLength', 'pattern', 'default',
        'format', 'const'
    ]

    KV_ARRAY = ['maxItems', 'minItems', 'uniqueItems']

    KV_OBJECT = ['maxProperties', 'minProperties']

    COMBINATORS = ['allOf', 'anyOf', 'oneOf']

    SINGLEOBJECTS = ['not']

    CONDITIONAL = ["if", "then", "else"]

    option_defaults = {
        'lift_title': True,
        'lift_description': False,
        'lift_definitions': False,
        'auto_target': False,
        'auto_reference': False
    }

    def __init__(self, state, lineno, source, options, app):
        super(WideFormat, self).__init__()
        self.app = app
        self.trans = None
        self.lineno = lineno
        self.state = state
        self.filename = self._get_filename(source)
        self.nesting = 0
        self.ref_titles = {}
        self.target_pointer = '#'

        self.options = deepcopy(self.option_defaults)
        self.options.update(app.config.jsonschema_options)
        self.options.update(options)

    def run(self, schema, pointer=''):
        # To set the correct auto target for a nested definitions we need to save
        # the current pointer that may be used inside recursive run append on
        before = self.target_pointer
        self.target_pointer += pointer

        result = []
        target = self._target(schema, self.target_pointer)
        section = self._section(schema)

        table, definitions = self.transform(schema)

        # restore to previous pointer
        self.target_pointer = before

        if target:
            result.append(target)

        if section:
            section += table
            section += definitions
            result.append(section)
        else:
            if table:
                result.append(table)
            if definitions:
                result.extend(definitions)

        return result

    def transform(self, schema):
        body, definitions = self._dispatch(schema)
        if len(body) > 0:
            cols, head, body = self._cover(schema, body)
            table = self.state.build_table((cols, head, body), self.lineno)
        else:
            table = None
        return table, definitions

    def _target(self, schema, pointer=''):
        # Wrap section and table in a target (anchor) node so
        # that it can be referenced from other sections.
        labels = self.app.env.domaindata['std']['labels']
        anonlabels = self.app.env.domaindata['std']['anonlabels']
        docname = self.app.env.docname

        targets = []

        if '$$target' in schema:
            if not isinstance(schema['$$target'], list):
                targets.append(schema['$$target'])
            else:
                targets.extend(schema['$$target'])
            del schema['$$target']

        if self.options['auto_target']:
            # When schema's multiple schema's are writen with content but without a pointer
            # you get multiple equal named targets, all $ref will link to the last created schema
            # The same applies if you would load files with equal name into your documentation
            if len(pointer) > 1:
                targets.append(self.filename + pointer)
            else:
                targets.append(self.filename)

        if targets:
            targetnode = nodes.target()
            anchorid = nodes.make_id((schema['title'] if 'title' in schema else targets[0]))
            targetnode['ids'].append(anchorid)
            targetnode['names'].append(anchorid)
            targetnode.line = self.lineno

            for target in targets:
                anchor = normalize_name(target)
                anonlabels[anchor] = docname, targetnode['ids'][0]
                labels[anchor] = docname, targetnode['ids'][0], (schema['title'] if 'title' in schema else target)

            return targetnode

        return None

    def _section(self, schema):
        if 'title' in schema and self.options['lift_title']:
            # Wrap the resulting table in a section giving it a caption and an
            # entry in the table of contents.
            # unable to use self.state.section() to make a section as style is unknown
            # all sections will be placed inside current section
            section_node = nodes.section()
            textnodes, title_messages = self.state.inline_text(schema['title'], self.lineno)
            titlenode = nodes.title(schema['title'], '', *textnodes)
            name = normalize_name(titlenode.astext())
            section_node['names'].append(name)
            section_node += titlenode
            section_node += title_messages
            self.state.document.note_implicit_target(section_node, section_node)

            if self.nesting == 0:
                self.ref_titles[self.nesting] = schema['title']

            if self.options['lift_description']:
                self._get_description(schema, section_node)

            del schema['title']
            return section_node

        return None

    def _get_description(self, schema, node):
        if 'description' in schema:
            self.state.nested_parse(self._convert_content(schema['description']), self.lineno, node)
            del schema['description']

        if '$$description' in schema:
            if isinstance(schema['$$description'], list):
                schema['$$description'] = '\n'.join(schema['$$description'])
            self.state.nested_parse(self._convert_content(schema['$$description']), self.lineno, node)
            del schema['$$description']

    def _cover(self, schema, body):
        # Patch up and finish the table.
        head = []

        # Outermost id becomes schema url
        # NB: disregards interior id's
        # to support both 'id' draft 4 only and '$id' from draft 6
        if 'id' in schema:
            body.insert(0, self._line(self._cell(schema['id'])))
            del schema['id']
        elif '$id' in schema:
            body.insert(0, self._line(self._cell(schema['$id'])))
            del schema['$id']

        # patch up if necessary, all rows should be of equal length
        nrcols = self._square(body)
        # assume len(head[n]) <= nrcols
        nrcols = self._square(head, nrcols)

        # create column spans and proper type casts
        self._calc_spans(head, nrcols)
        self._calc_spans(body, nrcols)

        # All columns have same width, to change alter the first element
        return [1] * nrcols, head, body

    def _dispatch(self, schema, label=None):
        # Main driver of the recursive schema traversal.
        rows = []
        self.nesting += 1

        definitions = []
        if self.options['lift_definitions']:
            if '$defs' in schema:
                definitions = self._definitions(schema, '$defs')
            elif 'definitions' in schema:
                definitions = self._definitions(schema, 'definitions')

        if 'type' in schema:
            # select processor for type
            if 'object' in schema['type']:
                rows = self._objecttype(schema)
            elif 'array' in schema['type']:
                rows = self._arraytype(schema)
        else:
            rows = self._objecttype(schema)
            self._check_description(schema, rows)
        rows.extend(self._simpletype(schema))

        if '$ref' in schema:
            rows.extend(self._reference(schema))

        rows.extend(self._complexstructures(schema))

        # definitions aren't really type equiv's but still best place for them
        rows.extend(self._objectproperties(schema, 'definitions'))
        rows.extend(self._objectproperties(schema, '$defs'))

        if label is not None:
            # prepend label column if required
            rows = self._prepend(label, rows)

        self.nesting -= 1
        return rows, definitions

    def _objecttype(self, schema):
        # create description and type rows
        rows = self._simpletype(schema)
        rows.extend(self._objectproperties(schema, 'properties'))
        rows.extend(self._objectproperties(schema, 'patternProperties'))
        rows.extend(self._bool_or_object(schema, 'additionalProperties'))
        rows.extend(self._dependencies(schema, 'dependencies'))
        rows.extend(self._kvpairs(schema, self.KV_OBJECT))
        return rows

    def _arraytype(self, schema):
        def oneline(label, item):
            if isinstance(item, dict):
                rows.extend(self._dispatch(item, label)[0])
            else:
                rows.append(self._line(label, self._cell(item)))

        # create description and type rows
        rows = self._simpletype(schema)

        if 'items' in schema:
            if type(schema['items']) == list:
                rows.append(self._line(self._cell('items')))
                for item in schema['items']:
                    label = self._cell('-')
                    oneline(label, item)
            else:
                oneline(self._cell('items'), schema['items'])
            del schema['items']

        rows.extend(self._bool_or_object(schema, 'additionalItems'))
        rows.extend(self._kvpairs(schema, self.KV_ARRAY))
        return rows

    def _simpletype(self, schema):
        rows = []

        if 'title' in schema and (not self.options['lift_title'] or self.nesting > 1):
            rows.append(self._line(self._cell('*' + schema['title'] + '*')))
            del schema['title']

        self._check_description(schema, rows)

        if 'type' in schema:
            rows.append(
                self._line(self._cell('type'),
                           self._decodetype(schema['type'])))
            del schema['type']

        if 'enum' in schema:
            rows.append(
                self._line(
                    self._cell('enum'),
                    self._cell(', '.join(
                        [str_unicode(e) for e in schema['enum']]))))
            del schema['enum']

        if 'examples' in schema:
            rows.extend(self._examples(schema['examples']))
            del schema['examples']

        rows.extend(self._kvpairs(schema, self.KV_SIMPLE))
        return rows

    def _objectproperties(self, schema, key):
        # process the `properties` key of the object type
        # used for `properties`, `patternProperties` and
        # `definitions`.
        rows = []

        if key in schema:
            rows.append(self._line(self._cell(key)))

            for prop in schema[key].keys():
                # insert spaces around the regexp OR operator
                # allowing the regexp to be split over multiple lines.
                proplist = prop.split('|')
                dispprop = self._escape(' | '.join(proplist))
                bold = ''
                if 'required' in schema:
                    if prop in schema['required']:
                        bold = '**'
                label = self._cell('- ' + bold + dispprop + bold)

                if isinstance(schema[key][prop], dict):
                    obj = schema[key][prop]
                    rows.extend(self._dispatch(obj, label)[0])
                else:
                    rows.append(self._line(label, self._cell(schema[key][prop])))
            del schema[key]
        return rows

    def _definitions(self, schema, defs_key):
        target = {}
        for name, item in schema[defs_key].items():
            # add title by the name of the object if title not defined
            if 'title' not in item:
                item['title'] = name

            target[name] = item['title']

        # To automate $ref to definitions titles save a copy
        # of the schema before continuing the recursive build
        # so reference can be set with the correct title
        if self.nesting in self.ref_titles:
            self.ref_titles[self.nesting].update(target)
        else:
            self.ref_titles[self.nesting] = target


        result = []
        for name, item in schema[defs_key].items():
            new_target = '/{defs_key}/{name}'.format(
                defs_key=defs_key, name=name
            )
            result.extend(self.run(item, new_target))

        del schema[defs_key]
        return result

    def _complexstructures(self, schema):
        rows = []

        for k in self.COMBINATORS:
            # combinators belong at this level as alternative to type
            if k in schema:
                items = []
                for s in schema[k]:
                    content = self._dispatch(s)[0]
                    if content:
                        items.extend(content)
                if items:
                    rows.extend(self._prepend(self._cell(k), items))
                del schema[k]

        for k in self.SINGLEOBJECTS:
            # combinators belong at this level as alternative to type
            if k in schema:
                rows.extend(self._dispatch(schema[k], self._cell(k))[0])
                del schema[k]

        if self.CONDITIONAL[0] in schema:
            # only if 'if' in schema there would be a needs to go through if, then & else
            items = []
            for k in self.CONDITIONAL:
                if k in schema:
                    content = self._dispatch(schema[k])[0]
                    if content:
                        items.append(self._prepend(self._cell(k), content))
                    del schema[k]
            if len(items) >= 2:
                for item in items:
                    rows.extend(item)

        return rows

    def _dependencies(self, schema, key):
        rows = []

        if key in schema:
            rows.append(self._line(self._cell(key)))

            for prop in schema[key].keys():
                label = self._cell('- ' + prop)
                obj = schema[key][prop]
                if type(obj) == list:
                    rows.append(
                        self._line(
                            label,
                            self._cell(str_unicode(', '.join(obj)))))
                else:
                    rows.extend(self._dispatch(obj, label)[0])
            del schema[key]
        return rows

    def _get_defined_reference(self, schema, key):
        target = '#/{defs}/'.format(defs=key)
        if schema['$ref'].startswith(target):
            reference = [r for r in schema['$ref'][2:].split('/') if r != key]
            return len(reference), reference[-1]

    def _reference(self, schema):
        if self.options['auto_reference'] and self.options['lift_title']:
            # first check if references is to own schema
            # when definitions is separated automated they will be linked to the title
            # otherwise it will only be a string
            reference = (
                self._get_defined_reference(schema, 'definitions') or
                self._get_defined_reference(schema, '$defs')
            )
            if schema['$ref'] == '#' or schema['$ref'] == '#/':
                if self.ref_titles.get(0, False):
                    row = (self._line(self._cell('`' + self.ref_titles[0] + '`_')))
                else:
                    row = (self._line(self._cell(schema['$ref'])))
            elif reference:
                ref_length, target_name = reference
                # Check if there are definitions available to make a reference
                if (self.ref_titles.get(ref_length, False) and
                        target_name in self.ref_titles[ref_length]):
                    ref_title = self.ref_titles[ref_length][target_name]
                    row = (self._line(self._cell('`' + ref_title + '`_')))
                else:
                    row = (self._line(self._cell(schema['$ref'])))
            elif schema['$ref'].startswith("#/"):
                # Other references to own schema should not be defined as :ref: only as string
                row = (self._line(self._cell(schema['$ref'])))
            elif schema['$ref'].startswith("http"):
                row = (self._line(self._cell(schema['$ref'])))
            elif "#/" in schema['$ref']:
                row = (self._line(self._cell(':ref:`' + self._get_filename(schema['$ref'], True) + '`')))
            else:
                row = (self._line(self._cell(':ref:`' + self._get_filename(schema['$ref']) + '`')))
        elif self.options['auto_reference'] and not self.options['lift_title']:
            # when using reference without titles we need to reference to our own targets
            # if auto_target is False linking won't work
            row = (self._line(self._cell(':ref:`' + self.filename + schema['$ref'] + '`')))
        else:
            row = (self._line(self._cell(':ref:`' + schema['$ref'] + '`')))
        del schema['$ref']
        return [row]

    def _bool_or_object(self, schema, key):
        # for those attributes that accept either a boolean or a schema.
        rows = []

        if key in schema:
            if type(schema[key]) == bool:
                rows.append(self._line(self._cell(key), self._cell(schema[key])))
                del schema[key]
            else:
                rows.extend(self._dispatch(schema[key], self._cell(key))[0])
                del schema[key]

        return rows

    def _kvpairs(self, schema, keys):
        # render key-value pairs
        rows = []

        for k in keys:
            if k in schema:
                value = schema[k]
                if k == 'pattern':
                    value = self._escape(value)
                if k == 'default':
                    rows.extend(self._prepend(self._cell(k), self._render_any_value(value)))
                else:
                    rows.append(self._line(self._cell(k), self._cell(value)))
                del schema[k]
        return rows

    def _prepend(self, prepend, rows):
        # prepend a label to a set of rows
        rcnt = len(rows)

        if rcnt == 0:
            # return a row with only the label
            return [self._line(prepend)]
        else:
            # add the label to the first row
            prepend[0] = rcnt - 1
            rows[0].insert(0, prepend)
            # following rows have an empty column prepended
            for r in range(1, rcnt):
                rows[r].insert(0, None)
            return rows

    def _decodetype(self, typ):
        # render (array of) simple type(s)
        if type(typ) == list:
            # construct list of basic types
            return self._cell(' / '.join(['*' + s + '*' for s in typ]))
        else:
            # create a single type
            return self._cell('*' + typ + '*')

    def _examples(self, examples):
        # Render examples as rows
        rows = []
        rows.extend(self._render_any_value(examples))
        rows = self._prepend(self._cell('examples'), rows)
        return rows

    def _check_description(self, schema, rows):
        if 'description' in schema:
            rows.append(self._line(self._cell(schema['description'])))
            del schema['description']
        if '$$description' in schema:
            if isinstance(schema['$$description'], list):
                rows.append(self._line(self._cell('\n'.join(schema['$$description']))))
                del schema['$$description']
            else:
                rows.append(self._line(self._cell(schema['$$description'])))
                del schema['$$description']

    def _render_any_value(self, value):
        # render a single value, an array of values or a dict with key/value pairs
        rows = []
        if isinstance(value, list):
            if len(value) == 0:
                rows.append(self._line(self._cell('')))
            else:
                for v in value:
                    rows.extend(self._render_any_value(v))
        elif isinstance(value, dict):
            if len(value) == 0:
                rows.append(self._line(self._cell('')))
            else:
                for k in value:
                    rows.extend(self._prepend(self._cell(k), self._render_any_value(value[k])))
        elif isinstance(value, str):
            rows.append(self._line(self._cell(self._escape(value) if value is not None else "null")))
        else:
            rows.append(self._line(self._cell(value if value is not None else "null")))
        return rows

    def _square(self, rows, nrcols=0):
        # determine max. number of columns
        if nrcols == 0:
            for row in rows:
                nrcols = max(nrcols, len(row))

        # extend each row to contain same number of columns
        for row in rows:
            if len(row) < nrcols:
                row += [None] * (nrcols - len(row))

        return nrcols

    def _calc_spans(self, rows, nrcols):
        # calculate colspan
        for row in rows:
            target = None
            for c in range(nrcols):
                if row[c] is not None:
                    # try to extend colspan on this cell
                    target = row[c]
                else:
                    if target is not None:
                        # extend colspan
                        target[1] += 1

        # convert arrays to tuples
        # arrays are needed to patch up colspan and rowspan
        # the table builder requires each cell to be a tuple, not an array
        for row in rows:
            for c in range(nrcols):
                if row[c] is not None:
                    row[c] = tuple(row[c])

    def _line(self, *cells):
        # turn a number of cells into a list
        return [c for c in cells]

    def _cell(self, text):
        # Table builder wants all cells as a tuple of 4 fields.
        # Returns a list since it needs to be mutable (tuple isn't).
        # Eventually, _calc_spans() will turn these lists into tuples.
        return [
            0,  # rowspan
            0,  # colspan
            self.lineno,  # source line number

            # turn string into multiline array of views on lists
            # required by table builder
            self._convert_content(text)
        ]

    def _convert_content(self, text):
        list_lines = statemachine.string2lines(str_unicode(text))
        # Adding a source and line number to each line text warnings may appear when writing if there are issues with a line
        # if left None the warnings would be counted but don't appear in the output then you don't know the source of it
        items = [(self.state.document.current_source, self.lineno)] * len(list_lines)

        return statemachine.StringList(list_lines, items=items)

    def _get_filename(self, path, include_pointer=False):
        # gets from filepath or url the name of file
        if '#/' in path:
            path, pointer = path.rsplit('#/', 1)
        if include_pointer:
            return Path(path).name + '#/' + pointer
        return Path(path).name

    def _escape(self, text):
        if text.startswith(NOESC):
            return text[len(NOESC):]
        text = text.replace('\\', '\\\\')
        text = text.replace('_', '\\_')
        text = text.replace('*', '\\*')
        return text
