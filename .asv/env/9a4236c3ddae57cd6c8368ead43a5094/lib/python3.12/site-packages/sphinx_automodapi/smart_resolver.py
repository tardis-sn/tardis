# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The classes in the astropy docs are documented by their API location,
which is not necessarily where they are defined in the source.  This
causes a problem when certain automated features of the doc build,
such as the inheritance diagrams or the `Bases` list of a class
reference a class by its canonical location rather than its "user"
location.

In the `autodoc-process-docstring` event, a mapping from the actual
name to the API name is maintained.  Later, in the `missing-reference`
event, unresolved references are looked up in this dictionary and
corrected if possible.
"""

from docutils.nodes import literal, reference


def process_docstring(app, what, name, obj, options, lines):
    if isinstance(obj, type):
        env = app.env
        if not hasattr(env, 'class_name_mapping'):
            env.class_name_mapping = {}
        mapping = env.class_name_mapping
        mapping[obj.__module__ + '.' + obj.__name__] = name


def merge_mapping(app, env, docnames, env_other):
    if not hasattr(env_other, 'class_name_mapping'):
        return
    if not hasattr(env, 'class_name_mapping'):
        env.class_name_mapping = {}
    env.class_name_mapping.update(env_other.class_name_mapping)


def missing_reference_handler(app, env, node, contnode):
    """
    Handler to be connect to the sphinx 'missing-reference' event.  The handler a
    resolves reference (node) and returns a new node when sphinx could not
    originally resolve the reference.

    see `missing-reference in sphinx documentation
    <https://www.sphinx-doc.org/en/master/extdev/appapi.html#event-missing-reference>`_

    :param app: The Sphinx application object
    :param env: The build environment (``app.builder.env`)
    :param node: The ``pending_xref`` node to be resolved. Its attributes reftype,
                 reftarget, modname and classname attributes determine the type and
                 target of the reference.
    :param contnode: The node that carries the text and formatting inside the
                     future reference and should be a child of the returned
                     reference node.
    """
    # a good example of how a missing reference handle works look to
    #  https://github.com/sphinx-doc/sphinx/issues/1572#issuecomment-68590981
    #
    # Important attributes of the "node":
    #
    #      example role:  :ref:`title <target>`
    #
    #  'reftype'     - role name (in the example above 'ref' is the reftype)
    #  'reftarget'   - target of the role, as given in the role content
    #                  (in the example 'target' is the reftarget
    #  'refexplicit' - the explicit title of the role
    #                  (in the example 'title' is the refexplicit)
    #  'refdoc'      - document in which the role appeared
    #  'refdomain'   - domain of the role, in our case emtpy

    if not hasattr(env, 'class_name_mapping'):
        env.class_name_mapping = {}
    mapping = env.class_name_mapping

    reftype = node['reftype']
    reftarget = node['reftarget']
    refexplicit = node.get('refexplicit')  # default: None
    refdoc = node.get('refdoc', env.docname)
    if reftype in ('obj', 'class', 'exc', 'meth'):
        suffix = ''
        if reftarget not in mapping:
            if '.' in reftarget:
                front, suffix = reftarget.rsplit('.', 1)
            else:
                front = None
                suffix = reftarget

            if suffix.startswith('_') and not suffix.startswith('__'):
                # If this is a reference to a hidden class or method,
                # we can't link to it, but we don't want to have a
                # nitpick warning.
                return node[0].deepcopy()

            if reftype in ('obj', 'meth') and front is not None:
                if front in mapping:
                    reftarget = front
                    suffix = '.' + suffix

            if (reftype in ('class', ) and '.' in reftarget and
                    reftarget not in mapping):

                if '.' in front:
                    reftarget, _ = front.rsplit('.', 1)
                    suffix = '.' + suffix
                reftarget = reftarget + suffix
                prefix = reftarget.rsplit('.')[0]
                inventory = getattr(env, 'intersphinx_named_inventory', {})
                if (reftarget not in mapping and
                        prefix in inventory):

                    if 'py:class' in inventory[prefix] and \
                            reftarget in inventory[prefix]['py:class']:
                        newtarget = inventory[prefix]['py:class'][reftarget][2]
                        if not refexplicit and '~' not in node.rawsource:
                            contnode = literal(text=reftarget)
                        newnode = reference('', '', internal=True)
                        newnode['reftitle'] = reftarget
                        newnode['refuri'] = newtarget
                        newnode.append(contnode)

                        return newnode

        if reftarget in mapping:
            newtarget = mapping[reftarget] + suffix
            if not refexplicit and '~' not in node.rawsource:
                contnode = literal(text=newtarget)
            newnode = env.domains['py'].resolve_xref(env, refdoc, app.builder, 'class',
                                                     newtarget, node, contnode)
            if newnode is not None:
                newnode['reftitle'] = reftarget
            return newnode


def setup(app):

    app.connect('autodoc-process-docstring', process_docstring)
    app.connect('missing-reference', missing_reference_handler)
    app.connect('env-merge-info', merge_mapping)

    return {'parallel_read_safe': True,
            'parallel_write_safe': True}
