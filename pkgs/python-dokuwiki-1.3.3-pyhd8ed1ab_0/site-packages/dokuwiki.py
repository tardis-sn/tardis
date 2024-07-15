# -*- coding: utf-8 -*-

"""This python module aims to manage
`DokuWiki <https://www.dokuwiki.org/dokuwiki>`_ wikis by using the
provided `XML-RPC API <https://www.dokuwiki.org/devel:xmlrpc>`_.  It is
compatible with python2.7 and python3+.

Installation
------------
It is on `PyPi <https://pypi.python.org/pypi/dokuwiki>`_ so you can use
the ``pip`` command to install it::

    pip install dokuwiki

Otherwise sources are in `github <https://github.com/fmenabe/python-dokuwiki>`_
"""

import re
import sys
import base64
import weakref
from xml.parsers.expat import ExpatError

PY_VERSION = sys.version_info[0]
if PY_VERSION == 3:
    from xmlrpc.client import ServerProxy, Binary, Fault, Transport, SafeTransport, ProtocolError
    from urllib.parse import quote
else:
    from xmlrpclib import ServerProxy, Binary, Fault, Transport, SafeTransport, ProtocolError
    from urllib import quote

from datetime import datetime, timedelta

ERR = 'XML or text declaration not at start of entity: line 2, column 0'

_URL_RE = re.compile(r'(?P<proto>https?)://(?P<host>[^/]*)(?P<uri>/.*)?')

def date(date):
    """DokuWiki returns dates of `xmlrpclib`/`xmlrpc.client` ``DateTime``
    type and the format changes between DokuWiki versions ... This function
    convert *date* to a `datetime` object.
    """
    date = date.value
    return (datetime.strptime(date[:-5], '%Y-%m-%dT%H:%M:%S')
            if len(date) == 24
            else datetime.strptime(date, '%Y%m%dT%H:%M:%S'))

def utc2local(date):
    """DokuWiki returns date with a +0000 timezone. This function convert *date*
    to the local time.
    """
    date_offset = (datetime.now() - datetime.utcnow())
    # Python < 2.7 don't have the 'total_seconds' method so calculate it by hand!
    date_offset = (date_offset.microseconds +
                   (date_offset.seconds + date_offset.days * 24 * 3600) * 1e6) / 1e6
    date_offset = int(round(date_offset / 60 / 60))
    return date + timedelta(hours=date_offset)


class DokuWikiError(Exception):
    """Exception raised by this module when there is an error."""
    pass


def CookiesTransport(proto='https'):
    """Generate transport class when using cookie based authentication."""
    _TransportClass_ = Transport if proto == 'http' else SafeTransport

    class CookiesTransport(_TransportClass_):
        """A Python3 xmlrpc.client.Transport subclass that retains cookies."""
        def __init__(self):
            _TransportClass_.__init__(self)
            self._cookies = dict()

        def send_headers(self, connection, headers):
            if self._cookies:
                cookies = map(lambda x: x[0] + '=' + x[1], self._cookies.items())
                connection.putheader('Cookie', '; '.join(cookies))
            _TransportClass_.send_headers(self, connection, headers)

        def parse_response(self, response):
            """parse and store cookie"""
            try:
                for header in response.msg.get_all("Set-Cookie"):
                    cookie = header.split(";", 1)[0]
                    cookieKey, cookieValue = cookie.split("=", 1)
                    self._cookies[cookieKey] = cookieValue
            finally:
                return _TransportClass_.parse_response(self, response)

    class CookiesTransport2(_TransportClass_):
        """A Python2 xmlrpclib.Transport subclass that retains cookies."""
        def __init__(self):
            _TransportClass_.__init__(self)
            self._cookies = dict()

        def send_request(self, connection, handler, request_body):
            _TransportClass_.send_request(self, connection, handler, request_body)
            # set cookie below handler
            if self._cookies:
                cookies = map(lambda x: x[0] + '=' + x[1], self._cookies.items())
                connection.putheader("Cookie", "; ".join(cookies))

        def parse_response(self, response):
            """parse and store cookie"""
            try:
                for header in response.getheader("set-cookie").split(", "):
                    # filter 'expire' information
                    if not header.startswith("D"):
                        continue
                    cookie = header.split(";", 1)[0]
                    cookieKey, cookieValue = cookie.split("=", 1)
                    self._cookies[cookieKey] = cookieValue
            finally:
                return _TransportClass_.parse_response(self, response)

    return CookiesTransport2() if PY_VERSION == 2 else CookiesTransport()

class DokuWiki(object):
    """Initialize a connection to a DokuWiki wiki. ``url``, ``user`` and
    ``password`` are respectively the URL, the login and the password for
    connecting to the wiki. ``kwargs`` are `xmlrpclib`/`xmlrpc.client`
    **ServerProxy** parameters.

    The exception `DokuWikiError` is raised if the authentication
    fails but others exceptions (like `socket.gaierror` for invalid domain,
    `xmlrpc.client.ProtocolError` for an invalid wiki, ...) are not catched.

    .. code::

        try:
            wiki = dokuwiki.DokuWiki('URL', 'USER', 'PASSWORD')
        except (DokuWikiError, Exception) as err:
            print('unable to connect: %s' % err)

    .. note::

        The URL format is: ``PROTO://FQDN[/PATH]`` (*https://www.example.com/dokuwiki*
        for example).

    To use cookie based authentication (use HTTP cookies instead of passing login
    and password in the URI), set ``cookieAuth`` parameter to *True*:

    .. code::

        wiki = dokuwiki.DokuWiki('URL', 'USER', 'PASSWORD', cookieAuth=True)
    """
    def __init__(self, url, user, password, **kwargs):
        """Initialize the object by connecting to the XMLRPC server."""
        # Parse input URL
        try:
            params = _URL_RE.search(url).groupdict()
        except AttributeError:
            raise DokuWikiError("invalid url '%s'" %  url)

        # Set auth string or transport for cookie based authentication.
        auth = '{:s}:{:s}@'.format(user, quote(password, safe=''))
        cookie_auth = kwargs.pop('cookieAuth', False)
        if cookie_auth:
            auth = ''
            kwargs['transport'] = CookiesTransport(params['proto'])

        xmlrpc_url = '%s://%s%s%s/lib/exe/xmlrpc.php' % (
            params['proto'], auth, params['host'], params['uri'] or '')
        self.proxy = ServerProxy(xmlrpc_url, **kwargs)

        # Force login for cookie based authentication.
        if cookie_auth and not self.login(user, password):
            raise DokuWikiError('invalid login or password!')

        # Dummy call to ensure the connection is up.
        try:
            self.version
        except ProtocolError as err:
            if err.errcode == 401:
                raise DokuWikiError('invalid login or password!')
            raise

        # Set "namespaces" for pages and medias functions.
        self.pages = _Pages(weakref.ref(self)())
        self.medias = _Medias(weakref.ref(self)())
        self.structs = _Structs(weakref.ref(self)())

    def send(self, command, *args, **kwargs):
        """Generic method for executing an XML-RPC *command*. *args* and
        *kwargs* are the arguments and parameters needed by the command.
        """
        args = list(args)
        if kwargs:
            args.append(kwargs)

        method = self.proxy
        for elt in command.split('.'):
            method = getattr(method, elt)

        try:
            return method(*args)
        except Fault as err:
            if err.faultCode == 121:
                return {}
            elif err.faultCode == 321:
                return []
            raise DokuWikiError(err)
        except ExpatError as err:
            if str(err) != ERR:
                raise DokuWikiError(err)

    @property
    def version(self):
        """Property that returns the DokuWiki version of the remote Wiki."""
        return self.send('dokuwiki.getVersion')

    @property
    def time(self):
        """Property that returns the current time at the remote wiki server as
        Unix timestamp.
        """
        return self.send('dokuwiki.getTime')

    @property
    def xmlrpc_version(self):
        """Property that returns the XML RPC interface version of the remote
        Wiki. This is DokuWiki implementation specific and independent of the
        supported standard API version returned by ``wiki.getRPCVersionSupported``.
        """
        return self.send('dokuwiki.getXMLRPCAPIVersion')

    @property
    def xmlrpc_supported_version(self):
        """Property that returns *2* with the supported RPC API version."""
        return self.send('wiki.getRPCVersionSupported')

    @property
    def title(self):
        """Property that returns the title of the wiki."""
        return self.send('dokuwiki.getTitle')

    def login(self, user, password):
        """Log to the wiki using *user* and *password* credentials. It returns
        a boolean that indicates if the user succesfully authenticate."""
        return self.send('dokuwiki.login', user, password)

    def add_acl(self, scope, user, permission):
        """Add an `ACL <https://www.dokuwiki.org/acl>`_ rule that restricts
        the page/namespace *scope* to *user* (use *@group* syntax for groups)
        with *permission* level. It returns a boolean that indicate if the rule
        was correctly added.
        """
        return self.send('plugin.acl.addAcl', scope, user, permission)

    def del_acl(self, scope, user):
        """Delete any ACL matching the given *scope* and *user* (or group if
        *@group* syntax is used). It returns a boolean that indicate if the rule
        was correctly removed.
        """
        return self.send('plugin.acl.delAcl', scope, user)


class _Pages(object):
    """This object regroup methods for managing pages of a DokuWiki. This object
    is accessible from the ``pages`` property of an `DokuWiki` instance::

        wiki = dokuwiki.DokuWiki('URL', 'User', 'Password')
        wiki.pages.list()
    """

    def __init__(self, dokuwiki):
        self._dokuwiki = dokuwiki

    def list(self, namespace='/', **options):
        """List all pages of the given *namespace*.

        Valid *options* are:

            * *depth*: (int) recursion level, 0 for all
            * *hash*: (bool) do an md5 sum of content
            * *skipacl*: (bool) list everything regardless of ACL
        """
        return self._dokuwiki.send('dokuwiki.getPagelist', namespace, options)

    def changes(self, timestamp):
        """Returns a list of changes since given *timestamp*.

        For example, for returning all changes since *2016-01-01*::

            from datetime import datetime
            wiki.pages.changes(datetime(2016, 1, 1).timestamp())
        """
        return self._dokuwiki.send('wiki.getRecentChanges', timestamp)

    def search(self, string):
        """Performs a fulltext search on *string* and returns the first 15
        results.
        """
        return self._dokuwiki.send('dokuwiki.search', string)

    def versions(self, page, offset=0):
        """Returns the available versions of *page*. *offset* can be used to
        list earlier versions in the history.
        """
        return self._dokuwiki.send('wiki.getPageVersions', page, offset)

    def info(self, page, version=None):
        """Returns informations of *page*. Informations of the last version
        is returned if *version* is not set.
        """
        return (self._dokuwiki.send('wiki.getPageInfoVersion', page, version)
                if version is not None
                else self._dokuwiki.send('wiki.getPageInfo', page))

    def get(self, page, version=None):
        """Returns the content of *page*. The content of the last version is
        returned if *version* is not set.
        """
        return (self._dokuwiki.send('wiki.getPageVersion', page, version)
                if version is not None
                else self._dokuwiki.send('wiki.getPage', page))


    def append(self, page, content, **options):
        """Appends *content* text to *page*.

        Valid *options* are:

            * *sum*: (str) change summary
            * *minor*: (bool) whether this is a minor change
        """
        return self._dokuwiki.send('dokuwiki.appendPage', page, content, options)

    def html(self, page, version=None):
        """Returns HTML content of *page*. The HTML content of the last version
        of the page is returned if *version* is not set.
        """
        return (self._dokuwiki.send('wiki.getPageHTMLVersion', page, version)
                if version is not None
                else self._dokuwiki.send('wiki.getPageHTML', page))

    def set(self, page, content, **options):
        """Set/replace the *content* of *page*.

        Valid *options* are:

            * *sum*: (str) change summary
            * *minor*: (bool) whether this is a minor change
        """
        try:
            return self._dokuwiki.send('wiki.putPage', page, content, options)
        except ExpatError as err:
            # Sometime the first line of the XML response is blank which raise
            # the 'ExpatError' exception although the change has been done. This
            # allow to ignore the error.
            if str(err) != ERR:
                raise DokuWikiError(err)

    def delete(self, page):
        """Delete *page* by setting an empty content."""
        return self.set(page, '')

    def lock(self, page):
        """Locks *page*."""
        result = self._dokuwiki.send('dokuwiki.setLocks',
                                     lock=[page], unlock=[])
        if result['lockfail']:
            raise DokuWikiError('unable to lock page')

    def unlock(self, page):
        """Unlocks *page*."""
        result = self._dokuwiki.send('dokuwiki.setLocks',
                                     lock=[], unlock=[page])
        if result['unlockfail']:
            raise DokuWikiError('unable to unlock page')

    def permission(self, page):
        """Returns the permission level of *page*."""
        return self._dokuwiki.send('wiki.aclCheck', page)

    def links(self, page):
        """Returns a list of all links contained in *page*."""
        return self._dokuwiki.send('wiki.listLinks', page)

    def backlinks(self, page):
        """Returns a list of all links referencing *page*."""
        return self._dokuwiki.send('wiki.getBackLinks', page)


class _Medias(object):
    """This object regroup methods for managing medias of a DokuWiki. This
    object is accessible from the ``medias`` property of an `DokuWiki`
    instance::

        wiki = dokuwiki.DokuWiki('URL', 'User', 'Password')
        wiki.medias.list()
    """
    def __init__(self, dokuwiki):
        self._dokuwiki = dokuwiki

    def list(self, namespace='/', **options):
        """Returns all medias of the given *namespace*.

        Valid *options* are:

            * *depth*: (int) recursion level, 0 for all
            * *skipacl*: (bool) skip acl checking
            * *pattern*: (str) check given pattern
            * *hash*: (bool) add hashes to result list
        """
        return self._dokuwiki.send('wiki.getAttachments', namespace, options)

    def changes(self, timestamp):
        """Returns the list of medias changed since given *timestamp*.

        For example, for returning all changes since *2016-01-01*::

            from datetime import datetime
            wiki.medias.changes(datetime(2016, 1, 1).timestamp())
        """
        return self._dokuwiki.send('wiki.getRecentMediaChanges', timestamp)

    def get(self, media, dirpath=None, filename=None, overwrite=False, b64decode=False):
        """Returns the binary data of *media* or save it to a file. If *dirpath*
        is not set the binary data is returned, otherwise the data is saved
        to a file. By default, the filename is the name of the media but it can
        be changed with *filename* parameter. *overwrite* parameter allow to
        overwrite the file if it already exists locally.
        """
        import os
        data = self._dokuwiki.send('wiki.getAttachment', media)
        data = base64.b64decode(data) if b64decode else data.data
        if dirpath is None:
            return data

        if filename is None:
            filename = media.replace('/', ':').split(':')[-1]
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        filepath = os.path.join(dirpath, filename)
        if os.path.exists(filepath) and not overwrite:
            raise FileExistsError("[Errno 17] File exists: '%s'" % filepath)

        with open(filepath, 'wb') as fhandler:
            fhandler.write(data)

    def info(self, media):
        """Returns informations of *media*."""
        return self._dokuwiki.send('wiki.getAttachmentInfo', media)

    def add(self, media, filepath, overwrite=True):
        """Set *media* from local file *filepath*. *overwrite* parameter specify
        if the media must be overwrite if it exists remotely.
        """
        with open(filepath, 'rb') as fhandler:
            self._dokuwiki.send('wiki.putAttachment', media,
                                Binary(fhandler.read()), ow=overwrite)

    def set(self, media, _bytes, overwrite=True, b64encode=False):
        """Set *media* from *_bytes*. *overwrite* parameter specify if the media
        must be overwrite if it exists remotely.
        """
        data = base64.b64encode(_bytes) if b64encode else Binary(_bytes)
        self._dokuwiki.send('wiki.putAttachment', media, data, ow=overwrite)

    def delete(self, media):
        """Delete *media*."""
        return self._dokuwiki.send('wiki.deleteAttachment', media)


class _Structs(object):
    def __init__(self, dokuwiki):
        """Get the structured data of a given page."""
        self._dokuwiki = dokuwiki

    def get_data(self, page, schema='', timestamp=0):
        """Get the structured data of a given page."""
        return self._dokuwiki.send('plugin.struct.getData', page, schema, timestamp)

    def save_data(self, page, data, summary='', minor=False):
        """Saves data for a given page (creates a new revision)."""
        return self._dokuwiki.send('plugin.struct.saveData', page, data, summary, minor)

    def get_schema(self, name=''):
        """Get info about existing schemas columns."""
        return self._dokuwiki.send('plugin.struct.getSchema', name)

    def get_aggregation_data(self, schemas, columns, data_filter=[], sort=''):
        """Get the data that would be shown in an aggregation."""
        return self._dokuwiki.send(
            'plugin.struct.getAggregationData', schemas, columns, data_filter, sort)


class Dataentry(object):
    """Object that manage `data entries <https://www.dokuwiki.org/plugin:data>`_."""

    @staticmethod
    def get(content, keep_order=False):
        """Get dataentry from *content*. *keep_order* indicates whether to
        return an ordered dictionary."""
        if keep_order:
            from collections import OrderedDict
            dataentry = OrderedDict()
        else:
            dataentry = {}

        found = False
        for line in content.split('\n'):
            if line.strip().startswith('---- dataentry'):
                found = True
                continue
            elif line == '----':
                break
            elif not found:
                continue

            line_split = line.split(':')
            key = line_split[0].strip()
            value = re.sub('#.*$', '', ':'.join(line_split[1:])).strip()
            dataentry.setdefault(key, value)

        if not found:
            raise DokuWikiError('no dataentry found')
        return dataentry

    @staticmethod
    def gen(name, data):
        """Generate dataentry *name* from *data*."""
        return '---- dataentry %s ----\n%s\n----' % (name, '\n'.join(
            '%s:%s' % (attr, value) for attr, value in data.items()))

    @staticmethod
    def ignore(content):
        """Remove dataentry from *content*."""
        page_content = []
        start = False
        for line in content.split('\n'):
            if line == '----' and not start:
                start = True
                continue
            if start:
                page_content.append(line)
        return '\n'.join(page_content) if page_content else content
