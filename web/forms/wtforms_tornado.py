# https://github.com/puentesarrin/wtforms-tornado
# This code is licensed under the Apache 2 License 
# (http://www.apache.org/licenses/LICENSE-2.0)

# -*- coding: utf-8 -*-
import re

from tornado import escape
from wtforms import form
from wtforms.compat import text_type

class TornadoInputWrapper(object):

    def __init__(self, multidict):
        self._wrapped = multidict

    def __iter__(self):
        return iter(self._wrapped)

    def __len__(self):
        return len(self._wrapped)

    def __contains__(self, name):
        return (name in self._wrapped)

    def __getitem__(self, name):
        return self._wrapped[name]

    def __getattr__(self, name):
        return self.__getitem__(name)

    def getlist(self, name):
        try:
            return [escape.to_unicode(v) for v in self._wrapped[name]]
        except KeyError:
            return []


class Form(form.Form):
    """
    A Form derivative which uses the locale module from Tornado.
    """

    def __init__(self, formdata=None, obj=None, prefix='', locale_code='en_US',
             **kwargs):
        self._locale_code = locale_code
        super(Form, self).__init__(formdata, obj, prefix, **kwargs)

    def process(self, formdata=None, obj=None, **kwargs):
        if formdata is not None and not hasattr(formdata, 'getlist'):
            formdata = TornadoInputWrapper(formdata)
        super(Form, self).process(formdata, obj, **kwargs)
