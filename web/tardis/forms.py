from django import forms


class DynamicForm(forms.Form):
    """
    Base class for a dynamic form generation.
    """

    def __init__(self, *args, **kwargs):
        self.legend, fields = kwargs.pop('legend'), kwargs.pop('fields')
        super(DynamicForm, self).__init__(*args, **kwargs)
        if isinstance(fields, forms.Field):
            fields = {'Input': fields}
        for field_name, field in fields.iteritems():
            self.fields[field_name] = field
