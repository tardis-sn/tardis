from django import forms

from crispy_forms.helper import FormHelper


class DynamicForm(forms.Form):
    """
    Base class for a dynamic form generation.
    """

    def __init__(self, *args, **kwargs):
        self.legend, fields = kwargs.pop('legend'), kwargs.pop('fields')
        super(DynamicForm, self).__init__(*args, **kwargs)
        self.helper = FormHelper()
        self.helper.form_tag = False
        self.helper.form_class = 'form-horizontal'
        self.helper.label_class = 'col-lg-2'
        self.helper.field_class = 'col-lg-8'
        if isinstance(fields, forms.Field):
            fields = {'Input': fields}
        for field_name, field in fields.iteritems():
            self.fields[field_name] = field
