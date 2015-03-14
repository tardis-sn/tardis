import yaml

from django.shortcuts import render
from django.http import HttpResponse

from .schema import generate_schema
from .forms import DynamicForm


def home(request):
    if request.method == 'POST':
        tardis_schema = generate_schema()
        config_version_form = DynamicForm(
            request.POST,
            legend='tardis_config_version',
            fields=tardis_schema['tardis_config_version'],
        )
        supernova_form = DynamicForm(
            request.POST,
            legend='supernova',
            fields=tardis_schema['supernova'],
        )
        atom_data_form = DynamicForm(
            request.POST,
            legend='Atom Data',
            fields=tardis_schema['atom_data'],
        )
        nlte_form = DynamicForm(
            request.POST,
            legend='NLTE',
            fields=tardis_schema['plasma'].pop('nlte'),
        )
        plasma_form = DynamicForm(
            request.POST,
            legend='Plasma',
            fields=tardis_schema['plasma'],
        )
        spectrum_form = DynamicForm(
            request.POST,
            legend='spectrum',
            fields=tardis_schema['spectrum']
        )
        _gen_schema = {}
        if config_version_form.is_valid():
            _gen_schema['tardis_config_version'] = {}
            print config_version_form.cleaned_data
            _gen_schema['tardis_config_version'] = \
                config_version_form.cleaned_data['Input']
        if supernova_form.is_valid():
            _gen_schema['supernova'] = {}
            for field in supernova_form.fields:
                _gen_schema['supernova'][field] = \
                    supernova_form.cleaned_data[field]
        if atom_data_form.is_valid():
            _gen_schema['atom_data'] = \
                atom_data_form.cleaned_data['Input']
        if plasma_form.is_valid():
            _gen_schema['plasma'] = {}
            for field in plasma_form.fields:
                _gen_schema['plasma'][field] = \
                    plasma_form.cleaned_data[field]
        if nlte_form.is_valid():
            _gen_schema['plasma']['nlte'] = {}
            for field in nlte_form.fields:
                _gen_schema['plasma']['nlte'][field] = \
                    nlte_form.cleaned_data[field]
        if spectrum_form.is_valid():
            _gen_schema['spectrum'] = \
                spectrum_form.cleaned_data['Input']

        response = HttpResponse(yaml.safe_dump(
            _gen_schema,
            default_flow_style=False
            )
        )
        response['Content-Disposition'] = 'attachment; filename="tardis.yml"'
        return response

    else:
        tardis_schema = generate_schema()
        config_version_form = DynamicForm(
            legend='tardis_config_version',
            fields=tardis_schema['tardis_config_version'],
        )
        supernova_form = DynamicForm(
            legend='supernova',
            fields=tardis_schema['supernova'],
        )
        atom_data_form = DynamicForm(
            legend='Atom Data',
            fields=tardis_schema['atom_data'],
        )
        nlte_form = DynamicForm(
            legend='NLTE',
            fields=tardis_schema['plasma'].pop('nlte'),
        )
        plasma_form = DynamicForm(
            legend='Plasma',
            fields=tardis_schema['plasma'],
        )
        spectrum_form = DynamicForm(
            legend='spectrum',
            fields=tardis_schema['spectrum']
        )
        return render(request, 'tardis/home.html', {
            'config_version_form': config_version_form,
            'supernova_form': supernova_form,
            'atom_data_form': atom_data_form,
            'plasma_form': plasma_form,
            'nlte_form': nlte_form,
            'spectrum_form': spectrum_form,
            }
        )
