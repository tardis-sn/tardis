from wtforms import validators
from wtforms.fields import *
from wtforms.widgets import *
from wtforms_tornado import Form
from wtforms.validators import ValidationError

mandatory = validators.DataRequired()

class SupernoveForm(Form):
    luminosity_requested = TextField('Luminosity Requested', 
                                     validators=[mandatory],
                                     description={
                                     'help_text':'requested output luminosity for simulation',
                                     'default':'1 solLum'})
    time_explosion = TextField('Time since Explosion', 
                               validators=[mandatory])
    distance = TextField('Distance from supernova')
    luminosity_wavelength_start = TextField('Luminosity Wavelength (Start)',
                                          description={
                                          'help_text':'start of the integral needed for getting the luminosity right',
                                          'default':'0 angstrom'})
    luminosity_wavelength_end = TextField('Luminosity Wavelength (End)',
                                          description={
                                          'help_text':'end of the integral needed for getting the luminosity right',
                                          'default':'inf angstrom'})

#------------------------------------------------------------------------------

class AtomForm(Form):
    atom_data = FileField(u'Atom Data', [mandatory])

#------------------------------------------------------------------------------

class PlasmaForm(Form):
    intial_t_inner = TextField('Initial Inner Temperature',
                               description={
                               'help_text':'initial temperature of the inner boundary black body. If set to -1 K will result in automatic calculation of boundary',
                               'default':'-1 K'})
    intial_t_rad = TextField('Initial Radiative Temperature', 
                             description={'help_text':'initial radiative temperature in all cells (if set)',
                             'default':'1000 K'})
    disable_electron_scattering = BooleanField(description={
                                               'help_text':'disable electron scattering process in montecarlo part - non-physical only for tests'})
    ionization_treatment_mode = SelectField(validators=[mandatory],
                                            choices = [('nebular','Nebular'),
                                                       ('lte','LTE')])
    excitation_treatment_mode = SelectField(validators=[mandatory],
                                            choices = [('lte','LTE'),
                                                       ('dilute-lte','Dilute LTE')])
    radiative_rates_treatment_mode = SelectField(validators=[mandatory],
                                                 choices = [('dilute-blackbody','Dilute Blackbody'),
                                                            ('detailed','Detailed')])
    line_interaction_mode = SelectField(validators=[mandatory],
                                        choices = [('scatter','Scatter'),
                                                   ('downbranch','Downbranch'),
                                                   ('macroatom','Macro Atom')])
    w_epsilon = TextField(description={
                          'help_text':'w to use when j_blues get numerically 0. - avoids numerical complications',
                          'default':'-1e-10'})
    delta_treatment = TextField(description={
                                'help_text':'In the saha calculation set delta equals to the number given in this configuration item. if set to None (default), normal delta treatment (as described in Mazzali &amp; Lucy 1993) will be applied'})
    nlte_options = BooleanField('Configure NLTE Options')
    nlte_options_species = TextField('- Species',
                                     description={'parent':'nlte_options',
                                     'help_text':'Species that are requested to be NLTE treated in the format [\'Si 2\', \'Ca 1\', etc.]'})
    nlte_options_coronal_approximation = BooleanField('- Coronal Approximation',
                                                      description={'parent':'nlte_options',
                                                      'help_text':'set all jblues=0.0',
                                                      'default':'False'})
    nlte_options_classical_nebular = BooleanField('- Classical Nebular',
                                                  description={'parent':'nlte_options',
                                                  'help_text':'sets all beta_sobolevs to 1',
                                                  'default':'False'})

#------------------------------------------------------------------------------

class StructureFileForm(Form):
    filename = FileField(u'Upload File',[mandatory])
    filetype = SelectField(validators=[mandatory],
                           choices = [('simple_ascii','Simple ASCII'),
                                      ('artis','Artis')])
    v_inner_boundary = TextField(description={
                                 'help_text':'location of the inner boundary chosen from the model',
                                 'default':'0 km/s'})
    v_outer_boundary = TextField(description={
                                 'help_text':'location of the outer boundary chosen from the model',
                                 'default':'inf km/s'})

class StructureSpecificForm(Form):
    velocity = HiddenField('Velocity')
    velocity_start = TextField('- Start',
                               validators = [mandatory])
    velocity_stop = TextField('- Stop',
                              validators = [mandatory])
    velocity_num = TextField('- Num',
                             validators = [mandatory])
    density = HiddenField('Density')
    density_type = SelectField('- Type',
                               validators=[mandatory],
                               choices = [('branch85_w7','branch85_w7'),
                                          ('exponential','Exponential'),
                                          ('power_law','Power Law'),
                                          ('uniform','Uniform')])

    density_time_0 = TextField('- Time(0)',
                               validators = [mandatory],
                               description={'help_text':'Time at which the pure model densities are right',
                                            'association':'density_type power_law exponential'})
    density_rho_0 = TextField('- Rho(0)',
                              validators = [mandatory],
                              description={'help_text':'density at Time(0)',
                                           'association':'density_type power_law exponential'})
    density_v_0 = TextField('- V(0)',
                            validators = [mandatory],
                            description={'help_text':'at what velocity the density Rho(0) applies',
                                         'association':'density_type power_law exponential'})
    density_exponent = TextField('- Exponent',
                                 validators = [mandatory],
                                 description={'help_text':'exponent for exponential density profile',
                                              'association':'density_type power_law'})
    density_value = TextField('- Value',
                              validators = [mandatory],
                              description={'help_text':'value for uniform density',
                                           'association':'density_type uniform'})


class StructureForm(Form):
    type = RadioField(choices=[('file','File Upload'),
                               ('specific','Specific Entries')], 
                      validators=[mandatory],
                      default='file')
    file = FormField(StructureFileForm)
    specific = FormField(StructureSpecificForm)

class AbundanceForm(Form):
    filename = FileField(u'Upload File')
    filetype = SelectField(choices = [('artis','Artis'),
                                      ('simple_ascii','Simple ASCII')])
    uniform_abundances = TextAreaField(description={
                                       'help_text':'Insert Uniform abundances of all the shells, in the format: C: 0.01 O: 0.01 etc...'})

class ModelForm(Form):
    model_structure = FormField(StructureForm)
    model_abundance = FormField(AbundanceForm)

#------------------------------------------------------------------------------

class MonteCarloForm(Form):
    Seed = TextField(description={'help_text':'Seed for the random number generator',
                     'default':'2311963'})
    no_of_packets = TextField(validators=[mandatory],
                              description={'help_text':'Seed for the random number generator'})
    iterations = TextField(validators=[mandatory],
                           description={'help_text':'Number of maximum iterations'})
    black_body = HiddenField('Blackbody Sampling',
                             description={'help_text':'Sampling of the black-body for energy packet creation (giving maximum and minimum packet frequency)'})
    black_body_sampling_start = TextField('- Start',
                                          description={'default':'50 angstrom'})
    black_body_sampling_stop = TextField('- Stop',
                                         description={'default':'200000 angstrom'})
    black_body_sampling_num = TextField('- Num',
                                        description={'default':'1000000'})
    last_no_of_packets = TextField(description={'help_text':'This can set the number of packets for the last run. If set negative it will remain the same as all other runs.',
                                   'default':'-1'})
    no_of_virtual_packets = TextField(description={'help_text':'Setting the number of virtual packets for the last iteration.',
                                      'default':'0'})
    enable_reflective_inner_boundary = BooleanField(description={'help_text':'experimental feature to enable a reflective boundary.', 
                                                    'default':'False'})
    inner_boundary_albedo = TextField(description={'help_text':'albedo of the reflective boundary',
                                      'default':'0.0'})
#------------------------------------------------------------------------------

class SpectrumForm(Form):
    start = TextField(description={'default':'50'})
    stop = TextField(description={'default':'200000'})
    num = TextField(description={'default':'1000000'})


    convergence_criteria = BooleanField()

    type = SelectField(choices=[('specific','Specific'),
                                ('damped','Damped')],
                       description={'parent':'convergence_criteria'})
    
    t_inner_update_exponent = TextField(description={
                                        'parent':'convergence_criteria',
                                        'help_text':'L=4*pi*r**2*T^y', 
                                        'default':'-0.5'})
    lock_t_inner_cycles = TextField(description={'parent':'convergence_criteria',
                                    'help_text':'The number of cycles to lock the update of the inner boundary temperature. This process helps with convergence. The default is to switch it off (1 cycle)', 
                                    'default':'1'})

    hold_iterations = TextField(validators=[mandatory],
                                description={'parent':'convergence_criteria',
                                'help_text':'The number of iterations that the convergence criteria need to be fulfilled before TARDIS accepts the simulation as converged', 
                                'default':'3',
                                'association':'spectrum_type specific'})
    fraction = TextField(validators=[mandatory],
                         description={'parent':'convergence_criteria',
                         'help_text':'he fraction of shells that have to converge to the given convergence threshold. For example, 0.8 means that 80% of shells have to converge to the threshold that convergence is established', 
                         'default':'0.8',
                         'association':'spectrum_type specific'})
    damping_constant_damped = TextField('Damping Constant',
                                        description={'parent':'convergence_criteria',
                                        'default':'0.5',
                                        'association':'spectrum_type damped'})
    damping_constant_specific = TextField('Damping Constant',
                                          validators=[mandatory], description={'parent':'convergence_criteria',
                                          'default':'0.5',
                                          'association':'spectrum_type specific'})
    threshold = TextField(validators=[mandatory],
                          description={'parent':'convergence_criteria','help_text':'specifies the threshold that is taken as convergence (i.e. 0.05 means that the value does not change more than 5%)', 
                          'default':'1',
                          'association':'spectrum_type specific'})

    t_inner = BooleanField(description={'parent':'convergence_criteria'})
    t_inner_damping_constant = TextField('Damping Constant',
                                         description={'parent':'t_inner',
                                         'default':'0.5'})
    t_inner_threshold = TextField('Threshold',
                                  description={'parent':'t_inner','help_text':'specifies the threshold that is taken as convergence (i.e. 0.05 means that the value does not change more than 5%)', 
                                  'default':'1'})
    t_rad = BooleanField(description={'parent':'convergence_criteria'})
    t_rad_damping_constant = TextField('Damping Constant',
                                       description={'parent':'t_rad',
                                       'default':'0.5'})
    t_rad_threshold = TextField('Threshold',
                                validators=[mandatory],
                                description={'parent':'t_rad','help_text':'specifies the threshold that is taken as convergence (i.e. 0.05 means that the value does not change more than 5%)', 
                                'default':'1'})
    t_w = BooleanField(description={'parent':'convergence_criteria'})
    t_w_damping_constant = TextField('Damping Constant',
                                     description={'parent':'t_w',
                                     'default':'0.5'})
    t_w_threshold = TextField('Threshold',
                              validators=[mandatory],
                              description={'parent':'t_w','help_text':'specifies the threshold that is taken as convergence (i.e. 0.05 means that the value does not change more than 5%)', 
                              'default':'1'})

#------------------------------------------------------------------------------

class TardisForm(Form):
    supernova = FormField(SupernoveForm)
    atom_data = FormField(AtomForm)
    plasma = FormField(PlasmaForm)
    model = FormField(ModelForm)
    montecarlo = FormField(MonteCarloForm)
    spectrum = FormField(SpectrumForm)
