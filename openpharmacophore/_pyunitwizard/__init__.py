import pyunitwizard as puw
puw.configure.load_library(['pint', "openmm.unit"])
puw.configure.set_default_form('pint')
puw.configure.set_standard_units(['angstroms', 'ps', 'K', 'mole', 'amu', 'e',
                                 'kJ/mol', 'kJ/(mol*nm**2)', 'N', 'degrees'])

