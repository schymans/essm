"""Generator for equation definitions."""

import re

from collections import defaultdict

from .variables import Variable

EQUATION_TPL = """
class {name}({parents}):
    \"\"\"{doc}\"\"\"
    {variables}
    expr = {expr}
"""

VARIABLE_TPL = """
class {name}(Variable):
    \"\"\"{doc}\"\"\"
    name = \'{name}\'
    unit = {units}
    latex_name = r\"{latexname}\"
    {default}   
"""


class VariableWriter(object):
    """Generate Variable definitions.

    Example:

    .. code-block:: python
        from essm._generator import VariableWriter


    """

    TPL = VARIABLE_TPL
    default_imports = {
        'essm.variables': {'Variable'},
    }

    def __init__(self, docstring=None):
        self.docstring = docstring
        self._imports = defaultdict(set)
        self._imports.update(**self.default_imports)
        self.vars = []

    @property
    def imports(self):
        for key, values in self._imports.items():
            yield 'from {key} import {names}'.format(
                key=key, names=', '.join(sorted(values)))

    def write(self, filename='temp/eqs_test.py'):
        with open(filename, 'w') as file_out:
            if self.docstring:
                file_out.write('"""' + self.docstring + '"""\n\n')
            file_out.write('\n'.join(self.imports) + '\n')
            file_out.write('\n\n'.join(
                self.TPL.format(**var).replace('^', '**')
                for var in self.vars))
            file_out.write(
                '\n\n__all__ = (\n{0}\n)'.format('\n'.join(
                    "    '{0}',".format(var['name']) for var in self.vars)))

    def var(self,
            name,
            doc='',
            units=None,
            domain1='real',
            latexname=None,
            value=None):
        if not latexname:
            latexname = name
        if value == None:
            default = ''
        else:
            default = 'default = ' + str(value)
            # Skip trailing zeroes from real numbers only
            if type(value) == type(0.1):
                default = 'default = ' + value.str(skip_zeroes=True).replace(
                    '^', '**')
        context = {
            "name": name,
            "doc": doc,
            "units": str(units).replace('^', '**') if units else '1/1',
            "domain1": domain1,
            "latexname": latexname,
            "default": default
        }
        self.vars.append(context)

        # register all imports of units
        if units:
            if units != 1:
                for arg in units.args():
                    self._imports['essm.variables.units'].add(str(arg))


class EquationWriter(object):
    """Generate Equation definitions.

    Example:

    .. code-block:: python
        from essm.variables import Variable
        from essm.equations import Equation
        from essm._generator import EquationWriter
        from essm.variables.units import second, meter, kelvin
        from essm.variables.physics.thermodynamics import R_s, D_va, T_a
        from essm.variables.leaf.unsorted import R_ll, H_l, E_l
        var('R_s R_ll H_l E_l D_va T_a p_Dva1 p_Dva2')
        writer = EquationWriter(docstring="Test.")
        writer.eq('eq_enbal', 0 == R_s - R_ll - H_l - E_l, doc='Energy balance.')
        writer.eq('eq_Rs_enbal', R_s == R_ll + H_l + E_l, doc='Calculate R_s from energy balance.', parents=['eq_enbal'])
        writer.eq('eq_Dva', D_va == p_Dva1*T_a - p_Dva2, doc='D_va as a function of air temperature'
                , variables = [{"name": "p_Dva1", "default": '1.49e-07', "units": meter^2/second/kelvin}, \
                {"name": "p_Dva2", "default": '1.96e-05', "units": meter^2/second}])
        writer.write(filename='temp/test.py')
    """

    TPL = EQUATION_TPL
    VAR_TPL = VARIABLE_TPL
    default_imports = {
        'essm.equations': {'Equation'},
    }

    def __init__(self, docstring=None):
        self.docstring = docstring
        self._imports = defaultdict(set)
        self._imports.update(**self.default_imports)
        self.eqs = []

    @property
    def imports(self):
        for key, values in self._imports.items():
            yield 'from {key} import {names}'.format(
                key=key, names=', '.join(sorted(values)))

    def write(self, filename='temp/eqs_test.py'):
        with open(filename, 'w') as file_out:
            if self.docstring:
                file_out.write('"""' + self.docstring + '"""\n\n')
            file_out.write('\n'.join(self.imports) + '\n')
            file_out.write('\n'.join(
                self.TPL.format(**eq).replace('^', '**')
                for eq in self.eqs))
            file_out.write('\n\n__all__ = (\n{0}\n)'.format(
                '\n'.join("    '{0}',".format(eq['name']) for eq in self.eqs)))

    def eq(self, name, expr, doc='', parents=None, variables=None):
        if parents:
            parents = ', '.join(parent + '.definition' for parent in parents)
        else:
            parents = 'Equation'

        if variables:
            for variable in variables:
                variable.setdefault('latexname', variable['name'])
                variable['doc'] = "Internal parameter of {0}.".format(
                    variable['name'])
                if 'default' in variable:
                    variable['default'] = 'default = {0}'.format(
                        variable['default'])
                else:
                    variable['default'] = ''
            variables = '\n'.join(re.sub(
                r'^', 4 * ' ',
                self.VAR_TPL.format(**variable),
                flags=re.MULTILINE,
            ) for variable in variables)
        else:
            variables = ''

        context = {
            "name": name,
            "doc": doc,
            "expr": expr,
            "parents": parents,
            "variables": variables,
        }
        self.eqs.append(context)

        # register all imports
        for arg in expr.args():
            if arg in Variable.__registry__:
                self._imports[Variable.__registry__[arg].__module__].add(
                    str(arg))
            else:
                self._imports['essm.variables'].add('Variable')
