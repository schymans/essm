# -*- coding: utf-8 -*-
"""Test variables."""

import pytest

from essm.variables import Variable
from essm.variables.units import derive_unit, joule, kilogram, meter, second


class demo_variable(Variable):
    """Test variable."""

    default = 1
    unit = meter


class demo_expression_variable(Variable):
    """Test expression variable."""

    expr = 2 * demo_variable


def test_variable_definition():
    """Test variable definition."""
    assert demo_variable.definition.default == 1
    assert demo_expression_variable.subs(Variable.__expressions__) \
        == 2 * demo_variable
    assert demo_expression_variable.definition.unit == meter


def test_local_definition():
    """Test local variable definition."""

    class local_definition(Variable):
        """Local definition."""

        class local_variable(Variable):
            """Local variable."""
            unit = meter
            default = 2

        expr = 1.5 * local_variable

    assert local_definition.__doc__ == local_definition.definition.__doc__
    assert local_definition.definition.expr.subs(Variable.__defaults__) == 3


def test_unit_check():
    """Test unit validation."""

    class valid_unit(Variable):
        expr = 3 * demo_variable
        unit = meter

    with pytest.raises(ValueError):

        class invalid_unit(Variable):
            expr = 4 * demo_variable
            unit = second


def test_derive_unit():
    """Test derive_unit from expression."""

    class lambda_E(Variable):
        unit = joule / kilogram

    class E_l(Variable):
        unit = joule / (meter ** 2 * second)

    assert derive_unit(2 * lambda_E * E_l) \
        == kilogram * meter ** 2 / second ** 5
    assert derive_unit(E_l/E_l) == 1

    class dimensionless(Variable):
        expr = demo_variable / demo_expression_variable

    assert dimensionless.definition.unit == 1


def test_remove_variable_from_registry():
    """Check is the variable is removed from registry."""

    class removable(Variable):
        """Should be removed."""

    with pytest.warns(UserWarning):
        del Variable[removable]

    assert removable not in Variable.__registry__

    with pytest.raises(KeyError):
        del Variable[removable]
