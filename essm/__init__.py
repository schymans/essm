# -*- coding: utf-8 -*-
#
# This file is part of essm.
# Copyright (C) 2017 ETH Zurich, Swiss Data Science Center.
#
# essm is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# essm is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with essm; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307, USA.
"""Variables module to deal with physical variables and units.

It allows attaching docstrings to variable names, defining their domains
(e.g. integer, real or complex), their units and LaTeX representations.
You can also provide a default value, which is particularly useful for
physical constants.

Creating variables
==================

To create custom variables, first import ``Variable``:

    >>> from essm.variables import Variable

To define units, you can either import these units from the library,
e.g.

``from essm.variables.units import joule, kelvin, meter``

or define the appropriate units from the SageMath units package, e.g.

::

    from sage.symbolic.units import units
    joule = units.energy.joule

.. doctest:: ipython2

    from sage.symbolic.units import units
    joule = units.energy.joule
    kelvin = units.temperature.kelvin
    meter = units.length.meter
    mole = units.amount_of_substance.mole
    pascal = units.pressure.pascal

Then you can define a custom variable with its name, description,
domain, latex\_name, unit, and an optional default value, e.g.:

.. doctest:: ipython2

    class R_mol(Variable):
        '''Molar gas constant.'''
        unit = joule/(kelvin*mole)
        latex_name = 'R_{mol}'
        default = 8.314472

The variables defined above hold information about their docstring,
units, latex representations and default values if any. Each can be
accessed by e.g.:

.. testcode:: ipython2

    print R_mol.__doc__
    print R_mol.definition.unit
    print R_mol.definition.latex_name
    print R_mol.definition.default


.. testoutput::

    Molar gas constant.
    joule/(kelvin*mole)
    R_{mol}
    8.31447200000000


We will now define a few additional variables.

.. doctest:: ipython2

    class P_g(Variable):
        '''Pressure of gas.'''
        unit = pascal
    
    class V_g(Variable):
        '''Volume of gas.'''
        unit = meter^3
        
    class n_g(Variable):
        '''Amount of gas.'''
        unit = mole
    
    class T_g(Variable):
        '''Temperature of gas.'''
        unit = kelvin
        
    class P_wa(Variable):
        '''Partial pressure of water vapour in air.'''
        unit = pascal
        latex_name = 'P_{wa}'

Variables with expressions as definitions
-----------------------------------------

.. doctest:: ipython2

    class Delta_Pwa(Variable):
        '''Slope of saturated vapour pressure, $\partial P_{wa} / \partial T_g'''
        expr = P_wa(T_g).diff(T_g)
        #unit = pascal/kelvin
        latex_name = r'\Delta'

.. doctest:: ipython2

    Delta_Pwa.definition.unit



.. parsed-literal::

    D[0](P_wa)(T_g*kelvin)/diff(P_wa(T_g), T_g)



.. doctest:: ipython2

    Delta_Pwa.definition.expr




.. parsed-literal::

    diff(P_wa(T_g), T_g)



Unfortunately, the units check for the differential equation does not
return a proper unit, as the variables do not cancel out:

.. doctest:: ipython2

    Delta_Pwa.definition.unit




.. parsed-literal::

    D[0](P_wa)(T_g*kelvin)/diff(P_wa(T_g), T_g)



The above should return instead:

::

    pascal/kelvin

This is a bug in ``essm.bases.expand_units()`` that needs to be fixed.

Creating equations
==================

To create custom equations, first import ``Equation``:

    >>> from essm.equations import Equation

We will now define an equation representing the ideal gas law, based on
the variables defined above:

.. doctest:: ipython2

    class eq_ideal_gas_law(Equation):
        '''Ideal gas law.'''
        
        expr = P_g*V_g == n_g*R_mol*T_g

Note that whenever an equation is defined, its units are checked for
consistency in the background and if they are not consistent, an error
message will be printed. To illustrate this, we will try to define the
above equation again, but omit temperature on the right hand side:

.. testcode:: python

   try:
      class eq_ideal_gas_law(Equation):
          '''Ideal gas law.'''

          expr = P_g*V_g == n_g*R_mol
   except Exception, exc1:
      print exc1

.. testoutput::

   Invalid expression units: kilogram*meter^2/second^2 == kilogram*meter^2/(kelvin*second^2)


The equation can be displayed in typesetted form, and the documentation
string can be accessed in a similar way as for Variable:

.. doctest:: python

   eq_ideal_gas_law.show()
   print eq_ideal_gas_law.__doc__



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{P_g} {V_g} = {R_{mol}} {T_g} {n_g}</script></html>


.. parsed-literal::

   Ideal gas law.


New equation based on manipulation of previous equations
--------------------------------------------------------

We can use the above equation just as any SageMath expression, and e.g.
solve it for pressure:

.. testcode:: ipython2

   soln = solve(eq_ideal_gas_law, P_g); print soln


.. testoutput::

   [
   P_g == R_mol*T_g*n_g/V_g
   ]


If we want to define a new equation based on a manipulation of
eq\_ideal\_gas\_law we can specify that the parent of the new equation
is ``eq_ideal_gas_law.definition``:

.. doctest:: ipython2

   class eq_Pg(eq_ideal_gas_law.definition):
       '''Calculate pressure of ideal gas.'''

       expr = soln[0]
   eq_Pg.show()



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{P_g} = \frac{{R_{mol}} {T_g} {n_g}}{{V_g}}</script></html>


To see the inheritance of the newly created equation:

.. testcode:: ipython2

   eq_Pg.definition.__bases__




.. testoutput::

   (<class '__main__.eq_ideal_gas_law'>,)



Empirical equations with internal variables
-------------------------------------------

Empirical equations not only contain variables but also numbers. As an
example, we will try to define the Clausius-Clapeyron equation for
saturation vapour pressure in the following example, after defining a
few additional variables used in this equation.

.. doctest:: ipython2

   kilogram = units.mass.kilogram

   class lambda_E(Variable):
       '''Latent heat of evaporation.'''
       unit = joule/kilogram
       latex_name = '\\lambda_E'
       default = 2.45e6   

   class M_w(Variable):
       '''Molar mass of water.'''
       unit = kilogram/mole
       default = 0.018

.. testcode:: ipython2

   try:
       class eq_Pwa_CC(Equation):
           '''Clausius-Clapeyron P_wa as function of T_g. 

           \cite[Eq. B3]{hartmann_global_1994}
           '''
    
           expr = P_wa == 611.*e**(-M_w*lambda_E*(1/T_g - 1/273.)/R_mol)
   except Exception, exc1:
       print exc1


.. testoutput::

   Invalid expression units: kilogram/(meter*second^2) == 1.0*e^(0.003663003663003663*M_w*kelvin*lambda_E/R_mol - 0.003663003663003663*M_w*lambda_E/R_mol)


.. doctest:: ipython2

   expr = P_wa == 611*e**(-M_w*lambda_E*(1/T_g - 1/273)/R_mol)
   expr.show()



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{P_{wa}} = 611 \, e^{\left(-\frac{{M_w} {\lambda_E} {\left(\frac{273}{{T_g}} - 1\right)}}{273 \, {R_{mol}}}\right)}</script></html>


The unit mismatch reported in the error message stems from the fact that
the numbers in the empirical equation actually need units. Since the
term in the exponent has to be non-dimensional, the units of ``611``
must be the same as those of ``P_wa``, i.e. pascal. The units of the
subtraction term in the exponent must match, meaning that ``273`` needs
units of kelvin. To avoid the error message, we can define the empirical
numbers as internal variables to the equation we want to define:

.. doctest:: ipython2

   class eq_Pwa_CC(Equation):
       '''Clausius-Clapeyron P_wa as function of T_g. 
   
       Eq. B3 in :cite{hartmann_global_1994}
       '''
           
       class p_CC1(Variable):
           '''Internal parameter of eq_Pwl.'''
           unit = pascal
           latex_name = '611'
           default = 611.   
       
        
        
       class p_CC2(Variable):
           '''Internal parameter of eq_Pwl.'''
           unit = kelvin
           latex_name = '273'
           default = 273.   
        
       expr = P_wa == p_CC1*e**(-M_w*lambda_E*(1/T_g - 1/p_CC2)/R_mol)

In the above, we defined the latex representation of the empirical
constants as their actual values, so the equation displays in the
familiar way:

.. doctest:: ipython2

   eq_Pwa_CC.show()



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{P_{wa}} = {611} e^{\left(-\frac{{M_w} {\lambda_E} {\left(\frac{1}{{T_g}} - \frac{1}{{273}}\right)}}{{R_{mol}}}\right)}</script></html>


All default values of variables defined along with the variable
definitions are stored in a dictionary that can be accessed as
``Variable.__defaults__``. We can substitute the values from this
dictionary into our empirical equation to plot saturation vapour
pressure as a function of temperature:

.. doctest:: ipython2

   print eq_Pwa_CC.subs(Variable.__defaults__)
   P = plot(eq_Pwa_CC.rhs().subs(Variable.__defaults__), (T_g, 273, 373), frame=True, axes=False)
   P.axes_labels(['Air temperature (K)', 'Saturation vapour pressure (Pa)'])
   P.show(dpi=60)


.. parsed-literal::

   P_wa == 611.000000000000*e^(-5304.00487246815/T_g + 19.4285892764401)



Use of functional notation
==========================

Above, we defined ``Delta_Pwa`` as a variable that represents the
partial derivative of ``P_wa`` with respect to ``T_g``:

::

   class Delta_Pwa(Variable):
       '''Slope of saturated vapour pressure, $\partial P_{ws} / \partial T_g'''
       expr = P_wa(T_g).diff(T_g)
       #unit = pascal/kelvin
       latex_name = r'\Delta'

This definition can be accessed by typing ``Delta_Pwa.definition.expr``.
Example:

.. doctest:: ipython2

   print Delta_Pwa.definition.expr
   show(Delta_Pwa == Delta_Pwa.definition.expr)


.. parsed-literal::

   diff(P_wa(T_g), T_g)



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{\Delta} = \frac{\partial}{\partial {T_g}}P_{{\rm wa}}\left({T_g}\right)</script></html>


We also defined the Clausius-Clapeyron approximation to
:math:`P_{wa}(T_g)` as ``eq_Pwa_CC``.

.. doctest:: ipython2

   eq_Pwa_CC.show()
   print eq_Pwa_CC.__doc__



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}{P_{wa}} = {611} e^{\left(-\frac{{M_w} {\lambda_E} {\left(\frac{1}{{T_g}} - \frac{1}{{273}}\right)}}{{R_{mol}}}\right)}</script></html>


.. parsed-literal::

   Clausius-Clapeyron P_wa as function of T_g. 
   
       Eq. B3 in :cite{hartmann_global_1994}
        


If we want to substitute this approximation into
``Delta_Pwa.definition.expr``, we need to use
``substitute_function(function, callable symbolic expression)``.
``P_wa.definition.function`` returns a symbolic function called
``P_wa``, while ``eq_Pwa_cc.rhs().function(T_g)`` returns the right hand
side of ``eq_Pwa_cc`` as a callable symbolic expression with argument
``T_g``:

.. testcode:: ipython2

   print Delta_Pwa.definition.function
   print 'type of Delta_Pwa: ' + str(type(Delta_Pwa.definition.function))
   print eq_Pwa_CC.rhs().function(T_g)

.. testoutput::

   Delta_Pwa
   type of Delta_Pwa: <class 'sage.symbolic.function_factory.NewSymbolicFunction'>
   T_g |--> p_CC1*e^(-M_w*lambda_E*(1/T_g - 1/p_CC2)/R_mol)


.. doctest:: ipython2

   exprDeltaPwa = Delta_Pwa.definition.expr.substitute_function(
       P_wa.definition.function, eq_Pwa_CC.rhs().function(T_g))
   show(exprDeltaPwa)



.. raw:: html

   <html><script type="math/tex; mode=display">\newcommand{\Bold}[1]{\mathbf{#1}}\frac{{M_w} {\lambda_E} {611} e^{\left(-\frac{{M_w} {\lambda_E} {\left(\frac{1}{{T_g}} - \frac{1}{{273}}\right)}}{{R_{mol}}}\right)}}{{R_{mol}} {T_g}^{2}}</script></html>


We can now plot ``expDeltaPwa`` for a range of values for ``T_g``:

.. doctest:: ipython2

   print exprDeltaPwa.subs(Variable.__defaults__)
   P = plot(exprDeltaPwa.subs(Variable.__defaults__), (T_g, 273, 373), frame=True, axes=False)
   P.axes_labels(['Air temperature (K)', 'Slope of saturation vapour pressure, $\Delta$\n (Pa K$^{-1}$)'])
   P.show(dpi=60)


.. parsed-literal::

   3.24074697707804e6*e^(-5304.00487246815/T_g + 19.4285892764401)/T_g^2





Importing variables and equations
=================================
You can import pre-defined variables and equations as e.g.:

    >>> from essm.variables.physics.thermodynamics import *
    >>> from essm.equations.physics.thermodynamics import *
"""

from __future__ import absolute_import

from .bases import convert, expand_units

__all__ = ('convert', 'expand_units')
