Thermodynamic Calculations
==========================
Overview of Thermodynamics Module
----------------------------------
eQuilibrator
^^^^^^^^^^^^
Built into Pickaxe is the ability to estimate the Gibbs free energy of compounds and reactions.
Pickaxe uses eQuilibrator to calculate thermodynamic values. More information about
`eQuilibrator can be found here <https://equilibrator.weizmann.ac.il/static/classic_rxns/about.html>`_.


Calculable Values
^^^^^^^^^^^^^^^^^
Pickaxe calculates the following values for compounds and reactions. More information
about these `conditions can be found here <https://equilibrator.weizmann.ac.il/static/classic_rxns/faq.html#what-are-rg-rg-and-rg>`_.

#. Compounds
    * : ΔfG’°: The Standard Gibbs Free Energy of Formation
        a. Uses pH = 7 and ionic strength = 0.1M
#. Reactions
    * ΔrG’°: The Standard Gibbs Free Energy of Reaction
        a. Uses pH = 7 and ionic strength = 0.1M
    * ΔrG’m: The Physiological Gibbs Free Energy of Reaction
        a. Uses pH = 7 and ionic strength = 0.1M
        b. Assumes concentrations are 1mM.
    * ΔrG’: Adjusted Gibbs Free Energy of Reaction
        a. User-specified conditions

Calculating Thermodynamics of a Pickaxe Run
-------------------------------------------
Set-up
^^^^^^
Thermodynamics.py uses the compound ids (c_id) and reaction ids (r_id) of pickaxe
runs to calculate values. This example assumes you have run a pickaxe run and have
it accessible either from a MongoDB or in memory in a pickaxe object. Pickaxe runs can
be stored later by using the pickleing functionality.

Additionally, an eQuilibrator database must be loaded. 

Compound Value Calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
If there is no eQuilibrator compounds.sqlite file present, generate one first.
::

    >>> from equilibrator_assets.local_compound_cache import LocalCompoundCache
    >>> lc = LocalCompoundCache()
    >>> lc.generate_local_cache_from_default_zenodo("compounds.sqlite")
    Copying default Zenodo compound cache to compounds.sqlite

Next, the thermodynamics class must be loaded and initialized, where `mongo_uri` is the
uri to your mongo server. Providing None will use the default localhost.
::

    >>> from minedatabase.thermodynamics import Thermodynamics
    >>> thermo = Thermodynamics()
    >>> thermo.load_thermo_from_sqlite("compounds.sqlite")
    Loading compounds from compounds.sqlite
    >>> thermo.load_mongo(mongo_uri=mongo_uri)
    
The following assumes you have a valid pickaxe object or a database to cross-reference
the c_id and r_id from. No c_id or r_id is given here, but example outputs are.

Calculating ∆Gf'°
::

    >>> thermo.standard_dg_formation_from_cid(c_id=c_id, pickaxe=pickaxe, db_name=db_name)
    -724.5684043965385

Calculating ΔrG’m
::

    >>> thermo.physiological_dg_prime_from_rid(r_id, pickaxe=pickaxe, db_name=db_name)
    <Measurement(-5.432945798382008, 3.37496192184388, kilojoule / mole)>

Calculating ΔrG’ at pH = 4, ionic strength = 0.05M
::

    >>> from equilibrator_api import Q_
    >>> p_h = Q_("4")
    >>> ionic_strength = Q_("0.05M)
    >>> thermo.dg_prime_from_rid(r_id=r_id, db_name=db_name, p_h=p_h, ionic_strength=ionic_strength)
    <Measurement(11.68189173633911, 3.37496192184388, kilojoule / mole)>

