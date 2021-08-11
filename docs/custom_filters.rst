Custom Filters
==============
Overview
--------
Pickaxe expansions can grow extremely quickly in size, resulting in more compounds and reactions than a computer can
efficiently handle, both during the expansion and during the subsequent analysis.
An option to limit the chemical space explored during an expansion is to create a filter that selects
only a subset of compounds to react in each generation.
For example, you could create a filter that only allows compounds below a certain molecular weight or only compounds with
specific structural features to react. A few filters have already been written by the Tyo lab, which you can find at :ref:`Built-In Filters`.
We recommend looking at these as examples of how to write a filter as you write your own.

By creating and using a custom filter, you can control the scope of your expansion, allowing you to also expand for more generations.
It also saves space in the database and should make downstream analysis faster.

Requirements
------------
Creating a custom filter requires a working knowledge of python. Default filters are created using [RDKit](https://rdkit.org/docs/api-docs.html), a python library
providing a collection of cheminformatic tools.

Ensure that that you have the [MINE-Database](https://github.com/tyo-nu/MINE-Database) github cloned on your machine.

The overall process for creating a filter is as follows:

#. Write custom Filter subclass in minedatabase/filters.py
#. Expose options for this filter subclass and add it to a pickaxe run in pickaxe_run.py
#. (optional) Write unit test(s) for this custom filter in tests/test_unit/test_filters.py

Writing Custom Filters
----------------------
To write a custom filter, you need to subclass the Filter class in filters.py. The Filter class specifies
the required functions your filter needs to implement as well as providing default methods that are inherited
by your custom filter.

There are three methods you must implement (at a minimum). These are the __init__ method and the methods that are decorated with the
@abc.abstractmethod decorator.

1. **__init__** - Initialize your filter's options and inputs here.

2. **filter_name** - This method just needs to return the filter name. This can be set to a constant string, set to a permanent self._filter_name (as in TanimotoSamplingFilter), or set to a custom self._filter_name (as in MetabolomicsFilter).

3. **_choose_cpds_to_filter** - This is the main method you need to implement, where you can loop through the compounds at each generation and decide which ones to keep and which ones to filter out. See the built-in filters' implementations of this method for examples. This method needs to return a set of compound IDs (e.g. "Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f") that you wish to **remove** from the expansion. Note that if a compound is a side-product of a reaction producing a *kept* compound, that compound will not be removed from the expansion, it just won't be reacted further.

There are two optional methods you may override as well. See the Filter ABC in filters.py (:doc:`api_ref/filters`) for more details.

1. **_pre_print** - This method prints to stdout just before the filter is applied. It should return None.

2. **_post_print** - This method prints to stdout just after the filter is applied. Useful for printing a summary of filtering results. It should return None.

Using Your Filter in a Pickaxe Run
----------------------------------
Now that you have a filter defined, the next step is to import it and use it in a Pickaxe run.
Refer to the example file, :download:`pickaxe_run.py <_static/files/pickaxe_run.py>`, which is detailed more in :doc:`pickaxe_run`
to see an example of a Pickaxe run that uses filters. If you open pickaxe_run.py, you will notice different sections for the various built-in filters. Initialize your filter with any options
that you have defined and then ensure you are appending your filter object to the pickaxe object.

You can find this in pickaxe_run.py by scrolling down to the comment that says "# Apply filters". The default filters all
have an if statement associated with them when the filter is defined earlier in the file. Either replicate this format, or simply append
your filter to the pickaxe object

.. code-block:: python

    pk.filters.append(my_filter)

That's it! Now, pickaxe will use that filter during any expansions.

If you have written tests for your filter and think it could be valuable to the community, feel free to submit a pull request at https://github.com/tyo-nu/MINE-Database to add your filter to the built-in set.

(Optional) Writing Tests for Your Filter
-----------------------------
While it is not necessary, it is a good idea to write filters for your test to ensure the behavior of your tests don't
change in the event of an update. There is an already existing file located at `tests/test_unit/test_filters.py` that you can add
your tests to. We utilize [pytest](https://docs.pytest.org/en/stable/) and have defined useful fixtures for use in the tests.
To run these tests run the following from the base MINE-Database directory

... codeblock::

    pytest tests/test_unit/test_filters.py
