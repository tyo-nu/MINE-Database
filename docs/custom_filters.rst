Custom Filters
==============

You can make your own custom filters to filter out compounds after each generation.
For example, you could create a filter to only keep compounds with certain groups or a filter to only keep molecules of a certain size, etc.
A few filters have already been written by the Tyo lab, which you can find at :ref:`Built-In Filters`.
We recommend looking at these as examples of how to write a filter as you write your own.

By creating and using a custom filter, you can limit the size of your expansion, which allows you to expand for more generations if desired.
It also saves space in the database and should make downstream analysis faster.

To create a custom filter, first, make sure that you have the MINE-Database repository cloned on your local machine (i.e. not installed with pip).

The steps are then:

1. Write custom Filter subclass in minedatabase/filters.py
2. (optional) Write unit test(s) for this custom filter in tests/test_unit/test_filters.py
3. Expose options for this filter subclass and add it to a pickaxe run in pickaxe_run.py

Write Custom Filter
-------------------
To write a custom filter, you need to subclass the Filter class in filters.py. This class is an
Abstract Base Class (ABC) which you can think of as a blueprint for other classes. In other words,
it specifies which methods you need to implement in your custom filter so that pickaxe can use it.

There are three methods you must implement (at a minimum). These are the __init__ method and the methods that are decorated with the
@abc.abstractmethod decorator.

1. **__init__** - Initialize your filter's options and inputs here.

2. **filter_name** - This method just needs to return the filter name. This can be set to a constant string, set to a permanent self._filter_name (as in TanimotoSamplingFilter), or set to a custom self._filter_name (as in MetabolomicsFilter).

3. **_choose_cpds_to_filter** - This is the main method you need to implement, where you can loop through the compounds at each generation and decide which ones to keep and which ones to filter out. See the built-in filters' implementations of this method for examples. This method needs to return a set of compound IDs (e.g. "Ccffda1b2e82fcdb0e1e710cad4d5f70df7a5d74f") that you wish to **remove** from the expansion. Note that if a compound is a side-product of a reaction producing a *kept* compound, that compound will not be removed from the expansion, even if it is in the set of compounds to remove.

There are two optional methods you may override as well. See the Filter ABC in filters.py (:doc:`api_ref/filters`) for more details.

1. **_pre_print** - This method prints to stdout just before the filter is applied. It should return None.

2. **_post_print** - This method prints to stdout just after the filter is applied. Useful for printing a summary of filtering results. It should return None.

(Optional) Write Filter Tests
-----------------------------
We recommend writing tests for your filter before using it. To do so, open tests/test_unit/test_filters.py.
Make sure to import your custom filter at the top of the file. Then, you can write tests for your filters (just follow other tests as examples).
To run your tests, make sure mongo is running locally and pytest is installed and run `pytest tests/test_unit/test_filters`.
Note that a database gets created and deleted for each test. To look at this database, you can temporarily
comment out the `delete_database('mongotest')` line in the test_db fixture in tests/conftest.py (and rerun pytest).

Expose Filter in pickaxe_run.py
-------------------------------
The last step is to expose your filter's settings and options in pickaxe_run.py, the main file used to run pickaxe.
If you open pickaxe_run.py, you will notice different sections for the various built-in filters.
Feel free to create your own section for your custom filter, with whatever parameters it requries as input (in its __init__ method).
Make sure to include a boolean parameter for whether to use your filter (similar to "tani_filter", "tani_sample", "mcs_filter" params).

Once this is done, scroll down to the comment that says "# Apply filters". You'll notice an if statement for each
filter. Under these, add your own if statement to check whether to use your filter with the boolean parameter
you previously specified.

Within this block, instantiate the filter with its required parameters and append it to pk.filters as done for the built-in filters.

That's it! Now, pickaxe will use that filter (as long as its boolean value is set to True) during any expansions.

If you have written tests for your filter and think it could be valuable to the community, feel free to submit a pull request at https://github.com/tyo-nu/MINE-Database to add your filter to the built-in set.