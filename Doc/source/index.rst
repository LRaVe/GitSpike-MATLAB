GitSpike-MATLAB documentation
=============================

.. contents:: Function Navigation
	:depth: 3
	:local:

Common Utilities
================

.. mat:module:: spike_common

add_auxiliary_spikes
--------------------

.. mat:autofunction:: add_auxiliary_spikes

autoMRTS
--------

.. mat:autofunction:: autoMRTS

coincidence_window
------------------

.. mat:autofunction:: coincidence_window

ISI Distance
=============

.. mat:module:: ISI-distance


f_ISI_distance
--------------

.. mat:autofunction:: f_ISI_distance

f_ISI_distance_adaptive_v1
--------------------------

.. mat:autofunction:: f_ISI_distance_adaptive_v1


SPIKE Distance
==============

.. mat:module:: SPIKE_distance

SPIKE_distances_all
-------------------

.. mat:autofunction:: SPIKE_distances_all

.. mat:module:: SPIKE_distance

auxiliary_delta
---------------

.. mat:autofunction:: auxiliary_delta

SPIKE_dist_2x2
--------------

.. mat:autofunction:: SPIKE_dist_2x2

SPIKE_dist_2x2_matlab
---------------------

.. mat:autofunction:: SPIKE_dist_2x2_matlab

SPIKE_dist_2x2_mex
------------------

.. mat:autofunction:: SPIKE_dist_2x2_mex

SPIKE_dist_N
------------

.. mat:autofunction:: SPIKE_dist_N

SPIKE_distances_plot
--------------------

.. mat:autofunction:: SPIKE_distances_plot

SPIKE Order
===========

.. mat:module:: SPIKE_order

order_spikes
------------

.. mat:autofunction:: order_spikes

pairwise_order
--------------

.. mat:autofunction:: pairwise_order

SPIKE Synchronization
=====================

.. mat:module:: SPIKE_synchro

f_adapt_spike_synchro
---------------------

.. mat:autofunction:: f_adapt_spike_synchro

f_spike_synchro
---------------

.. mat:autofunction:: f_spike_synchro

f_spike_synchro_multi
---------------------

.. mat:autofunction:: f_spike_synchro_multi

SPIKE Train Order
=================

.. mat:module:: SPIKE_train_order

compute_spike_train_order_value
-------------------------------

.. mat:autofunction:: compute_spike_train_order_value

order_trains
------------

.. mat:autofunction:: order_trains

pairwise_train_order
--------------------

.. mat:autofunction:: pairwise_train_order

Latency Correction
==================

.. mat:module:: Latency_Correction

f_Cost_matrix
-------------

.. mat:autofunction:: f_Cost_matrix

f_first_diagonal
----------------

.. mat:autofunction:: f_first_diagonal

f_lc_simulated_annealing
------------------------

.. mat:autofunction:: f_lc_simulated_annealing

f_plot_trains_with_correction
-----------------------------

.. mat:autofunction:: f_plot_trains_with_correction

f_row
-----

.. mat:autofunction:: f_row

f_synfire
---------

.. mat:autofunction:: f_synfire

f_TD_matrix
-----------

.. mat:autofunction:: f_TD_matrix

plot_shifts_row
---------------

.. mat:autofunction:: plot_shifts_row

plot_synfire_trains
-------------------

.. mat:autofunction:: plot_synfire_trains

values_to_colors
----------------

.. mat:autofunction:: values_to_colors

Discriminative neuronal subpopulation 
=====================================

This section describes the algorithms implemented to find the most discriminative neuronal subpopulation. Several hypotheses and methods are explained.

.. seealso::
   Kreuz T, et al. *Using spike train distances to identify the most discriminative neuronal subpopulation*, J Neurosci Methods, 2017-2018.

   `Link to see the publication <https://www.thomaskreuz.org/publications/journal-articles#h.j3c8pcr54mx6>`_

Datasets Generation
-------------------

.. mat:module:: ALGO_Simple.Gen_data

generate_SP_dataset
~~~~~~~~~~~~~~~~~~~

.. mat:autofunction:: generate_SP_dataset

generate_LL_dataset
~~~~~~~~~~~~~~~~~~~

.. mat:autofunction:: generate_LL_dataset


Summed Population
-----------------

Raster Plot
~~~~~~~~~~~

.. mat:module:: ALGO_Simple.Summed_Population.Raster_plot

plot_SP_figure
^^^^^^^^^^^^^^

.. mat:autofunction:: plot_SP_figure

pool_neurons
^^^^^^^^^^^^

.. mat:autofunction:: pool_neurons


Performance Matrices
~~~~~~~~~~~~~~~~~~~~

.. mat:module:: ALGO_Simple.Summed_Population.Performance_matrices

build_trials
^^^^^^^^^^^^

.. mat:autofunction:: build_trials

compute_discrimination_performance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: compute_discrimination_performance

compute_population_distance_matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: compute_population_distance_matrix

plot_distance_matrix
^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: plot_distance_matrix

Algorithms
~~~~~~~~~~

.. mat:module:: ALGO_Simple.Summed_Population.Algorithms

evaluate_population
^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: evaluate_population

f_brute_force_V2
^^^^^^^^^^^^^^^^

.. mat:autofunction:: f_brute_force_V2

f_bottom_up_V2
^^^^^^^^^^^^^^

.. mat:autofunction:: f_bottom_up_V2

top_down_gradient
^^^^^^^^^^^^^^^^^

.. mat:autofunction:: top_down_gradient

plot_top_down_gradient
^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: plot_top_down_gradient

Simulated Annealing
~~~~~~~~~~~~~~~~~~~

.. mat:module:: ALGO_Simple.Summed_Population.Simulated_Annealing

evaluate_population_cached
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: evaluate_population_cached

initialize_temperature
^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: initialize_temperature

metropolis_acceptance
^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: metropolis_acceptance


random_neighbor
^^^^^^^^^^^^^^^

.. mat:autofunction:: random_neighbor

simulated_annealing
^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: simulated_annealing

plot_simulated_annealing
^^^^^^^^^^^^^^^^^^^^^^^^

.. mat:autofunction:: plot_simulated_annealing

Labeled Line 
-------------

.. mat:module:: ALGO_Simple.Labeled_Line

evaluate_LL_population
~~~~~~~~~~~~~~~~~~~~~~

.. mat:autofunction:: evaluate_LL_population

plot_LL_results
~~~~~~~~~~~~~~~

.. mat:autofunction:: plot_LL_results

plot_LL_figure4
~~~~~~~~~~~~~~~

.. mat:autofunction:: plot_LL_figure4

Managing the data
-----------------

.. mat:module:: ALGO_Simple

export_spikes_to_txt
~~~~~~~~~~~~~~~~~~~~

.. mat:autofunction:: export_spikes_to_txt

import_spikes_from_txt
~~~~~~~~~~~~~~~~~~~~~~

.. mat:autofunction:: import_spikes_from_txt