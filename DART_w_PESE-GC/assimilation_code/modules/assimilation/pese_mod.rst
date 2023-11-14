MODULE pese_mod
===============

Module author: Man-Yau ("Joseph") Chan

Overview
--------

The code in this module increases the number of ensemble members through creating virtual members. These virtual members
will not be written to disk -- i.e., virtual members only exist in the computer's Random Access Memory (RAM) when running
DART's ```filter``` executable. 

This virtual member creation procedure is called "Probit-space Ensemble Size Expansion for Gaussian Copulas" or "PESE-GC"
("peace gee see") for short. PESE-GC is documented in Man-Yau Chan's manuscript on PESE-GC.




Caveats for this version of ``pese_mod``
----------------------------------------

This version of PESE-GC has only been tested with the Normal Distribution and the Bounded Normal Rank Histogram (BNRH)
Distribution. 

Furthermore, because PESE-GC has only been tested on the Lorenz 1996 model, for simplicity, it is assumed that every 
model variable has the same family of marginal distributions. It is also assumed that every prior observation follows
the same family of marginal distributions. In other words, this verion of ``pese_mod`` only supports having two different
distribution families: one family for model variables and another family for observation variables.

To be clear, the theoretical formulation of PESE-GC has no "same family" assumption built in! PESE-GC can support any 
combination of marginal distributions! The current version of ``pese_mod`` can be (easily) upgraded to use multiple 
marginal distribution families for the model variables (e.g., using Normal for zonal wind, BNRH for specific humidity,
Gamma for cloud water mixing ratios at the same time).





Other modules used
------------------

::

   types_mod 
   io_filenames_mod 
   random_seq_mod 
   utilities_mod 
   ensemble_manager_mod 
   probit_transform_mod 
   algorithm_info_mod 
   location_mod 
   distribution_params_mod 
   assim_model_mod 



Public interfaces
-----------------

====================== ===================================
*use pese_mod, only :* init_pese_info
\                      generate_pese_resampling_weights
\                      pese_type
\                      pese_modify_file_info_io_flags
\                      apply_pese
\                      pese_modify_all_file_info_io_flags
====================== ===================================

| 

.. container:: routine

   *call init_pese_info( expanded_ens_size, dynamic_ens_size )*
   ::

      integer,                intent(in)        :: expanded_ens_size
      integer,                intent(in)        :: dynamic_ens_size
|

.. container:: indent1

   Initializes and allocates variables needed for PESE-GC. This includes populating the structured type `pese_type`
   and computing the coefficients needed for resampling.

   +-----------------------+-----------------------------------------------------------+
   | ``expanded_ens_size`` | Number of ensemble members after applying PESE-GC         |
   +-----------------------+-----------------------------------------------------------+
   | ``dynamic_ens_size``  | Number of ensemble members before applying PESE-GC (i.e., |
   |                       | number of forecast ensemble members.)                     |
   +-----------------------+-----------------------------------------------------------+

| 

.. container:: routine

   *call generate_pese_resampling_weights( )*


.. container:: indent1

   Generates Gaussian resampling coefficients. These are used to construct virtual members in probit space.
   Note that this subroutine reads `pese_mod`'s structured type variable `pese_info` to figure out how many
   forecast ensemble members exist and how many virtual ensemble members to create.
   

| 
.. container:: routine

   *call pese_modify_file_info_io_flags( file_info, ENS_MEM_START, ens_size_dynamic, ens_size_expanded )*
   ::
      type(file_info_type),         intent(inout)  :: file_info
      integer,                      intent(in)     :: ENS_MEM_START
      integer,                      intent(in)     :: ens_size_dynamic
      integer,                      intent(in)     :: ens_size_expanded


.. container:: indent1

   Modifies DART's IO flags to prevent the virtual ensemble members from being written out.

   +-----------------------+-----------------------------------------------------------+
   | ``file_info``         | Structured type containing settings about what to write   |
   |                       | for a particular stage of DART.                           |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+
   | ``ENS_MEM_START``     | Copy number of the first ensemble member in `ens_handle`  |
   +-----------------------+-----------------------------------------------------------+
   | ``expanded_ens_size`` | Number of ensemble members after applying PESE-GC         |
   +-----------------------+-----------------------------------------------------------+
   | ``dynamic_ens_size``  | Number of ensemble members before applying PESE-GC (i.e., |
   |                       | number of forecast ensemble members.)                     |
   +-----------------------+-----------------------------------------------------------+


| 
.. container:: routine

   *call apply_pese( ens_handle, print_pese_debug )*
   ::
      type(ensemble_type),          intent(inout)  :: ens_handle
      logical,                      intent(in)     :: print_pese_debug


.. container:: indent1

   This subroutine generate virtual ensemble members from the forecast ensemble members.

   +-----------------------+-----------------------------------------------------------+
   | ``ens_handle``        | Handle for the model state variable ensemble              |
   +-----------------------+-----------------------------------------------------------+
   | ``print_pese_debug``  | Logical controlling whether to print debug messages from  |
   |                       | ``pese_mod``                                              |
   +-----------------------+-----------------------------------------------------------+


| 
.. container:: routine

   *call pese_modify_all_file_info_io_flags( dynamic_ens_size, ens_size, ENS_MEM_START, file_info_input, 
      file_info_forecast, file_info_preassim, file_info_postassim, file_info_analysis, file_info_output )*
   ::
      integer,                            intent(in)     :: dynamic_ens_size
      integer,                            intent(in)     :: ens_size
      integer,                            intent(in)     :: ENS_MEM_START
      type( file_info_type ),             intent(inout)  :: file_info_input
      type( file_info_type ),             intent(inout)  :: file_info_forecast
      type( file_info_type ),             intent(inout)  :: file_info_preassim
      type( file_info_type ),             intent(inout)  :: file_info_postassim
      type( file_info_type ),             intent(inout)  :: file_info_analysis
      type( file_info_type ),             intent(inout)  :: file_info_output

.. container:: indent1

   Modifies DART's IO flags to prevent the virtual ensemble members from being written out.
   This subroutine calls ``pese_modify_file_info_io_flags`` on the file handle of every stage
   in DART's ``filter`` (``input``, ``forecast``, ``preassim``, ``analysis``, ``output``).

   +-----------------------+-----------------------------------------------------------+
   | ``dynamic_ens_size``  | Number of ensemble members before applying PESE-GC (i.e., |
   |                       | number of forecast ensemble members.)                     |
   +-----------------------+-----------------------------------------------------------+
   | ``ens_size``          | Number of ensemble members after applying PESE-GC         |
   +-----------------------+-----------------------------------------------------------+
   | ``ENS_MEM_START``     | Copy number of the first ensemble member in `ens_handle`  |
   +-----------------------+-----------------------------------------------------------+
   | ``file_info_input``   | Structured type containing settings about what to write   |
   |                       | for the ``input`` stage of DART.                          |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+
   | ``file_info_forecast``| Structured type containing settings about what to write   |
   |                       | for the ``forecast`` stage of DART.                       |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+
   | ``file_info_preassim``| Structured type containing settings about what to write   |
   |                       | for the ``preassim`` stage of DART.                       |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+
   | ``file_info_analysis``| Structured type containing settings about what to write   |
   |                       | for the ``analysis`` stage of DART.                       |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+
   | ``file_info_output``  | Structured type containing settings about what to write   |
   |                       | for the ``output`` stage of DART.                         |
   |                       | For more information, see ``../io/io_filenames_mod.f90``. |
   +-----------------------+-----------------------------------------------------------+

| 

Namelist
--------

This module is controlled by two namelists: ``filter_nml`` and ``algorithm_info_mod_nml``

+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&filter_nml :: use_pese``                               | Logical controlling whether PESE-GC is used.              |
|                                                           | PESE-GC is used if ``use_pese = .true.``, and is unused   |
|                                                           | if ``use_pese = .false.``.                                |
+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&filter_nml :: ens_size_expanded``                      | Size of the expanded ensemble (i.e., number of forecast   |
|                                                           | members plus number of virtual members)                   |
+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&filter_nml :: print_pese_debug``                       | Logical controlling whether to print debug messages from  |
|                                                           | ``pese_mod``                                              |
+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&algorithm_info_mod_nml :: filt_type_override``         | Just set this to the same value as                        |
|                                                           | ``assim_tools_mod :: filter_kind``                        |
+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&algorithm_info_mod_nml :: dist_type_override_model``   | Marginal distribution for model variables to use with     |
|                                                           | PESE-GC                                                   |
+-----------------------------------------------------------+-----------------------------------------------------------+
| ``&algorithm_info_mod_nml :: dist_type_override_obs``     | Marginal distribution for observation variables to use    |
|                                                           | with PESE-GC                                              |
+-----------------------------------------------------------+-----------------------------------------------------------+



References
----------

Chan, M.-Y. (in prep): Improving Ensemble Data Assimilation through Probit-space Ensemble Size Expansion for Gaussian
Copulas (PESE-GC). Will be submitted to Nonlinear Processes in Geophysics.


Private components
------------------

N/A
