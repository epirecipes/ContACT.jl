module ContACT

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.UndirectedWiringDiagrams: AbstractUWD
using Catlab.Programs: @relation
using DataFrames
using LinearAlgebra
using Random

# Core types
export AbstractPartition, IntervalPartition, CategoricalPartition, ProductPartition,
       AgePartition, ContactSurvey, ContactMatrix, SourceStratifiedContactMatrix,
       UnitSemantics, MeanContacts, ContactCounts, PerCapitaRate,
       dimension, group_labels, age_limits, age_labels, n_groups, matrix, population,
       target_partition, source_partition, n_target_groups, n_source_groups,
       target_group_labels, source_group_labels, source_total_contacts

# Survey operations
export filter_survey, subset_survey

# Functor: Survey → ContactMatrix
export compute_matrix, compute_source_stratified_matrix,
       coarsen_sources, align_source_stratified_matrix

# Coarsening & refinement
export coarsen, refine, PartitionMap, AgeMap
export ActivityRefinement, ActivityMixingKernel,
       AssortativeMixing, DisassortativeMixing, ProportionateMixing,
       activity_partition, activity_mixing_plan, activity_refine
export GeneralizedLift, GeneralizedMixingKernel,
       RandomMixing, BlockMixing, AssortativeDimensionMixing,
       product_population, generalize, generalized_lift
export ConstrainedGeneralizedLift, full_partition, intermediate_matrix,
       source_map, structural_zeros, constrained_generalize,
       BlockAssortativityParams, ParameterizedConstrainedLift,
       is_feasible, sample_constrained_lifts

# Composition
export compose_matrices

# Stratification
export stratify

# Symmetrisation
export symmetrise

# Operators (type \name<TAB> in REPL)
export ⊕, ⊗, ↓, ↑, ⤊, ⊠, ▷, ↔, ρ, ×, RefinementPrior

# Utilities
export to_per_capita, to_counts, spectral_radius,
       next_generation_matrix, basic_reproduction_number,
       calibrate_transmissibility, R0, R₀,
       marginal_matrix, assortativity_index, type_reproduction_number,
       control_threshold, control_effort

# ACSet schemas & categorical machinery
export SchContactSurvey, ContactSurveyACSet,
       SchContactMatrix, LabelledContactMatrix,
       OpenContactMatrixOb, OpenContactMatrix,
       migrate_coarsen,
       ContactSharer, compose_uwd

include("types.jl")
include("survey.jl")
include("functor.jl")
include("coarsening.jl")
include("source_stratified.jl")
include("constrained_lift.jl")
include("refinement.jl")
include("composition.jl")
include("stratification.jl")
include("symmetrise.jl")
include("activity_refinement.jl")
include("generalized_lift.jl")
include("operators.jl")
include("schemas.jl")
include("utils.jl")

end # module
