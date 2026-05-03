module ContACT

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.WiringDiagrams
using Catlab.WiringDiagrams.UndirectedWiringDiagrams: AbstractUWD
using Catlab.Programs: @relation
using DataFrames
using LinearAlgebra

# Core types
export AbstractPartition, IntervalPartition, CategoricalPartition, ProductPartition,
       AgePartition, ContactSurvey, ContactMatrix,
       UnitSemantics, MeanContacts, ContactCounts, PerCapitaRate,
       dimension, group_labels, age_limits, age_labels, n_groups, matrix, population

# Survey operations
export filter_survey, subset_survey

# Functor: Survey → ContactMatrix
export compute_matrix

# Coarsening & refinement
export coarsen, refine, PartitionMap, AgeMap
export ActivityRefinement, ActivityMixingKernel,
       AssortativeMixing, DisassortativeMixing, ProportionateMixing,
       activity_partition, activity_mixing_plan, activity_refine

# Composition
export compose_matrices

# Stratification
export stratify

# Symmetrisation
export symmetrise

# Operators (type \name<TAB> in REPL)
export ⊕, ⊗, ↓, ↑, ⤊, ▷, ↔, ρ, ×, RefinementPrior

# Utilities
export to_per_capita, to_counts, spectral_radius

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
include("refinement.jl")
include("composition.jl")
include("stratification.jl")
include("symmetrise.jl")
include("activity_refinement.jl")
include("operators.jl")
include("schemas.jl")
include("utils.jl")

end # module
