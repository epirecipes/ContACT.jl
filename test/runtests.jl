using Test
using ContACT
using CSV
using DataFrames
using LinearAlgebra
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Programs: @relation

# Resolve ambiguity: use ContACT's operators explicitly
import ContACT: ⊕, ⊗, ↓, ↑, ⤊, ⊠, ▷, ↔, ρ, ×

# ---------------------------------------------------------------------------
# Test data: synthetic survey
# ---------------------------------------------------------------------------

function make_test_survey()
    # 10 participants, ages 5–70
    participants = DataFrame(
        part_id = 1:10,
        part_age = [5.0, 8.0, 15.0, 25.0, 30.0, 45.0, 50.0, 60.0, 70.0, 72.0],
        country = ["UK","UK","UK","UK","UK","DE","DE","DE","DE","DE"]
    )
    # Each participant makes some contacts
    contacts = DataFrame(
        part_id = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10],
        cnt_age = [6.0,25.0, 7.0,30.0, 14.0,40.0, 22.0,55.0, 28.0,65.0,
                   44.0,10.0, 48.0,20.0, 62.0,5.0, 68.0,35.0, 71.0,15.0]
    )
    ContactSurvey(participants, contacts)
end

# ---------------------------------------------------------------------------
@testset "ContACT.jl" begin
# ---------------------------------------------------------------------------

@testset "AgePartition" begin
    p = AgePartition([0, 18, 65])
    @test n_groups(p) == 3
    @test age_limits(p) == [0.0, 18.0, 65.0]
    @test age_labels(p) == ["[0,18)", "[18,65)", "65+"]
    @test age_limits(AgePartition([65, 0, 18])) == [0.0, 18.0, 65.0]

    # Custom labels
    p2 = AgePartition([0, 5, 18]; labels=["child", "youth", "adult"])
    @test age_labels(p2) == ["child", "youth", "adult"]
    @test_throws ArgumentError AgePartition([18, 0, 65]; labels=["adult", "child", "senior"])

    # Validation
    @test_throws ArgumentError AgePartition(Int[])
    @test_throws ArgumentError AgePartition([5, 5, 10])
    @test_throws ArgumentError AgePartition([-1, 0, 18])
    @test_throws ArgumentError AgePartition([0, Inf])
end

@testset "ContactSurvey" begin
    survey = make_test_survey()
    @test nrow(survey.participants) == 10
    @test nrow(survey.contacts) == 20

    # Filter by country
    uk = filter_survey(survey; country="UK")
    @test nrow(uk.participants) == 5
    @test all(uk.participants.country .== "UK")
    @test all(row -> row.part_id in Set(uk.participants.part_id), eachrow(uk.contacts))

    # Filter by function
    adults = filter_survey(survey; part_age = a -> a >= 18)
    @test all(adults.participants.part_age .>= 18)

    # Subset
    sub = subset_survey(survey, [1, 3, 5])
    @test nrow(sub.participants) == 3
end

@testset "General partitions" begin
    participants = DataFrame(
        part_id = 1:4,
        part_age = [8.0, 30.0, 12.0, 50.0],
        part_sex = ["F", "M", "F", "M"],
        part_region = ["North", "North", "South", "South"],
    )
    contacts = DataFrame(
        part_id = [1, 1, 2, 3, 4],
        cnt_age = [35.0, 7.0, 31.0, 10.0, 9.0],
        cnt_sex = ["M", "F", "M", "F", "F"],
        cnt_region = ["North", "South", "North", "South", "North"],
    )
    survey = ContactSurvey(participants, contacts)

    sex = CategoricalPartition(:sex;
        participant_col=:part_sex,
        contact_col=:cnt_sex,
        levels=["F", "M"],
        labels=["female", "male"],
    )
    @test dimension(sex) == :sex
    @test group_labels(sex) == ["female", "male"]

    cm_sex = survey ▷ sex
    @test matrix(cm_sex) ≈ [1.0 0.5; 0.5 0.5]
    @test population(cm_sex) == [2.0, 2.0]
    @test group_labels(cm_sex) == ["female", "male"]
    @test_throws ArgumentError age_limits(cm_sex)

    all_sex = CategoricalPartition(:sex;
        participant_col=:part_sex,
        contact_col=:cnt_sex,
        levels=["all"],
    )
    f_all = PartitionMap(sex, all_sex, Dict("F" => "all", "M" => "all"))
    cm_all = cm_sex ↓ f_all
    @test n_groups(cm_all) == 1
    @test population(cm_all) == [4.0]
    @test matrix(cm_all) ≈ [1.25;;]

    region = CategoricalPartition(:region;
        participant_col=:part_region,
        contact_col=:cnt_region,
        levels=["North", "South"],
    )
    sex_region = sex × region
    @test dimension(sex_region) == (:sex, :region)
    @test group_labels(sex_region) == ["female:North", "female:South", "male:North", "male:South"]

    cm_product = survey ▷ sex_region
    @test n_groups(cm_product) == 4
    @test size(matrix(cm_product)) == (4, 4)
    @test population(cm_product) == [1.0, 1.0, 1.0, 1.0]
    @test_throws ArgumentError age_limits(cm_product)

    @test_throws ArgumentError survey ▷ CategoricalPartition(:occupation;
        participant_col=:part_occupation,
        contact_col=:cnt_occupation,
        levels=["worker"],
    )
end

@testset "Source-stratified matrices" begin
    participants = DataFrame(
        part_id = 1:4,
        part_age = [10.0, 10.0, 30.0, 30.0],
        part_sep = ["low", "high", "low", "high"],
        part_edu = ["low", "low", "high", "high"],
        weight = [2.0, 1.0, 1.0, 1.0],
    )
    contacts = DataFrame(
        part_id = [1, 1, 2, 3, 3, 4, 4],
        cnt_age = [10.0, 30.0, 10.0, 30.0, 30.0, 10.0, 30.0],
    )
    survey = ContactSurvey(participants, contacts)

    age = AgePartition([0, 18]; labels=["child", "adult"])
    sep = CategoricalPartition(:sep;
        participant_col=:part_sep,
        contact_col=:cnt_sep,
        levels=["low", "high"],
    )
    edu = CategoricalPartition(:edu;
        participant_col=:part_edu,
        contact_col=:cnt_edu,
        levels=["low", "high"],
    )
    source = age × sep × edu

    partial = compute_source_stratified_matrix(survey, age, source)
    @test partial isa SourceStratifiedContactMatrix
    @test size(partial) == (2, 8)
    @test target_partition(partial) == age
    @test source_partition(partial) == source
    @test n_target_groups(partial) == 2
    @test n_source_groups(partial) == 8
    @test target_group_labels(partial) == ["child", "adult"]
    @test source_group_labels(partial)[1:4] == [
        "child:low:low", "child:low:high", "child:high:low", "child:high:high"
    ]
    @test population(partial) ≈ [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    @test matrix(partial) ≈ [
        1.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0;
        1.0 0.0 0.0 0.0 0.0 2.0 0.0 1.0
    ]
    @test source_total_contacts(partial) ≈ matrix(partial) * Diagonal(population(partial))

    weighted = compute_source_stratified_matrix(survey, age, source; weights=:weight)
    @test population(weighted) ≈ [2.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
    @test matrix(weighted) ≈ matrix(partial)

    provided_pop = [10.0, 0.0, 8.0, 0.0, 0.0, 12.0, 0.0, 9.0]
    with_pop = compute_source_stratified_matrix(survey, age, source; population=provided_pop)
    @test population(with_pop) == provided_pop

    @test_throws DimensionMismatch SourceStratifiedContactMatrix(
        zeros(2, 7), age, source, provided_pop
    )
    @test_throws DimensionMismatch SourceStratifiedContactMatrix(
        zeros(2, 8), age, source, provided_pop[1:7]
    )
    @test_throws ArgumentError compute_source_stratified_matrix(
        survey, age, source; weights=:missing_weight
    )

    source_to_age = PartitionMap(source, age)
    age_source = coarsen_sources(partial, source_to_age)
    @test target_partition(age_source) == age
    @test source_partition(age_source) == age
    @test population(age_source) ≈ [2.0, 2.0]
    @test matrix(age_source) ≈ [1.0 0.5; 0.5 1.5]

    base = ContactMatrix([2.0 1.0; 1.0 3.0], age, [2.0, 2.0])
    aligned = align_source_stratified_matrix(partial, base, source_to_age)
    @test matrix(coarsen_sources(aligned, source_to_age)) ≈ matrix(base)
    @test source_total_contacts(aligned) ≈ 2 .* source_total_contacts(partial)

    @test_throws ArgumentError align_source_stratified_matrix(
        partial,
        ContactMatrix([2.0 1.0; 1.0 3.0], age, [3.0, 1.0]),
        source_to_age,
    )
    @test_throws ArgumentError align_source_stratified_matrix(
        partial,
        ContactMatrix([2.0 2.0; 0.0 3.0], age, [2.0, 2.0]),
        source_to_age,
    )
    @test_throws ArgumentError align_source_stratified_matrix(
        SourceStratifiedContactMatrix(zeros(2, 8), age, source, population(partial)),
        base,
        source_to_age,
    )

    constrained = ConstrainedGeneralizedLift(aligned; source_map=source_to_age)
    @test full_partition(constrained) == source
    @test intermediate_matrix(constrained) == aligned
    @test source_map(constrained) == source_to_age
    @test structural_zeros(constrained) == (population(aligned) .== 0)

    @test_throws DimensionMismatch ConstrainedGeneralizedLift(
        aligned; source_map=source_to_age, structural_zeros=[false, true]
    )
    @test_throws ArgumentError ConstrainedGeneralizedLift(
        aligned; source_map=source_to_age, structural_zeros=trues(8)
    )

    full = base ⊠ constrained
    @test full isa ContactMatrix
    @test full.partition == source
    @test population(full) ≈ population(aligned)
    @test matrix(full ↓ age) ≈ matrix(base)
    @test matrix(base ↑ constrained) ≈ matrix(full)
    @test matrix(constrained_generalize(base, constrained)) ≈ matrix(full)

    full_counts = matrix(full) * Diagonal(population(full))
    @test full_counts ≈ transpose(full_counts)
    aligned_counts = source_total_contacts(aligned)
    fmap = collect(source_to_age.mapping)
    reconstructed_intermediate = zeros(size(aligned_counts))
    for source_group in axes(full_counts, 2)
        for target_group in axes(full_counts, 1)
            reconstructed_intermediate[fmap[target_group], source_group] +=
                full_counts[target_group, source_group]
        end
    end
    @test reconstructed_intermediate ≈ aligned_counts

    for group in findall(structural_zeros(constrained))
        @test all(iszero, full_counts[:, group])
        @test all(iszero, full_counts[group, :])
    end
end

@testset "Parameterized constrained lift (q-parameters)" begin
    # Setup: 2 age groups × 2 SEP groups
    age = AgePartition([0, 18]; labels=["child", "adult"])
    sep = CategoricalPartition(:sep;
        participant_col=:part_sep,
        contact_col=:cnt_sep,
        levels=["low", "high"],
    )
    prod = age × sep
    # 4 groups: (child,low), (child,high), (adult,low), (adult,high)
    prod_pop = [100.0, 100.0, 150.0, 150.0]

    # Build a source-stratified intermediate: 2 target (age) × 4 source (age×sep)
    # Mean contacts reported by each source group into each target age group
    interm_M = [3.0 2.0 1.0 1.5;   # child contacted by (child,low), (child,high), (adult,low), (adult,high)
                1.5 1.0 2.5 2.0]   # adult contacted by ...
    interm = SourceStratifiedContactMatrix(interm_M, age, prod, prod_pop)

    # Reciprocal base matrix over age
    base_pop = [200.0, 300.0]
    base_M = [5.0 2.5; 2.5 4.5]
    base_counts = base_M * Diagonal(base_pop)
    # Ensure reciprocal: symmetrize total contacts
    base_counts_sym = (base_counts + base_counts') / 2
    base_M_sym = base_counts_sym * Diagonal(1.0 ./ base_pop)
    base = ContactMatrix(base_M_sym, age, base_pop)

    source_to_age = PartitionMap(prod, age)
    spec = ConstrainedGeneralizedLift(interm; source_map=source_to_age)

    @testset "q=0 matches proportionate solver" begin
        params = BlockAssortativityParams(q=Dict(:sep => 0.0))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        result_q0 = constrained_generalize(base, pspec)
        result_prop = constrained_generalize(base, spec)
        @test matrix(result_q0) ≈ matrix(result_prop)
    end

    @testset "assortative q > 0 increases same-SEP contacts" begin
        params = BlockAssortativityParams(q=Dict(:sep => 0.5))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        result_assort = constrained_generalize(base, pspec)
        result_prop = constrained_generalize(base, spec)
        # Assortativity index for SEP should be higher
        ai_assort = assortativity_index(result_assort, :sep)
        ai_prop = assortativity_index(result_prop, :sep)
        @test ai_assort > ai_prop
    end

    @testset "disassortative q < 0 decreases same-SEP contacts" begin
        params = BlockAssortativityParams(q=Dict(:sep => -0.5))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        result_dis = constrained_generalize(base, pspec)
        result_prop = constrained_generalize(base, spec)
        ai_dis = assortativity_index(result_dis, :sep)
        ai_prop = assortativity_index(result_prop, :sep)
        @test ai_dis < ai_prop
    end

    @testset "coarsening invariant preserved" begin
        params = BlockAssortativityParams(q=Dict(:sep => 0.3))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        result = constrained_generalize(base, pspec)
        @test matrix(result ↓ age) ≈ matrix(base) atol=1e-10
    end

    @testset "reciprocity preserved" begin
        params = BlockAssortativityParams(q=Dict(:sep => 0.4))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        result = constrained_generalize(base, pspec)
        C = matrix(result) * Diagonal(population(result))
        @test C ≈ C' atol=1e-10
    end

    @testset "operator dispatch" begin
        params = BlockAssortativityParams(q=Dict(:sep => 0.2))
        pspec = ParameterizedConstrainedLift(spec; default_params=params)
        @test matrix(base ⊠ pspec) ≈ matrix(constrained_generalize(base, pspec))
        @test matrix(base ↑ pspec) ≈ matrix(constrained_generalize(base, pspec))
    end

    @testset "is_feasible" begin
        feasible_params = BlockAssortativityParams(q=Dict(:sep => 0.3))
        pspec_f = ParameterizedConstrainedLift(spec; default_params=feasible_params)
        @test is_feasible(base, pspec_f) == true

        # q is constrained to [-1, 1]. Every interior value is feasible here —
        # even q=0.99 leaves a strictly positive minimum entry — while both
        # endpoints fail, because the transport cannot then be balanced against
        # the base marginals.
        near_params = BlockAssortativityParams(q=Dict(:sep => 0.99))
        pspec_n = ParameterizedConstrainedLift(spec; default_params=near_params)
        @test is_feasible(base, pspec_n) == true
        @test minimum(matrix(constrained_generalize(base, pspec_n))) > 0

        for q_end in (1.0, -1.0)
            pspec_e = ParameterizedConstrainedLift(spec;
                default_params=BlockAssortativityParams(q=Dict(:sep => q_end)))
            @test is_feasible(base, pspec_e) == false
            @test_throws ArgumentError constrained_generalize(base, pspec_e)
        end
    end

    @testset "sample_constrained_lifts" begin
        using Random
        rng = Random.MersenneTwister(42)
        samples = sample_constrained_lifts(base, spec, 5;
            dimensions=[:sep], bounds=(-0.5, 0.5), rng=rng)
        # Every |q| < 1 is feasible for this base/spec (see "is_feasible" above), and
        # (-0.5, 0.5) is strictly interior, so every draw is accepted regardless of seed.
        @test length(samples) == 5
        for (params, cm) in samples
            @test cm isa ContactMatrix
            @test all(matrix(cm) .>= -1e-10)
            # Coarsening invariant
            @test matrix(cm ↓ age) ≈ matrix(base) atol=1e-10
        end
    end

    @testset "per-block params" begin
        block_params = Dict(
            (1, 1) => BlockAssortativityParams(q=Dict(:sep => 0.5)),
            (2, 2) => BlockAssortativityParams(q=Dict(:sep => -0.3)),
        )
        pspec = ParameterizedConstrainedLift(spec;
            block_params=block_params,
            default_params=BlockAssortativityParams(q=Dict(:sep => 0.0)))
        result = constrained_generalize(base, pspec)
        @test result isa ContactMatrix
        @test matrix(result ↓ age) ≈ matrix(base) atol=1e-10
    end

    @testset "QParameterSpace" begin
        space = QParameterSpace(base, spec; dimensions=[:sep])
        # 2 age groups → 3 blocks: (1,1), (1,2), (2,2); 1 dimension each
        @test space.n_params == 3
        @test length(space.block_keys) == 3
        @test space.dimensions == [:sep]

        # Round-trip vector ↔ block_params
        θ = [0.1, 0.2, -0.3]
        bp = ContACT._vector_to_block_params(space, θ)
        @test bp[(1,1)].q[:sep] ≈ 0.1
        @test bp[(1,2)].q[:sep] ≈ 0.2
        @test bp[(2,2)].q[:sep] ≈ -0.3

        θ_back = ContACT._block_params_to_vector(space, bp)
        @test θ_back ≈ θ
    end

    @testset "sample_perblock_lifts" begin
        using Random
        rng = Random.MersenneTwister(123)
        samples = sample_perblock_lifts(base, spec, 3;
            dimensions=[:sep], bounds=(-0.4, 0.4), rng=rng)
        # Same reasoning as sample_constrained_lifts above: (-0.4, 0.4) is strictly
        # interior to the feasible range, so every draw is accepted regardless of seed.
        @test length(samples) == 3
        for (bp, cm) in samples
            @test cm isa ContactMatrix
            @test all(matrix(cm) .>= -1e-10)
            @test matrix(cm ↓ age) ≈ matrix(base) atol=1e-10
            # Verify it's truly per-block: different blocks can have different q
            @test bp isa Dict{Tuple{Int,Int},BlockAssortativityParams}
        end
    end

    @testset "mcmc_constrained_lifts" begin
        using Random
        rng = Random.MersenneTwister(99)

        # Uniform (flat prior) MCMC
        result = mcmc_constrained_lifts(base, spec, 10;
            dimensions=[:sep], bounds=(-0.4, 0.4),
            proposal_scale=0.05, burnin=20, thin=1, rng=rng)
        @test result isa MCMCResult
        @test length(result.chain) == 10
        @test length(result.matrices) == 10
        @test length(result.log_densities) == 10
        @test 0.0 <= result.acceptance_rate <= 1.0
        @test result.space.n_params == 3

        # All samples preserve invariants
        for cm in result.matrices
            @test all(matrix(cm) .>= -1e-10)
            @test matrix(cm ↓ age) ≈ matrix(base) atol=1e-10
        end

        # Omitting log_density means a flat target, recorded as zero throughout.
        @test all(iszero, result.log_densities)

        # A chain run against a log-density records, for each sample, the density
        # of the state it actually stored. Both are exact: the recorded value is
        # the one computed during the step, and the stored parameters regenerate
        # the stored matrix.
        log_dens = (cm, _) -> assortativity_index(cm, :sep)
        result_targeted = mcmc_constrained_lifts(base, spec, 10;
            dimensions=[:sep], bounds=(-0.4, 0.4),
            log_density=log_dens,
            proposal_scale=0.05, burnin=20, thin=1, rng=Random.MersenneTwister(77))
        @test length(result_targeted.chain) == 10
        for k in 1:10
            cm, bp = result_targeted.matrices[k], result_targeted.chain[k]
            @test result_targeted.log_densities[k] == log_dens(cm, bp)
            regenerated = constrained_generalize(base,
                ParameterizedConstrainedLift(spec; block_params=bp))
            @test matrix(regenerated) == matrix(cm)
            @test all(matrix(cm) .>= -1e-10)
            @test matrix(cm ↓ age) ≈ matrix(base) atol=1e-10
        end

        # A constant density gives every proposal an acceptance ratio of 1, which
        # is what omitting log_density does. On a shared seed the two chains are
        # therefore identical step for step.
        for c in (0.0, 3.7, -12.5)
            r_flat = mcmc_constrained_lifts(base, spec, 10;
                dimensions=[:sep], bounds=(-0.4, 0.4),
                proposal_scale=0.05, burnin=20, thin=1, rng=Random.MersenneTwister(4))
            r_const = mcmc_constrained_lifts(base, spec, 10;
                dimensions=[:sep], bounds=(-0.4, 0.4), log_density=(cm, bp) -> c,
                proposal_scale=0.05, burnin=20, thin=1, rng=Random.MersenneTwister(4))
            @test all(matrix(r_const.matrices[k]) == matrix(r_flat.matrices[k]) for k in 1:10)
            @test r_const.acceptance_rate == r_flat.acceptance_rate
            @test all(==(c), r_const.log_densities)
        end
    end
end

@testset "ContactMatrix construction" begin
    p = AgePartition([0, 18, 65])
    M = [2.0 1.0 0.5; 1.5 3.0 1.0; 0.3 0.8 1.5]
    pop = [1000.0, 3000.0, 500.0]
    cm = ContactMatrix(M, p, pop)

    @test n_groups(cm) == 3
    @test matrix(cm) == M
    @test population(cm) == pop
    @test cm.semantics isa MeanContacts

    # Dimension mismatch
    @test_throws DimensionMismatch ContactMatrix(M, AgePartition([0, 10]), pop)
    @test_throws DimensionMismatch ContactMatrix(M, p, [1.0, 2.0])

    # Value validation
    @test_throws ArgumentError ContactMatrix([-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], p, pop)
    @test_throws ArgumentError ContactMatrix([NaN 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], p, pop)
    @test_throws ArgumentError ContactMatrix(M, p, [1000.0, -1.0, 500.0])
    @test_throws ArgumentError ContactMatrix(M, p, [1000.0, Inf, 500.0])
end

@testset "compute_matrix (functor)" begin
    survey = make_test_survey()
    partition = AgePartition([0, 18, 65])
    cm = compute_matrix(survey, partition)

    @test n_groups(cm) == 3
    @test size(matrix(cm)) == (3, 3)
    # All entries should be non-negative
    @test all(matrix(cm) .>= 0)
    # Should have some contacts in each group
    @test sum(matrix(cm)) > 0

    weighted_participants = DataFrame(part_id=[1, 2], part_age=[5.0, 25.0], weight=[2.0, 1.0])
    weighted_contacts = DataFrame(part_id=[1, 1, 2], cnt_age=[5.0, 25.0, 25.0])
    weighted = ContactSurvey(weighted_participants, weighted_contacts)
    weighted_cm = compute_matrix(weighted, AgePartition([0, 18]); weights=:weight)
    @test matrix(weighted_cm) ≈ [1.0 0.0; 1.0 1.0]
    @test_throws ArgumentError compute_matrix(weighted, AgePartition([0, 18]); weights=:missing_weight)
    weighted_participants.bad_weight = [1.0, -1.0]
    @test_throws ArgumentError compute_matrix(
        ContactSurvey(weighted_participants, weighted_contacts),
        AgePartition([0, 18]);
        weights=:bad_weight,
    )
end

@testset "dropped contacts are excluded, not miscounted" begin
    partition = AgePartition([0, 18, 65])
    # participant 3's age and contact row 2's age both fall below the partition's
    # first limit, so assign_group returns nothing for each (src/types.jl:214-221).
    participants = DataFrame(part_id=[1, 2, 3], part_age=[10.0, 30.0, -5.0])
    contacts = DataFrame(part_id=[1, 1, 2, 3], cnt_age=[10.0, -1.0, 30.0, 30.0])
    survey = ContactSurvey(participants, contacts)

    cm = @test_logs(
        (:warn, r"dropped 2 contact\(s\): 1 with missing/unmapped contact group; 1 with missing/unmapped participant group"),
        compute_matrix(survey, partition))
    @test population(cm) == [1.0, 1.0, 0.0]
    @test matrix(cm) ≈ [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
    # the two dropped rows add nothing to either counts or weights
    @test sum(matrix(cm) * Diagonal(population(cm))) ≈ 2.0

    # missing and non-finite contact ages both drop, by different branches of assign_group
    missing_survey = ContactSurvey(
        DataFrame(part_id=[1], part_age=[10.0]),
        DataFrame(part_id=[1, 1, 1, 1], cnt_age=[10.0, missing, NaN, Inf]),
    )
    cm_missing = @test_logs (:warn, r"dropped 3 contact\(s\): 3 with missing/unmapped contact group") compute_matrix(
        missing_survey, partition)
    @test matrix(cm_missing) ≈ [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]

    # an age above the last limit is NOT dropped: the top group is open-ended
    top_survey = ContactSurvey(
        DataFrame(part_id=[1], part_age=[10.0]),
        DataFrame(part_id=[1], cnt_age=[999.0]),
    )
    cm_top = @test_logs compute_matrix(top_survey, partition)   # asserts: no log records at all
    @test matrix(cm_top) ≈ [0.0 0.0 0.0; 0.0 0.0 0.0; 1.0 0.0 0.0]

    # The source-stratified functor mirrors both reachable drop reasons. Distinct
    # target/source partitions keep the roles observable: the target partition is
    # applied to contacts (rows), the source partition to participants (columns).
    coarse = AgePartition([0, 18])
    ssm = @test_logs(
        (:warn, r"^compute_source_stratified_matrix dropped 2 contact\(s\): 1 with missing/unmapped contact group; 1 with missing/unmapped participant group"),
        compute_source_stratified_matrix(survey, partition, coarse))
    @test size(matrix(ssm)) == (3, 2)
    @test population(ssm) == [1.0, 1.0]
    @test matrix(ssm) ≈ [1.0 0.0; 0.0 1.0; 0.0 0.0]
end

@testset "unknown participant id is dropped after in-place mutation" begin
    # ContactSurvey's constructor validates every contact's part_id against
    # participants, so this drop reason can't be triggered through the constructor,
    # filter_survey, or subset_survey — the only ways to build a survey. It exists to
    # protect compute_matrix/compute_source_stratified_matrix against a survey mutated
    # after construction: DataFrame fields are public and mutable, so external code
    # (e.g. a future package integration holding a ContactSurvey reference, such as
    # CategoricalInterventions.jl's planned ContACT extension) could append a contact
    # row directly to survey.contacts without going through validation. This test pins
    # that the drop-and-warn behaviour still holds if that ever happens.
    partition = AgePartition([0, 18, 65])
    survey = ContactSurvey(
        DataFrame(part_id=[1], part_age=[10.0]),
        DataFrame(part_id=[1], cnt_age=[10.0]),
    )
    push!(survey.contacts, (part_id=99, cnt_age=10.0))
    cm = @test_logs (:warn, r"^compute_matrix dropped 1 contact\(s\): 1 with unknown participant id") compute_matrix(survey, partition)
    @test matrix(cm) ≈ [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    @test population(cm) == [1.0, 0.0, 0.0]   # the phantom row adds no participant weight

    # the same branch exists in the source-stratified functor
    ssm = @test_logs(
        (:warn, r"^compute_source_stratified_matrix dropped 1 contact\(s\): 1 with unknown participant id"),
        compute_source_stratified_matrix(survey, partition, partition))
    @test matrix(ssm) ≈ [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
end

@testset "Coarsening" begin
    fine = AgePartition([0, 5, 18, 65])
    M = [4.0 1.0 0.5 0.2;
         1.0 3.0 1.0 0.3;
         0.5 1.0 2.5 0.8;
         0.2 0.3 0.8 1.5]
    pop = [500.0, 1000.0, 3000.0, 800.0]
    cm = ContactMatrix(M, fine, pop)

    # Coarsen to 2 groups
    coarse = AgePartition([0, 18])
    cm_coarse = coarsen(cm, coarse)
    @test n_groups(cm_coarse) == 2
    @test size(matrix(cm_coarse)) == (2, 2)

    # Population should be preserved
    @test sum(population(cm_coarse)) ≈ sum(pop)

    # Operator syntax
    cm_coarse2 = cm ↓ coarse
    @test matrix(cm_coarse2) ≈ matrix(cm_coarse)
end

@testset "Composition (⊕)" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M_home = [2.0 0.5 0.3; 0.5 1.0 0.2; 0.3 0.2 1.5]
    M_work = [0.0 0.0 0.0; 0.0 2.0 0.5; 0.0 0.5 0.5]
    cm_home = ContactMatrix(M_home, p, pop)
    cm_work = ContactMatrix(M_work, p, pop)

    total = cm_home ⊕ cm_work
    @test matrix(total) ≈ M_home + M_work

    # Associativity
    M_school = [1.0 0.1 0.0; 0.1 0.5 0.0; 0.0 0.0 0.0]
    cm_school = ContactMatrix(M_school, p, pop)
    @test matrix((cm_home ⊕ cm_work) ⊕ cm_school) ≈ matrix(cm_home ⊕ (cm_work ⊕ cm_school))

    # Commutativity
    @test matrix(cm_home ⊕ cm_work) ≈ matrix(cm_work ⊕ cm_home)

    # Numeric types promote naturally
    cm_int = ContactMatrix([1 2 3; 4 5 6; 7 8 9], p, [1000, 3000, 500])
    mixed = cm_int ⊕ cm_home
    @test eltype(matrix(mixed)) == Float64
    @test matrix(mixed) ≈ [1 2 3; 4 5 6; 7 8 9] .+ M_home

    @testset "mismatched inputs are rejected" begin
        p2 = AgePartition([0, 10, 65])
        cm_diff_partition = ContactMatrix(M_home, p2, pop)
        @test_throws ArgumentError cm_home ⊕ cm_diff_partition

        cm_diff_pop = ContactMatrix(M_home, p, pop .+ 1.0)
        @test_throws ArgumentError cm_home ⊕ cm_diff_pop

        # Matching semantics is the third documented precondition; it is enforced by the
        # shared type parameter on compose_matrices, so the rejection is a MethodError.
        @test_throws MethodError cm_home ⊕ to_counts(cm_home)

        @test_throws ArgumentError compose_matrices(ContactMatrix[])
    end
end

@testset "Stratification (⊗)" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ContactMatrix(M, p, pop)

    # Uniform mixing across 2 strata
    coupling = [0.7 0.3; 0.3 0.7]
    cm_strat = cm ⊗ coupling
    @test n_groups(cm_strat) == 6  # 2 strata × 3 ages
    @test cm_strat.partition isa ProductPartition
    @test group_labels(cm_strat)[1:3] == ["S1:[0,18)", "S1:[18,65)", "S1:65+"]

    # Block structure check: diagonal blocks should have coupling[i,i] * M
    M_strat = matrix(cm_strat)
    @test M_strat[1:3, 1:3] ≈ 0.7 .* M
    @test M_strat[4:6, 1:3] ≈ 0.3 .* M

    stratum_pop = [600.0 400.0;
                   2000.0 1000.0;
                   300.0 200.0]
    cm_strat_pop = stratify(cm, coupling;
        stratum_populations=stratum_pop,
        stratum_labels=["North", "South"],
    )
    @test population(cm_strat_pop) == vec(stratum_pop)
    @test group_labels(cm_strat_pop)[1:3] == ["North:[0,18)", "North:[18,65)", "North:65+"]

    @test_throws ArgumentError cm ⊗ [1.0 -0.1; 0.0 1.0]
    @test_throws ArgumentError stratify([cm, ContactMatrix(M, AgePartition([0, 10, 65]), pop)], coupling)
    @test_throws ArgumentError stratify([cm, ContactMatrix(M, p, pop, PerCapitaRate())], coupling)
end

@testset "Symmetrisation" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.5 3.0 1.0; 0.3 0.8 1.5]
    cm = ContactMatrix(M, p, pop)

    cm_sym = symmetrise(cm)
    M_sym = matrix(cm_sym)

    # Reciprocity: M_sym[i,j] * N_j = M_sym[j,i] * N_i
    for i in 1:3
        for j in 1:3
            @test M_sym[i, j] * pop[j] ≈ M_sym[j, i] * pop[i] atol=1e-10
        end
    end

    # Idempotence: symmetrise(symmetrise(M)) == symmetrise(M)
    cm_sym2 = symmetrise(cm_sym)
    @test matrix(cm_sym2) ≈ M_sym atol=1e-10
    @test matrix(↔(cm)) ≈ M_sym

    zero_pop_ok = ContactMatrix([0.0 0.0; 0.0 4.0], AgePartition([0, 18]), [0.0, 10.0])
    zero_pop_sym = symmetrise(zero_pop_ok)
    @test matrix(zero_pop_sym) == [0.0 0.0; 0.0 4.0]
    @test_throws ArgumentError symmetrise(
        ContactMatrix([1.0 2.0; 3.0 4.0], AgePartition([0, 18]), [0.0, 10.0])
    )
end

@testset "Demographic reprojection" begin
    p = AgePartition([0, 18])
    N = [100.0, 200.0]
    M = [0.0 1.0; 2.0 0.0]              # already reciprocal at N
    cm = ContactMatrix(M, p, N)

    Np = [50.0, 400.0]
    cm_proj = reproject(cm, Np)
    M_proj = matrix(cm_proj)

    # Hand-computed: M'[i,j] = (M[i,j]*N'[j] + M[j,i]*N'[i]) / (2*N'[j])
    @test M_proj ≈ [0.0 0.625; 5.0 0.0]
    @test population(cm_proj) == Np

    # Reciprocity: M'[i,j] * N'[j] == M'[j,i] * N'[i], for any input (reciprocal or not)
    M_nonrecip = [0.0 1.0; 3.0 0.0]     # not reciprocal at N
    cm_nonrecip = ContactMatrix(M_nonrecip, p, N)
    M_proj2 = matrix(reproject(cm_nonrecip, Np))
    for i in 1:2, j in 1:2
        @test M_proj2[i, j] * Np[j] ≈ M_proj2[j, i] * Np[i] atol=1e-10
    end

    # Identity: reprojecting a matrix already reciprocal at N onto N itself is a no-op
    @test matrix(reproject(cm, N)) ≈ M

    # Zero target population: zero contacts pass through, nonzero contacts throw
    zero_target_ok = reproject(ContactMatrix([0.0 0.0; 0.0 4.0], p, N), [0.0, 10.0])
    @test matrix(zero_target_ok) == [0.0 0.0; 0.0 4.0]
    @test_throws ArgumentError reproject(ContactMatrix([1.0 2.0; 3.0 4.0], p, N), [0.0, 10.0])

    # Target population length must match the partition
    @test_throws DimensionMismatch reproject(cm, [1.0, 2.0, 3.0])

    # Not functorial in the population: routing through an intermediate population
    # differs from reprojecting straight to the final one.
    Npp = [1.0, 1.0]
    two_step = matrix(reproject(reproject(cm, Np), Npp))
    one_step = matrix(reproject(cm, Npp))
    @test two_step ≉ one_step

    # Does not commute with coarsening: reproject-then-coarsen differs from
    # coarsen-then-reproject, even onto matching populations.
    fine = AgePartition([0, 10, 20, 30])
    coarse = AgePartition([0, 20])
    N4 = [100.0, 200.0, 150.0, 50.0]
    M4 = [0.0 1.0 0.5 0.2; 2.0 0.0 0.3 0.1; 0.4 0.6 0.0 0.8; 0.1 0.2 0.9 0.0]
    cm4 = ContactMatrix(M4, fine, N4)
    N4p = [80.0, 300.0, 60.0, 90.0]
    N_coarse_target = [N4p[1] + N4p[2], N4p[3] + N4p[4]]

    reproject_then_coarsen = coarsen(reproject(cm4, N4p), coarse)
    coarsen_then_reproject = reproject(coarsen(cm4, coarse), N_coarse_target)
    @test population(reproject_then_coarsen) ≈ population(coarsen_then_reproject)
    @test matrix(reproject_then_coarsen) ≉ matrix(coarsen_then_reproject)
end

@testset "Population transport" begin
    p = AgePartition([0, 18])
    N = [100.0, 200.0]
    M = [0.0 1.0; 2.0 0.0]              # reciprocal at N
    cm = ContactMatrix(M, p, N)
    Np = [50.0, 400.0]

    # Identity: transporting onto its own population is exact, unlike reproject
    @test matrix(transport_population(cm, N)) ≈ M
    @test population(transport_population(cm, N)) == N

    # Composition is exact (diagonal-group action), unlike reproject
    Npp = [10.0, 20.0]
    two_step = matrix(transport_population(transport_population(cm, Np), Npp))
    one_step = matrix(transport_population(cm, Npp))
    @test two_step ≈ one_step

    # Total contacts preserved entrywise: M'[i,j]*N'[j] == M[i,j]*N[j]
    transported = matrix(transport_population(cm, Np))
    for i in 1:2, j in 1:2
        @test transported[i, j] * Np[j] ≈ M[i, j] * N[j]
    end

    # Reciprocity is preserved when the source is already reciprocal at its own
    # population (conditional, unlike reproject's unconditional repair)
    for i in 1:2, j in 1:2
        @test transported[i, j] * Np[j] ≈ transported[j, i] * Np[i] atol=1e-10
    end

    # Non-reciprocal input: transport carries the imbalance through unchanged in
    # total-contacts terms rather than repairing it (that's `reproject`'s job)
    M_nonrecip = [0.0 1.0; 3.0 0.0]      # not reciprocal at N
    cm_nonrecip = ContactMatrix(M_nonrecip, p, N)
    transported_nonrecip = matrix(transport_population(cm_nonrecip, Np))
    @test transported_nonrecip[1, 2] * Np[2] ≉ transported_nonrecip[2, 1] * Np[1]
    @test transported_nonrecip[1, 2] * Np[2] ≈ M_nonrecip[1, 2] * N[2]
    @test transported_nonrecip[2, 1] * Np[1] ≈ M_nonrecip[2, 1] * N[1]

    # Invertible on the reciprocal fibre: round trip recovers the original exactly
    roundtrip = matrix(transport_population(transport_population(cm, Np), N))
    @test roundtrip ≈ M

    # Target population length must match the partition
    @test_throws DimensionMismatch transport_population(cm, [1.0, 2.0, 3.0])

    # Zero (or negative) target population is rejected outright: unlike
    # `reproject`, transport has no zero-population branch to fall back on
    @test_throws ArgumentError transport_population(cm, [0.0, 400.0])
    @test_throws ArgumentError transport_population(cm, [-1.0, 400.0])

    # Identity is bitwise exact, so these assert `==` rather than `≈`. A column
    # whose population is unchanged is copied rather than rescaled, because
    # `0.1*3.0/3.0` is `0.10000000000000002`.
    p_exact = AgePartition([0, 18])
    cm_exact = ContactMatrix([0.1 1.0; 0.7 2.0], p_exact, [3.0, 5.0])
    identity_transport = transport_population(cm_exact, [3.0, 5.0])
    @test matrix(identity_transport) == matrix(cm_exact)
    @test population(identity_transport) == population(cm_exact)

    # The copy is per column, so a partial demographic update leaves the groups
    # whose population did not change bitwise unchanged.
    partial = transport_population(cm_exact, [3.0, 10.0])
    @test matrix(partial)[:, 1] == matrix(cm_exact)[:, 1]
    @test matrix(partial)[:, 2] ≈ [0.5, 1.0]

    # An empty source group is admitted when it carries no contacts: there is
    # nothing to carry forward, and the transported column stays zero.
    empty_ok = ContactMatrix([0.0 1.0; 0.0 2.0], p_exact, [0.0, 10.0])
    transported_empty = transport_population(empty_ok, [5.0, 20.0])
    @test matrix(transported_empty)[:, 1] == [0.0, 0.0]
    @test matrix(transported_empty)[:, 2] ≈ [0.5, 1.0]

    # An empty source group carrying nonzero contacts is rejected. Tagging the
    # same numbers differently must not change whether the call succeeds, since
    # transport accepts exactly the matrices `reinterpret_units` accepts.
    empty_bad = [5.0 7.0; 2.0 1.0]
    @test_throws ArgumentError transport_population(
        ContactMatrix(empty_bad, p_exact, [0.0, 10.0]), [5.0, 10.0])
    @test_throws ArgumentError transport_population(
        ContactMatrix(empty_bad, p_exact, [0.0, 10.0], ContactCounts()), [5.0, 10.0])

    # `reinterpret_units` only checks the axes it rescales, and `PerCapitaRate`
    # rescales rows alone. A matrix whose empty group has a zero row but a
    # nonzero column therefore reaches transport in every representation, and
    # each must reject it on the column.
    zero_row_nonzero_col = [0.0 0.0; 2.0 0.0]
    for semantics in (MeanContacts(), ContactCounts(), PerCapitaRate())
        @test_throws ArgumentError transport_population(
            ContactMatrix(zero_row_nonzero_col, p_exact, [0.0, 10.0], semantics),
            [5.0, 10.0])
    end
end

@testset "transport: coarsening equivariance" begin
    fine = AgePartition([0, 10, 20, 30])
    coarse = AgePartition([0, 18])
    N4 = [100.0, 200.0, 150.0, 50.0]
    N4p = [80.0, 300.0, 60.0, 90.0]
    M4 = [0.0 1.0 0.5 0.2; 2.0 0.0 0.3 0.1; 0.4 0.6 0.0 0.8; 0.1 0.2 0.9 0.0]
    cm4 = ContactMatrix(M4, fine, N4)

    # Coarsening adds up the contacts of the groups it merges, and transport
    # leaves those totals alone, so doing the two in either order gives the same
    # matrix — provided the coarse target population is the fine one added up the
    # same way. This holds however the fine groups are assigned to coarse ones,
    # which is why the loop also covers assignments that interleave them rather
    # than only ones that keep them adjacent. Reprojection has no matching law;
    # the "Demographic reprojection" testset pins that failure.
    for assignment in ([1, 1, 2, 2], [1, 2, 1, 2], [2, 1, 2, 1])
        f = PartitionMap(fine, coarse, assignment)
        coarse_target = zeros(2)
        for j in 1:4
            coarse_target[assignment[j]] += N4p[j]
        end

        transport_then_coarsen = coarsen(transport_population(cm4, N4p), f)
        coarsen_then_transport = transport_population(coarsen(cm4, f), coarse_target)
        @test population(transport_then_coarsen) ≈ population(coarsen_then_transport)
        @test matrix(transport_then_coarsen) ≈ matrix(coarsen_then_transport) atol=1e-10
        @test transport_then_coarsen.semantics isa MeanContacts
    end
end

@testset "transport: symmetrisation commutativity" begin
    p = AgePartition([0, 18])
    N = [100.0, 200.0]
    Np = [50.0, 400.0]
    cm = ContactMatrix([0.0 1.0; 3.0 0.0], p, N)     # not reciprocal at N

    # Symmetrisation replaces total contacts by their symmetric part, transport
    # leaves total contacts alone, so the order does not matter. Stronger than
    # transport's conditional reciprocity preservation, which says nothing about
    # input that is not already reciprocal.
    @test matrix(transport_population(symmetrise(cm), Np)) ≈
          matrix(symmetrise(transport_population(cm, Np))) atol=1e-10
end

@testset "transport and reprojection: additivity over ⊕" begin
    p = AgePartition([0, 18])
    N = [100.0, 200.0]
    Np = [50.0, 400.0]
    home = ContactMatrix([0.0 1.0; 2.0 0.0], p, N)
    work = ContactMatrix([1.0 0.5; 0.25 2.0], p, N)

    # Both are linear in the matrix at a fixed target population, so settings can
    # be composed before or after either one. `⊕` requires its arguments to share
    # a population, which is why both parts are built at `N`; reprojection itself
    # never reads the source population.
    @test matrix(transport_population(home ⊕ work, Np)) ≈
          matrix(transport_population(home, Np) ⊕ transport_population(work, Np)) atol=1e-10
    @test matrix(reproject(home ⊕ work, Np)) ≈
          matrix(reproject(home, Np) ⊕ reproject(work, Np)) atol=1e-10
end

@testset "reprojection: idempotence at a fixed target" begin
    p = AgePartition([0, 18])
    N = [100.0, 200.0]
    Np = [50.0, 400.0]
    cm = ContactMatrix([0.0 1.0; 3.0 0.0], p, N)     # not reciprocal at N

    # Reprojection always returns a matrix reciprocal at the target population,
    # and is the identity on such matrices, so repeating it at the same target
    # changes nothing. This is the case of the composition law where both targets
    # coincide; the general case fails, as the "Demographic reprojection" testset
    # asserts.
    once = reproject(cm, Np)
    twice = reproject(once, Np)
    @test matrix(twice) ≈ matrix(once) atol=1e-10
    @test population(twice) == Np
    @test twice.semantics isa MeanContacts
end

# ---------------------------------------------------------------------------
# POLYMOD UK: cross-validation against socialmixr's reprojection formula
# ---------------------------------------------------------------------------
# Uses a pre-computed R reference matrix committed as a CSV fixture
# (test/fixtures/reprojection/). build_raw_matrix.jl builds the raw matrix
# below from the bundled POLYMOD UK data; compute_reference.R applies
# socialmixr's own normalise_weighted_matrix formula to it and writes
# reprojection_polymod_matrix_R.csv. This validates that ContACT.jl's
# `reproject` produces identical results to socialmixr's own formula.
#
# Data is bundled locally (data/polymod_uk_*.csv); no download is needed.
@testset "reproject matches socialmixr's normalise_weighted_matrix (POLYMOD UK)" begin
    data_dir = joinpath(@__DIR__, "..", "data")
    fixture_dir = joinpath(@__DIR__, "fixtures", "reprojection")

    parse_age(x::AbstractString) = x == "NA" ? missing : parse(Float64, x)
    parse_age(x::Real) = Float64(x)
    parse_age(::Missing) = missing

    participants = CSV.read(joinpath(data_dir, "polymod_uk_participants.csv"), DataFrame)
    contacts = CSV.read(joinpath(data_dir, "polymod_uk_contacts.csv"), DataFrame)
    rename!(participants, :part_age_exact => :part_age)
    participants.part_age = Float64.(participants.part_age)
    cnt_exact = parse_age.(contacts.cnt_age_exact)
    cnt_min = parse_age.(contacts.cnt_age_est_min)
    cnt_max = parse_age.(contacts.cnt_age_est_max)
    contacts.cnt_age = coalesce.(cnt_exact, (cnt_min .+ cnt_max) ./ 2)
    select!(contacts, [:part_id, :cnt_age])
    dropmissing!(contacts, :cnt_age)
    valid_ids = Set(participants.part_id)
    filter!(row -> row.part_id in valid_ids, contacts)

    survey = ContactSurvey(participants, contacts)
    age_partition = AgePartition([0, 18, 65])
    uk_pop_age = [11000.0, 33000.0, 9500.0]
    cm_age = compute_matrix(survey, age_partition; population=uk_pop_age)

    target_pop = [10500.0, 34000.0, 13500.0]
    proj = reproject(cm_age, target_pop)

    ref_df = CSV.read(joinpath(fixture_dir, "reprojection_polymod_matrix_R.csv"), DataFrame)
    ref = Matrix{Float64}(ref_df[:, 2:end])

    @test isapprox(matrix(proj), ref; atol=1e-4)
end

@testset "Refinement" begin
    # Start with a coarse matrix, refine it
    coarse = AgePartition([0, 18, 65])
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    pop_coarse = [2000.0, 4000.0, 1000.0]
    cm = ContactMatrix(M, coarse, pop_coarse)

    fine = AgePartition([0, 5, 18, 45, 65])
    fine_pop = [800.0, 1200.0, 2000.0, 2000.0, 1000.0]
    cm_fine = refine(cm, fine, fine_pop)

    @test n_groups(cm_fine) == 5
    @test population(cm_fine) == fine_pop
end

@testset "Activity refinement" begin
    participants = DataFrame(
        part_id = 1:4,
        part_age = [10.0, 10.0, 30.0, 30.0],
        part_sex = ["F", "F", "M", "M"],
        score = [1.0, 3.0, 1.0, 3.0],
    )
    contacts = DataFrame(
        part_id = [1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4],
        cnt_age = [10.0, 30.0, 10.0, 10.0, 30.0, 30.0,
                   10.0, 30.0, 10.0, 10.0, 30.0, 30.0],
        cnt_sex = ["F", "M", "F", "F", "M", "M",
                   "F", "M", "F", "F", "M", "M"],
    )
    survey = ContactSurvey(participants, contacts)
    base = AgePartition([0, 18])
    cm = survey ▷ base
    @test matrix(cm) ≈ [1.5 1.5; 1.5 1.5]

    spec = ActivityRefinement(survey; n=2, mixing=:assortative, score_col=:score)
    refined = activity_refine(cm, spec)
    @test refined.partition isa ProductPartition
    @test group_labels(refined) == ["[0,18):low", "[0,18):high", "18+:low", "18+:high"]
    @test population(refined) == [1.0, 1.0, 1.0, 1.0]
    @test matrix(refined ↓ base) ≈ matrix(cm)
    @test matrix(cm ↑ spec) ≈ matrix(refined)
    @test matrix(cm ⤊ spec) ≈ matrix(refined)

    C = matrix(refined) * Diagonal(population(refined))
    @test C ≈ transpose(C)

    for mixing in (:disassortative, :proportionate)
        refined_mixing = activity_refine(survey, cm; n=2, mixing=mixing, score_col=:score)
        @test matrix(refined_mixing ↓ base) ≈ matrix(cm)
        C_mixing = matrix(refined_mixing) * Diagonal(population(refined_mixing))
        @test C_mixing ≈ transpose(C_mixing)
    end

    sex = CategoricalPartition(:sex;
        participant_col=:part_sex,
        contact_col=:cnt_sex,
        levels=["F", "M"],
    )
    sex_activity = activity_partition(survey, sex; cutpoints=[2.0], score_col=:score)
    @test sex_activity isa ProductPartition
    @test group_labels(sex_activity) == ["F:low", "F:high", "M:low", "M:high"]
    refined_sex = activity_refine(survey, survey ▷ sex;
        cutpoints=[2.0],
        mixing=:proportionate,
        score_col=:score,
    )
    @test matrix(refined_sex ↓ sex) ≈ matrix(survey ▷ sex)

    row = [1.0, 3.6]
    col = [1.2, 3.4]
    @test activity_mixing_plan(row, col, :assortative) ≈ [1.0 0.0; 0.2 3.4]
    @test activity_mixing_plan(row, col, :disassortative) ≈ [0.0 1.0; 1.2 2.4]
    @test activity_mixing_plan(row, col, :proportionate) ≈ row * transpose(col) ./ sum(row)

    @testset "activity_mixing_plan rejects invalid marginals" begin
        @test_throws DimensionMismatch activity_mixing_plan([1.0, 2.0], [1.0, 1.0, 1.0])
        @test_throws ArgumentError activity_mixing_plan([-1.0, 2.0], [1.0, 0.0])  # row branch
        @test_throws ArgumentError activity_mixing_plan([1.0, 0.0], [-1.0, 2.0])  # column branch
        @test_throws ArgumentError activity_mixing_plan([NaN, 1.0], [1.0, 0.0])   # non-finite
        @test_throws ArgumentError activity_mixing_plan([1.0, 2.0], [1.0, 1.0])   # totals 3 vs 2
    end

    @test_throws ArgumentError activity_refine(
        ContactMatrix([1.0 2.0; 0.0 1.0], base, [2.0, 2.0]),
        spec,
    )
end

@testset "Generalized contact matrices" begin
    age = AgePartition([0, 18])
    pop = [100.0, 200.0]
    base = ContactMatrix([2.0 0.5; 1.0 1.5], age, pop)
    ses = CategoricalPartition(:ses;
        levels=["low", "middle", "high"],
        labels=["low", "middle", "high"],
    )
    dist = [0.35, 0.45, 0.20]

    random_spec = GeneralizedLift(ses; distribution=dist, mixing=:random)
    g_random = base ⊠ random_spec
    @test g_random.partition isa ProductPartition
    @test group_labels(g_random) == [
        "[0,18):low", "[0,18):middle", "[0,18):high",
        "18+:low", "18+:middle", "18+:high",
    ]
    @test population(g_random) ≈ product_population(pop, dist)
    @test matrix(g_random ↓ age) ≈ matrix(base)
    @test matrix(base ↑ random_spec) ≈ matrix(g_random)
    @test matrix(generalize(base, ses; distribution=dist)) ≈ matrix(g_random)
    @test ρ(g_random) ≈ ρ(base)

    activity = [0.2, 0.4, 0.4]
    assortativity = [0.60, 0.50, 0.65]
    split = [0.60, 0.60, 0.50]
    assortative_spec = GeneralizedLift(ses;
        distribution=dist,
        mixing=AssortativeDimensionMixing(activity, assortativity; offdiag_split=split),
    )
    g_assortative = base ⊠ assortative_spec
    @test matrix(g_assortative ↓ age) ≈ matrix(base)
    total_assortative = matrix(g_assortative) * Diagonal(population(g_assortative))
    @test total_assortative ≈ transpose(total_assortative)
    @test !(ρ(g_assortative) ≈ ρ(base))

    explicit_block = BlockMixing(dist * transpose(dist))
    pop_matrix = [35.0 45.0 20.0;
                  70.0 90.0 40.0]
    g_joint = generalize(base, ses, pop_matrix; mixing=explicit_block)
    @test population(g_joint) ≈ [35.0, 45.0, 20.0, 70.0, 90.0, 40.0]
    @test matrix(g_joint ↓ age) ≈ matrix(base)
    @test n_groups(g_joint ↓ ses) == 3

    K = next_generation_matrix(base; transmissibility=0.5, recovery_rate=0.25)
    @test K ≈ 2.0 .* matrix(base)
    @test R₀(base; transmissibility=0.5, recovery_rate=0.25) ≈ 2.0 * ρ(base)
    beta = calibrate_transmissibility(base, 2.7; recovery_rate=0.25)
    @test R0(base; transmissibility=beta, recovery_rate=0.25) ≈ 2.7

    @test_throws ArgumentError GeneralizedLift(ses; distribution=[0.5, 0.5, 0.5])
    @test_throws ArgumentError generalize(base, ses, [10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    @test_throws ArgumentError BlockMixing([0.5 0.5; 0.5 0.5])

    @testset "assortative mixing kernel guards" begin
        ses2 = CategoricalPartition(:ses2; levels=["low", "high"])
        ses3 = CategoricalPartition(:ses3; levels=["low", "middle", "high"])
        ses4 = CategoricalPartition(:ses4; levels=["a", "b", "c", "d"])

        # The constructor itself only checks the simplex and [0,1] ranges ...
        unequal = AssortativeDimensionMixing([0.5, 0.5], [0.9, 0.1])
        @test unequal isa AssortativeDimensionMixing
        # ... the feasibility condition (1-r₁)·a₁ == (1-r₂)·a₂ is enforced during the lift.
        @test_throws ArgumentError generalize(base, ses2; distribution=[0.5, 0.5], mixing=unequal)

        # Balanced cross-group mass is accepted, and the lift still coarsens back.
        balanced = AssortativeDimensionMixing([0.5, 0.5], [0.6, 0.6])
        g2 = generalize(base, ses2; distribution=[0.5, 0.5], mixing=balanced)
        @test matrix(g2 ↓ age) ≈ matrix(base)

        # Four groups: automatic assortative diagonal blocks are not implemented.
        @test_throws ArgumentError generalize(base, ses4;
            distribution=fill(0.25, 4),
            mixing=AssortativeDimensionMixing(fill(0.25, 4), fill(0.5, 4)))

        # offdiag_split is the three-group Manna construction only. Reaching this guard needs
        # a diagonal block that *does* build (n=2, balanced), so the failure is the split itself.
        @test_throws ArgumentError generalize(base, ses2;
            distribution=[0.5, 0.5],
            mixing=AssortativeDimensionMixing([0.5, 0.5], [0.6, 0.6]; offdiag_split=[0.5, 0.5]))

        # activity length must match the added partition (src/generalized_lift.jl:339-340)
        @test_throws DimensionMismatch generalize(base, ses3;
            distribution=[0.3, 0.3, 0.4],
            mixing=AssortativeDimensionMixing([0.5, 0.5], [0.6, 0.6]))

        # n=3 parameters implying negative off-diagonal mass (src/generalized_lift.jl:407-410)
        @test_throws ArgumentError generalize(base, ses3;
            distribution=[0.3, 0.3, 0.4],
            mixing=AssortativeDimensionMixing([0.1, 0.1, 0.8], [0.0, 0.0, 0.0]))
    end
end

@testset "Utilities" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ContactMatrix(M, p, pop)

    # Per capita conversion divides by the *contactee* (row) population, matching
    # socialmixr's `per_capita` once the row/column conventions are aligned.
    cm_pc = to_per_capita(cm)
    @test cm_pc.semantics isa PerCapitaRate
    @test matrix(cm_pc) ≈ [M[i, j] / pop[i] for i in 1:3, j in 1:3]

    # Counts conversion
    cm_counts = to_counts(cm)
    @test cm_counts.semantics isa ContactCounts
    @test matrix(cm_counts) ≈ [M[i, j] * pop[j] for i in 1:3, j in 1:3]

    # Conversions dispatch on the semantics of their input, so counts → means is
    # `to_mean_contacts`, and counts → per-capita divides by both populations.
    @test matrix(to_mean_contacts(cm_counts)) ≈ M
    @test matrix(to_per_capita(cm_counts)) ≈ [M[i, j] / pop[i] for i in 1:3, j in 1:3]

    # Zero-population guards are per-direction: `to_counts` rescales columns,
    # `to_per_capita` rescales rows.
    @test_throws ArgumentError to_counts(
        ContactMatrix([0.0 1.0; 0.0 0.0], AgePartition([0, 18]), [10.0, 0.0])
    )
    @test_throws ArgumentError to_per_capita(
        ContactMatrix([0.0 0.0; 1.0 0.0], AgePartition([0, 18]), [10.0, 0.0])
    )

    # Spectral radius
    sr = spectral_radius(cm)
    @test sr > 0
    @test sr ≈ maximum(abs.(eigvals(M)))
end

# ─────────────────────────────────────────────────────────────────────────────
# Unit semantics as a groupoid of representations, and naturality of the
# operations defined on it. Formalised in proofs/ContACTProofs/UnitSemantics.lean
# (`rescale_one`/`rescale_rescale`/`rescale_inv`, `symmetrise_naturality`,
# `symmetriseAt_counts`, `symmetriseAt_perCapita`,
# `symmetrise_not_natural_naive`).
# ─────────────────────────────────────────────────────────────────────────────
@testset "Unit semantics" begin
    p = AgePartition([0, 18, 65])
    # Deliberately unequal populations: every failure mode here is invisible
    # when all N_j agree.
    N = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ContactMatrix(M, p, N)

    reps = (MeanContacts(), ContactCounts(), PerCapitaRate())

    # Exponent the reciprocity law actually sees: e = b - a.
    recip_exp(s) = population_exponents(s)[2] - population_exponents(s)[1]

    @testset "representations form a groupoid" begin
        # Exponents (a,b) on total contacts: R[i,j] = T[i,j]*N_i^-a*N_j^-b.
        # Admissible representations are corners of the unit square, since a
        # contact involves one participant and one contactee.
        @test population_exponents(ContactCounts()) == (0, 0)
        @test population_exponents(MeanContacts()) == (0, 1)
        @test population_exponents(PerCapitaRate()) == (1, 1)
        for s in reps
            a, b = population_exponents(s)
            @test a in (0, 1) && b in (0, 1)
        end

        for s in reps
            cs = reinterpret_units(cm, s)
            @test cs.semantics isa typeof(s)
            # Identity: converting to the representation already held.
            @test matrix(reinterpret_units(cs, s)) ≈ matrix(cs)
            # Inverse.
            @test matrix(reinterpret_units(cs, MeanContacts())) ≈ M
            # Composition: any route between two representations agrees.
            for t in reps, u in reps
                direct = reinterpret_units(cs, u)
                viaT = reinterpret_units(reinterpret_units(cs, t), u)
                @test matrix(direct) ≈ matrix(viaT)
            end
        end
    end

    @testset "reciprocity is representation-dependent" begin
        # The general form is R[i,j]*N_j^e == R[j,i]*N_i^e, with e = b - a: the
        # reciprocity law sees only the difference of the two exponents.
        cm_sym = symmetrise(cm)
        for s in reps
            R = matrix(reinterpret_units(cm_sym, s))
            e = recip_exp(s)
            for i in 1:3, j in 1:3
                @test R[i, j] * N[j]^e ≈ R[j, i] * N[i]^e
            end
        end

        C = matrix(to_counts(cm_sym))
        P = matrix(to_per_capita(cm_sym))
        S = matrix(cm_sym)
        for i in 1:3, j in 1:3
            @test S[i, j] * N[j] ≈ S[j, i] * N[i]   # MeanContacts: e = 1
            @test C[i, j] ≈ C[j, i]                 # ContactCounts: e = 0
            @test P[i, j] ≈ P[j, i]                 # PerCapitaRate: e = 0
        end
        # The property PerCapitaRate is named for: a reciprocal matrix is
        # symmetric in this representation. This is Lomas et al.'s β_ij.
        @test P ≈ transpose(P)
    end

    @testset "symmetrise: naturality square" begin
        # reinterpret_units ∘ symmetrise == symmetrise ∘ reinterpret_units,
        # for every pair of representations.
        for s in reps, t in reps
            cs = reinterpret_units(cm, s)
            convert_then_sym = symmetrise(reinterpret_units(cs, t))
            sym_then_convert = reinterpret_units(symmetrise(cs), t)
            @test matrix(convert_then_sym) ≈ matrix(sym_then_convert)
            @test typeof(convert_then_sym.semantics) == typeof(sym_then_convert.semantics)
        end
    end

    @testset "symmetrise: closed form per representation" begin
        # e = 0: plain averaging, no population in the formula.
        C = matrix(to_counts(cm))
        @test matrix(symmetrise(to_counts(cm))) ≈ (C .+ transpose(C)) ./ 2

        P = matrix(to_per_capita(cm))
        @test matrix(symmetrise(to_per_capita(cm))) ≈ (P .+ transpose(P)) ./ 2

        # e = 1: the familiar MeanContacts formula.
        @test matrix(symmetrise(cm)) ≈
            [(M[i, j] * N[j] + M[j, i] * N[i]) / (2 * N[j]) for i in 1:3, j in 1:3]

        # Two representations can share a formula (both e = 0) and still be
        # different objects.
        @test recip_exp(ContactCounts()) == recip_exp(PerCapitaRate()) == 0
        @test !(C ≈ P)
    end

    @testset "PerCapitaRate agrees with the activity lift's β" begin
        # Lomas et al. (2025): β_ij = a_i·a_j/D. ContACT's reciprocal lift is
        # M[i,j] = a_i·N_i·a_j/D, so β is exactly the PerCapitaRate view of it.
        act = [0.8, 1.3, 0.5]
        D = sum(act .* N)
        M_lift = [act[i] * N[i] * act[j] / D for i in 1:3, j in 1:3]
        cm_lift = ContactMatrix(M_lift, p, N)
        β = [act[i] * act[j] / D for i in 1:3, j in 1:3]

        @test matrix(to_per_capita(cm_lift)) ≈ β
        @test β ≈ transpose(β)
        # The lift is already reciprocal, so symmetrisation is a no-op there.
        @test matrix(symmetrise(to_per_capita(cm_lift))) ≈ β
    end

    @testset "symmetrise: idempotent in every representation" begin
        for s in reps
            cs = reinterpret_units(cm, s)
            @test matrix(symmetrise(symmetrise(cs))) ≈ matrix(symmetrise(cs))
        end
    end

    @testset "the naive non-dispatching formula breaks naturality" begin
        # Regression guard: applying the MeanContacts formula to a ContactCounts
        # matrix is what `symmetrise` used to do. It is a different matrix, so a
        # future refactor that drops the dispatch will fail here rather than
        # silently return a mislabelled result.
        C = matrix(to_counts(cm))
        naive = [(C[i, j] * N[j] + C[j, i] * N[i]) / (2 * N[j]) for i in 1:3, j in 1:3]
        correct = matrix(symmetrise(to_counts(cm)))
        @test !(naive ≈ correct)
        # And the naive result is not even reciprocal in the counts sense.
        @test !(naive ≈ transpose(naive))
    end

    @testset "zero populations" begin
        p2 = AgePartition([0, 18])

        # Guards are per-direction, because the conversions rescale different
        # axes. `to_counts` scales columns; `to_per_capita` scales rows.
        bad_col = ContactMatrix([0.0 1.0; 0.0 0.0], p2, [10.0, 0.0])
        @test_throws ArgumentError to_counts(bad_col)
        @test matrix(to_per_capita(bad_col)) ≈ [0.0 0.1; 0.0 0.0]  # rows are fine

        bad_row = ContactMatrix([0.0 0.0; 1.0 0.0], p2, [10.0, 0.0])
        @test_throws ArgumentError to_per_capita(bad_row)
        @test matrix(to_counts(bad_row)) ≈ [0.0 0.0; 10.0 0.0]     # columns are fine

        # `symmetrise` is *more* permissive than conversion, and deliberately
        # so: its guard is on total contacts `M[i,j]·N_j`, which vanish for
        # `bad_col`, whereas conversion needs the raw entry to vanish to stay
        # invertible. So the naturality square above is stated on the
        # conversions' domain, which is the smaller of the two.
        @test matrix(symmetrise(bad_col)) ≈ zeros(2, 2)
        # A nonzero *total* against an empty group is inconsistent for both.
        @test_throws ArgumentError symmetrise(bad_row)

        # An empty group with no contacts is fine in every representation.
        ok = ContactMatrix([1.0 0.0; 0.0 0.0], p2, [10.0, 0.0])
        for s in reps
            os = reinterpret_units(ok, s)
            @test matrix(reinterpret_units(os, MeanContacts())) ≈ matrix(ok)
            @test matrix(symmetrise(os)) ≈ matrix(symmetrise(os))
        end
        # ContactCounts symmetrisation needs no division, so it goes through.
        @test matrix(symmetrise(to_counts(ok))) ≈ [10.0 0.0; 0.0 0.0]
    end
end

# Every structural morphism commutes with the unit-semantics isomorphisms:
#     reinterpret_units(op(cm), s) == op(reinterpret_units(cm, s))
# Any new morphism belongs in this testset.
@testset "Naturality of unit semantics" begin
    fine = AgePartition([0, 18, 65])
    coarse = AgePartition([0, 18])
    N = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ↔(ContactMatrix(M, fine, N))
    cmc = coarsen(cm, coarse)
    cm2 = ↔(ContactMatrix(M .* 0.5, fine, N))
    coupling = [0.7 0.3; 0.3 0.7]

    reps = (MeanContacts(), ContactCounts(), PerCapitaRate())

    Nreproj = [800.0, 1200.0, 900.0]

    morphisms = (
        ("symmetrise (↔)", c -> symmetrise(c),                cm),
        ("coarsen (↓)",    c -> coarsen(c, coarse),           cm),
        ("refine (↑)",     c -> refine(c, fine, N),           cmc),
        ("stratify (⊗)",   c -> stratify(c, coupling),        cm),
        ("reproject",      c -> reproject(c, Nreproj),        cm),
        ("transport_population", c -> transport_population(c, Nreproj), cm),
    )

    for (name, op, obj) in morphisms
        @testset "$name commutes with reinterpret_units" begin
            for s in reps
                lhs = op(reinterpret_units(obj, s))
                rhs = reinterpret_units(op(obj), s)
                @test matrix(lhs) ≈ matrix(rhs)
                @test typeof(lhs.semantics) == typeof(rhs.semantics)
                @test population(lhs) ≈ population(rhs)
                @test group_labels(lhs) == group_labels(rhs)
            end
        end
    end

    @testset "compose (⊕) commutes with reinterpret_units" begin
        for s in reps
            lhs = compose_matrices(reinterpret_units(cm, s), reinterpret_units(cm2, s))
            rhs = reinterpret_units(compose_matrices(cm, cm2), s)
            @test matrix(lhs) ≈ matrix(rhs)
            @test typeof(lhs.semantics) == typeof(rhs.semantics)
        end
    end

    @testset "heterogeneous stratify (⊗) commutes" begin
        for s in reps
            cs = [reinterpret_units(cm, s), reinterpret_units(cm2, s)]
            lhs = stratify(cs, coupling)
            rhs = reinterpret_units(stratify([cm, cm2], coupling), s)
            @test matrix(lhs) ≈ matrix(rhs)
            @test typeof(lhs.semantics) == typeof(rhs.semantics)
        end
    end

    @testset "transport is the identity in the home representation" begin
        # ↔ and ↓ are both written in MeanContacts, and cm is MeanContacts, so _via
        # routes straight to the home implementation. These pin that routing: a wrong
        # home tag, or a formula belonging to another representation, changes the
        # numbers outright. == rather than ≈ holds the arithmetic itself fixed —
        # rewriting a morphism in another representation reassociates it.
        # (_via's short-circuit is not what this pins: reinterpret_units returns early
        # on an identity conversion, so the long path is the same numbers on the same
        # domain, two intermediate wrappers more.)
        @test matrix(symmetrise(cm)) == matrix(ContACT._symmetrise_mean(cm))
        @test matrix(coarsen(cm, coarse)) ==
              matrix(ContACT._coarsen_mean(cm, PartitionMap(fine, coarse)))
        @test matrix(reproject(cm, Nreproj)) ==
              matrix(ContACT._reproject_mean(cm, Nreproj))
    end

    @testset "morphisms compose across representations" begin
        pipeline(c) = symmetrise(refine(coarsen(c, coarse), fine, N))
        for s in reps
            @test matrix(pipeline(reinterpret_units(cm, s))) ≈
                  matrix(reinterpret_units(pipeline(cm), s))
        end
    end
end

# Scalars leave the category, so transport does not apply to them. Each must
# either agree across representations or reject the ones it is not defined on.
@testset "Representation-dependent scalars" begin
    p = AgePartition([0, 18, 65])
    N = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ↔(ContactMatrix(M, p, N))

    @testset "assortativity_index requires MeanContacts" begin
        @test assortativity_index(cm) > 0
        @test_throws ArgumentError assortativity_index(to_counts(cm))
        @test_throws ArgumentError assortativity_index(to_per_capita(cm))
        # Guarded because it is not invariant: the counts reading differs.
        Mc = matrix(to_counts(cm))
        @test sum(Mc[i, i] / sum(Mc[i, :]) for i in 1:3) ≉ assortativity_index(cm)
    end

    @testset "next-generation quantities require MeanContacts" begin
        for f in (next_generation_matrix, basic_reproduction_number, R0, R₀)
            @test_throws ArgumentError f(to_counts(cm))
            @test_throws ArgumentError f(to_per_capita(cm))
        end
        @test_throws ArgumentError r0_bounds(to_counts(cm))
        @test_throws ArgumentError calibrate_transmissibility(to_per_capita(cm), 2.0)
    end

    @testset "spectral_radius is raw and representation-dependent" begin
        @test spectral_radius(cm) ≈ R0(cm)
        @test spectral_radius(to_counts(cm)) ≉ spectral_radius(cm)
        @test spectral_radius(to_per_capita(cm)) ≉ spectral_radius(cm)
    end
end

@testset "ACSet schemas" begin
    @testset "ContactSurveyACSet" begin
        survey = make_test_survey()
        partition = AgePartition([0, 18, 65])
        acs = ContactSurveyACSet(survey, partition)

        @test nparts(acs, :G) == 3
        @test nparts(acs, :P) == 10
        @test nparts(acs, :C) == 20

        # Every contact links to a participant
        for k in 1:nparts(acs, :C)
            reporter_idx = subpart(acs, k, :reporter)
            @test 1 <= reporter_idx <= nparts(acs, :P)
        end

        invalid_age = ContactSurvey(
            DataFrame(part_id=[1], part_age=[-1.0]),
            DataFrame(part_id=[1], cnt_age=[10.0]),
        )
        @test_throws ArgumentError ContactSurveyACSet(invalid_age, partition)

        # Survey with contacts referencing unknown participants is rejected at construction
        @test_throws ArgumentError ContactSurvey(
            DataFrame(part_id=[1], part_age=[10.0]),
            DataFrame(part_id=[2], cnt_age=[10.0]),
        )
    end

    @testset "LabelledContactMatrix ACSet" begin
        p = AgePartition([0, 18, 65])
        pop = [1000.0, 3000.0, 500.0]
        M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
        cm = ContactMatrix(M, p, pop)

        acs = LabelledContactMatrix(cm)
        @test nparts(acs, :G) == 3
        @test nparts(acs, :E) == 9  # 3×3 entries
        @test subpart(acs, 1, :gname) == "[0,18)"
    end

    @testset "migrate_coarsen (functorial data migration)" begin
        survey = make_test_survey()
        fine = AgePartition([0, 18, 45, 65])
        acs_fine = ContactSurveyACSet(survey, fine)

        coarse = AgePartition([0, 45])
        f = AgeMap(fine, coarse)
        acs_coarse = migrate_coarsen(acs_fine, f)

        @test nparts(acs_coarse, :G) == 2
        @test nparts(acs_coarse, :P) == nparts(acs_fine, :P)
        @test nparts(acs_coarse, :C) == nparts(acs_fine, :C)
    end

    @testset "ContactSharer & UWD composition" begin
        M_home = [2.0 0.5; 0.5 1.5]
        M_work = [0.0 0.0; 0.0 2.0]

        diagram = @relation (age,) begin
            home(age)
            work(age)
        end

        sharers = Dict(
            :home => ContactSharer(M_home),
            :work => ContactSharer(M_work)
        )

        result = compose_uwd(diagram, sharers)
        @test result ≈ M_home + M_work
    end
end

@testset "Extended operators" begin
    @testset "ρ (spectral radius)" begin
        p = AgePartition([0, 18, 65])
        M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
        pop = [1000.0, 3000.0, 500.0]
        cm = ContactMatrix(M, p, pop)
        @test ρ(cm) == spectral_radius(cm)
        @test ρ(cm) > 0
    end

    @testset "↓ with AgeMap" begin
        fine = AgePartition([0, 5, 18, 65])
        M = [4.0 1.0 0.5 0.2; 1.0 3.0 1.0 0.3; 0.5 1.0 2.5 0.8; 0.2 0.3 0.8 1.5]
        pop = [500.0, 1000.0, 3000.0, 800.0]
        cm = ContactMatrix(M, fine, pop)
        coarse = AgePartition([0, 18])
        f = AgeMap(fine, coarse)
        # Operator with AgeMap should equal function call
        @test matrix(cm ↓ f) ≈ matrix(coarsen(cm, f))
        @test matrix(cm ↓ f) ≈ matrix(cm ↓ coarse)
    end

    @testset "↑ with RefinementPrior" begin
        coarse = AgePartition([0, 18, 65])
        M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
        pop_coarse = [2000.0, 4000.0, 1000.0]
        cm = ContactMatrix(M, coarse, pop_coarse)

        fine = AgePartition([0, 5, 18, 45, 65])
        fine_pop = [800.0, 1200.0, 2000.0, 2000.0, 1000.0]
        prior = RefinementPrior(fine, fine_pop)

        cm_fine = cm ↑ prior
        @test n_groups(cm_fine) == 5
        @test population(cm_fine) == fine_pop
        @test matrix(cm_fine) ≈ matrix(refine(cm, fine, fine_pop))

        # Dimension mismatch in RefinementPrior
        @test_throws DimensionMismatch RefinementPrior(fine, [1.0, 2.0])
        @test_throws ArgumentError RefinementPrior(fine, [800.0, -1.0, 2000.0, 2000.0, 1000.0])
    end

    @testset "⤊ with ActivityRefinement" begin
        participants = DataFrame(part_id=1:2, part_age=[10.0, 10.0], score=[1.0, 3.0])
        contacts = DataFrame(part_id=[1, 2, 2], cnt_age=[10.0, 10.0, 10.0])
        survey = ContactSurvey(participants, contacts)
        cm = ↔(survey ▷ AgePartition([0]))
        spec = ActivityRefinement(survey; n=2, score_col=:score)
        @test matrix(cm ⤊ spec) ≈ matrix(activity_refine(cm, spec))
    end

    @testset "▷ (functor application)" begin
        survey = make_test_survey()
        partition = AgePartition([0, 18, 65])
        # Operator should equal function call
        cm_op = survey ▷ partition
        cm_fn = compute_matrix(survey, partition)
        @test matrix(cm_op) ≈ matrix(cm_fn)
        @test n_groups(cm_op) == 3
    end

    @testset "∘ (AgeMap composition)" begin
        fine = AgePartition([0, 5, 18, 45, 65])
        medium = AgePartition([0, 18, 65])
        coarse = AgePartition([0, 65])

        f = AgeMap(fine, medium)   # fine → medium
        g = AgeMap(medium, coarse) # medium → coarse
        h = g ∘ f                  # fine → coarse

        @test h.domain.limits == fine.limits
        @test h.codomain.limits == coarse.limits

        # Functoriality: coarsen(cm, g ∘ f) == coarsen(coarsen(cm, f), g)
        M = [4.0 1.0 0.5 0.2 0.1;
             1.0 3.0 1.0 0.3 0.1;
             0.5 1.0 2.5 0.8 0.2;
             0.2 0.3 0.8 1.5 0.4;
             0.1 0.1 0.2 0.4 1.0]
        pop = [500.0, 1000.0, 2000.0, 1500.0, 800.0]
        cm = ContactMatrix(M, fine, pop)

        lhs = cm ↓ h
        rhs = (cm ↓ f) ↓ g
        @test matrix(lhs) ≈ matrix(rhs)

        # Incompatible composition should throw
        wrong = AgeMap(AgePartition([0, 30, 65]), coarse)
        @test_throws ArgumentError wrong ∘ f
    end
end

@testset "Integrated algebraic pipeline" begin
    survey = make_test_survey()
    fine = AgePartition([0, 18, 45, 65])
    medium = AgePartition([0, 18, 65])
    coarse = AgePartition([0, 65])
    coupling = [0.8 0.2; 0.1 0.9]

    f = AgeMap(fine, medium)
    g = AgeMap(medium, coarse)
    h = g ∘ f

    cm_fine = survey ▷ fine
    cm_medium = cm_fine ↓ f
    cm_coarse = cm_fine ↓ h
    cm_two_step = cm_medium ↓ g
    @test matrix(cm_coarse) ≈ matrix(cm_two_step)

    reciprocal = ↔(cm_medium)
    N = population(reciprocal)
    R = matrix(reciprocal)
    @test all(R[i, j] * N[j] ≈ R[j, i] * N[i] for i in 1:n_groups(reciprocal), j in 1:n_groups(reciprocal))

    regional = reciprocal ⊗ coupling
    @test n_groups(regional) == 2 * n_groups(reciprocal)
    @test ρ(regional) > 0

    intervention = reciprocal ⊕ ContactMatrix(zeros(n_groups(reciprocal), n_groups(reciprocal)),
                                               reciprocal.partition,
                                               population(reciprocal))
    @test matrix(intervention) ≈ matrix(reciprocal)
end

@testset "SEP metrics" begin
    # 2×2 product partition: age × SEP
    age = IntervalPartition{:age,Float64}([0.0, 50.0])
    sep = CategoricalPartition(:sep, ["low", "high"])
    prod_part = age × sep
    # 4 groups: (young,low), (young,high), (old,low), (old,high)
    pop = [300.0, 200.0, 250.0, 250.0]
    # Build a matrix with known structure:
    # - diagonal-heavy in SEP dimension (assortative by SEP)
    M = [4.0 1.0 1.0 0.5;   # (young,low) contacted by ...
         1.0 3.0 0.5 1.0;   # (young,high) contacted by ...
         1.0 0.5 4.0 1.0;   # (old,low) contacted by ...
         0.5 1.0 1.0 3.0]   # (old,high) contacted by ...
    cm = ContactMatrix(M, prod_part, pop)

    @testset "marginal_matrix" begin
        # Marginalize to SEP dimension
        cm_sep = marginal_matrix(cm, :sep)
        @test n_groups(cm_sep) == 2
        @test group_labels(cm_sep) == ["low", "high"]
        # Marginal should be 2×2 with age summed out
        @test all(matrix(cm_sep) .> 0)

        # Marginalize to age dimension
        cm_age = marginal_matrix(cm, :age)
        @test n_groups(cm_age) == 2
        @test group_labels(cm_age) == ["[0,50)", "50+"]
    end

    @testset "assortativity_index" begin
        # Purely proportionate mixing 2×2: assortativity = 1.0
        prop_M = [0.6 0.4; 0.6 0.4]  # row-normalized constant → diag/row sums to 1.0
        prop_cm = ContactMatrix(prop_M, sep, [500.0, 500.0])
        ai_prop = assortativity_index(prop_cm)
        @test ai_prop ≈ 1.0

        # Perfect assortative 2×2: assortativity = 2.0
        assort_M = [1.0 0.0; 0.0 1.0]
        assort_cm = ContactMatrix(assort_M, sep, [500.0, 500.0])
        ai_assort = assortativity_index(assort_cm)
        @test ai_assort ≈ 2.0

        # Anti-assortative 2×2: assortativity = 0.0
        anti_M = [0.0 1.0; 1.0 0.0]
        anti_cm = ContactMatrix(anti_M, sep, [500.0, 500.0])
        ai_anti = assortativity_index(anti_cm)
        @test ai_anti ≈ 0.0

        # Product matrix SEP-dimension assortativity
        ai_sep = assortativity_index(cm, :sep)
        @test ai_sep > 1.0  # our matrix is SEP-assortative
    end

    @testset "type_reproduction_number" begin
        # Simple 2×2 next-gen matrix with known Tg
        # K = transmissibility/recovery * M * diag(pop)/sum(pop)
        simple_M = [2.0 1.0; 1.0 2.0]
        simple_pop = [500.0, 500.0]
        simple_cm = ContactMatrix(simple_M, sep, simple_pop)

        # Target group 1 only
        Tg = type_reproduction_number(simple_cm, [1])
        @test Tg > 0
        @test isfinite(Tg)

        # Target all groups → should equal R0
        Tg_all = type_reproduction_number(simple_cm, [1, 2])
        R0_val = basic_reproduction_number(simple_cm)
        @test Tg_all ≈ R0_val

        # Boolean mask target
        Tg_bool = type_reproduction_number(simple_cm, [true, false])
        @test Tg_bool ≈ Tg

        # Label-based target
        Tg_label = type_reproduction_number(simple_cm, ["low"])
        @test Tg_label ≈ Tg
    end

    @testset "control_threshold" begin
        @test control_threshold(0.0) == 0.0
        @test control_threshold(1.0) == 0.0
        @test control_threshold(2.0) ≈ 0.5
        @test control_threshold(4.0) ≈ 0.75
        @test_throws ArgumentError control_threshold(-1.0)
        @test_throws ArgumentError control_threshold(Inf)
    end

    @testset "control_effort" begin
        simple_M = [2.0 1.0; 1.0 2.0]
        simple_pop = [600.0, 400.0]
        simple_cm = ContactMatrix(simple_M, sep, simple_pop)

        # Manual computation
        Tg = type_reproduction_number(simple_cm, [1])
        thresh = control_threshold(Tg)
        effort = control_effort(simple_cm, [1], thresh)
        @test effort ≈ thresh * 600.0 / 1000.0
        @test effort >= 0

        # Keyword form should match
        effort_kw = control_effort(simple_cm, [1])
        @test effort_kw ≈ effort
    end

    @testset "invalid target groups are rejected" begin
        simple_cm = ContactMatrix([2.0 1.0; 1.0 2.0], sep, [500.0, 500.0])

        @test_throws ArgumentError type_reproduction_number(simple_cm, Int[])            # empty collection
        @test_throws ArgumentError type_reproduction_number(simple_cm, [false, false])   # mask selects nothing
        @test_throws ArgumentError type_reproduction_number(simple_cm, [1, 1])           # duplicate index
        @test_throws ArgumentError type_reproduction_number(simple_cm, [3])              # index above n_groups
        @test_throws ArgumentError type_reproduction_number(simple_cm, [0])              # index below 1
        @test_throws ArgumentError type_reproduction_number(simple_cm, ["nonexistent"])  # label matches no group
        @test_throws ArgumentError type_reproduction_number(simple_cm, (1, 2))           # not index/label/mask
        @test_throws DimensionMismatch type_reproduction_number(simple_cm, [true, false, true])

        # control_effort routes through the same helper
        @test_throws ArgumentError control_effort(simple_cm, [3], 0.5)
        @test_throws DimensionMismatch control_effort(simple_cm, [true, false, true], 0.5)
    end

    @testset "PartitionMap product→product projection" begin
        # age × sep → sep projection
        proj = PartitionMap(prod_part, sep)
        cm_sep = coarsen(cm, proj)
        @test n_groups(cm_sep) == 2
        @test group_labels(cm_sep) == ["low", "high"]

        # Verify total contacts preserved
        C_full = matrix(cm) * Diagonal(pop)
        C_sep = matrix(cm_sep) * Diagonal(population(cm_sep))
        @test sum(C_full) ≈ sum(C_sep)

        # age × sep → age projection
        proj_age = PartitionMap(prod_part, age)
        cm_age = coarsen(cm, proj_age)
        @test n_groups(cm_age) == 2
        C_age = matrix(cm_age) * Diagonal(population(cm_age))
        @test sum(C_age) ≈ sum(C_full)
    end
end

@testset "Epidemic bounds" begin
    @testset "Scalar solvers" begin
        # Subcritical: t_α = 0 for α ≤ 1
        @test solve_final_size_scalar(0.5) == 0.0
        @test solve_final_size_scalar(1.0) == 0.0

        # Known values: α=2 → τ ≈ 0.7968
        τ2 = solve_final_size_scalar(2.0)
        @test abs(1 - τ2 - exp(-2.0 * τ2)) < 1e-12
        @test 0.79 < τ2 < 0.80

        # α=3 → τ ≈ 0.9401
        τ3 = solve_final_size_scalar(3.0)
        @test abs(1 - τ3 - exp(-3.0 * τ3)) < 1e-12
        @test 0.93 < τ3 < 0.95

        # Extended: γ=0 reduces to basic
        @test solve_final_size_ext(2.0, 0.0) ≈ τ2 atol=1e-10

        # Extended: γ > 0 always has positive solution
        t_ext = solve_final_size_ext(0.5, 1.0)
        @test t_ext > 0.0
        @test abs(t_ext - (1 - exp(-0.5 * t_ext - 1.0))) < 1e-12

        # Extended: larger γ → larger final size
        t1 = solve_final_size_ext(2.0, 0.5)
        t2 = solve_final_size_ext(2.0, 1.0)
        @test t2 > t1 > τ2
    end

    @testset "Vector final-size solver" begin
        # Homogeneous: K = [m], π = [1] → τ = t_m
        K1 = [2.5;;]
        π1 = [1.0]
        τ1 = solve_final_size_vector(K1, π1)
        @test length(τ1) == 1
        @test τ1[1] ≈ solve_final_size_scalar(2.5) atol=1e-10

        # Subcritical: R₀ ≤ 1
        K_sub = [0.4 0.2; 0.3 0.5]
        π_sub = [0.3, 0.7]
        @test all(solve_final_size_vector(K_sub, π_sub) .== 0.0)

        # 2-type symmetric: K = [a b; b a], π = [0.5, 0.5]
        K2 = [1.5 0.5; 0.5 1.5]
        π2 = [0.5, 0.5]
        τ2 = solve_final_size_vector(K2, π2)
        @test τ2[1] ≈ τ2[2] atol=1e-10  # symmetric → same final sizes
        @test all(τ2 .> 0)  # R₀ = 2 > 1
    end

    @testset "R₀ bounds - general" begin
        # Homogeneous: all row sums equal → bounds collapse
        K_homo = [1.0 1.0; 1.0 1.0]
        b = r0_bounds(K_homo; info=:row)
        @test b.lower ≈ 2.0
        @test b.upper ≈ 2.0

        # Asymmetric
        K_asym = [1.5 0.5; 0.3 1.7]
        b_row = r0_bounds(K_asym; info=:row)
        b_col = r0_bounds(K_asym; info=:col)
        R0_actual = maximum(abs.(eigvals(K_asym)))
        @test b_row.lower ≤ R0_actual ≤ b_row.upper
        @test b_col.lower ≤ R0_actual ≤ b_col.upper

        # Row sums: [2.0, 2.0] → bounds [2,2]
        @test b_row.lower ≈ 2.0
        @test b_row.upper ≈ 2.0

        # Both info
        b_both = r0_bounds(K_asym; info=:both)
        @test b_both.lower ≥ b_row.lower
        @test b_both.upper ≤ b_row.upper
    end

    @testset "R₀ bounds - detailed balance" begin
        # Symmetric K with equal π → detailed balance holds
        K_db = [2.0 0.5; 0.5 1.5]
        π_db = [0.5, 0.5]  # πᵢKᵢⱼ = πⱼKⱼᵢ ✓
        b = r0_bounds_detailed_balance(K_db, π_db; info=:row)
        R0_actual = maximum(abs.(eigvals(K_db)))

        # Lower bound ≥ general lower bound
        b_gen = r0_bounds(K_db; info=:row)
        @test b.lower ≥ b_gen.lower - 1e-10  # DB lower is tighter

        # Bounds contain actual R₀
        @test b.lower ≤ R0_actual + 1e-10
        @test b.upper ≥ R0_actual - 1e-10

        # ContactMatrix wrapper
        part = AgePartition([0, 50])
        cm = ContactMatrix([7.0 3.0; 3.0 5.0], part, [5000.0, 5000.0])
        b_cm = r0_bounds_detailed_balance(cm; info=:row)
        @test b_cm.lower > 0
        @test b_cm.upper ≥ b_cm.lower
    end

    @testset "Final size bounds - column sums" begin
        # 2-type example
        K = [1.8 0.4; 0.6 1.2]
        π = [0.3, 0.7]
        bounds = final_size_bounds(K, π; info=:col)

        # Actual final sizes for comparison
        τ_actual = solve_final_size_vector(K, π)

        # Bounds must contain actual values
        for i in 1:2
            @test bounds.lower[i] ≤ τ_actual[i] + 1e-8
            @test bounds.upper[i] ≥ τ_actual[i] - 1e-8
        end

        # Lower ≤ upper
        @test all(bounds.lower .≤ bounds.upper .+ 1e-10)
    end

    @testset "Final size bounds - row sums" begin
        K = [1.8 0.4; 0.6 1.2]
        π = [0.3, 0.7]
        bounds = final_size_bounds(K, π; info=:row)

        τ_actual = solve_final_size_vector(K, π)

        # Lower bounds are 0 (always valid)
        @test all(bounds.lower .== 0.0)

        # Upper bounds must contain actual
        for i in 1:2
            @test bounds.upper[i] ≥ τ_actual[i] - 1e-8
        end
    end

    @testset "Total final size bounds" begin
        K = [2.0 0.3; 0.5 1.5]
        π = [0.4, 0.6]
        τ_actual = solve_final_size_vector(K, π)
        τ_bar_actual = sum(π .* τ_actual)

        # Row sums
        b_row = total_final_size_bounds(K, π; info=:row)
        @test b_row.lower ≤ τ_bar_actual + 1e-8
        @test b_row.upper ≥ τ_bar_actual - 1e-8

        # Column sums
        b_col = total_final_size_bounds(K, π; info=:col)
        @test b_col.lower ≤ τ_bar_actual + 1e-8
        @test b_col.upper ≥ τ_bar_actual - 1e-8

        # Subcritical → all zero
        K_sub = [0.5 0.2; 0.1 0.4]
        b_sub = total_final_size_bounds(K_sub, π; info=:row)
        @test b_sub.lower == 0.0
        @test b_sub.upper == 0.0
    end

    @testset "info= argument validation" begin
        K = [1.5 0.5; 0.3 1.7]
        π = [0.5, 0.5]

        @test_throws ArgumentError r0_bounds(K; info=:bogus)
        @test_throws ArgumentError r0_bounds_detailed_balance(K, π; info=:bogus)
        @test_throws ArgumentError final_size_bounds(K, π; info=:bogus)
        @test_throws ArgumentError total_final_size_bounds(K, π; info=:bogus)

        # :both is valid for the R₀ bounds but NOT for the final-size bounds, whose
        # bodies are `if info == :col ... else <row>` — a widened guard would silently
        # return row bounds. This is the case that pins the two allowed sets apart.
        @test_throws ArgumentError final_size_bounds(K, π; info=:both)
        @test_throws ArgumentError total_final_size_bounds(K, π; info=:both)

        # The :both branch of the detailed-balance bounds. Intersecting the row and column
        # intervals is only meaningful when detailed balance actually holds — a precondition
        # the function does not check — so the pair is chosen to satisfy it exactly and the
        # property is asserted rather than assumed. Non-uniform π keeps the row and column
        # intervals distinct, so returning either one under :both fails here.
        K_db = [1.5 0.75; 0.25 1.25]
        π_db = [0.25, 0.75]
        @test π_db[1] * K_db[1, 2] == π_db[2] * K_db[2, 1]
        b_both = r0_bounds_detailed_balance(K_db, π_db; info=:both)
        @test b_both.lower ≤ b_both.upper
        @test b_both.lower ≤ maximum(abs.(eigvals(K_db))) ≤ b_both.upper
        @test b_both.lower > r0_bounds_detailed_balance(K_db, π_db; info=:row).lower
    end

    @testset "ContactMatrix convenience" begin
        part = AgePartition([0, 18, 65])
        pop = [11000.0, 33000.0, 9500.0]
        M = [7.0 2.5 1.0; 2.0 8.0 2.0; 0.5 2.0 4.0]
        cm = ContactMatrix(M, part, pop)

        b = r0_bounds(cm; info=:row)
        @test b.lower > 0
        @test b.upper ≥ b.lower

        fs = final_size_bounds(cm; info=:col)
        @test length(fs.lower) == 3
        @test all(fs.lower .≥ 0)
        @test all(fs.upper .≤ 1.0)

        tfs = total_final_size_bounds(cm; info=:row)
        @test 0 ≤ tfs.lower ≤ tfs.upper ≤ 1.0
    end

    @testset "Epidemic uncertainty over fiber" begin
        # Create a small set of contact matrices
        part = AgePartition([0, 50])
        pop = [5000.0, 5000.0]
        cms = [
            ContactMatrix([6.0 2.0; 2.0 5.0], part, pop),
            ContactMatrix([7.0 1.5; 1.5 4.5], part, pop),
            ContactMatrix([5.0 3.0; 3.0 6.0], part, pop),
        ]
        result = epidemic_uncertainty(cms)
        @test result.r0.lower ≤ result.r0.upper
        @test result.final_size.lower ≤ result.final_size.upper
        @test result.r0.lower > 0
        @test result.final_size.lower ≥ 0
    end
end

@testset "Review fixes regression" begin
    @testset "Duplicate ProductPartition dimensions rejected" begin
        age = AgePartition([0, 18])
        @test_throws ArgumentError ProductPartition(age, age)
    end

    @testset "ContactSurvey ID validation" begin
        # Contacts referencing unknown participant
        @test_throws ArgumentError ContactSurvey(
            DataFrame(part_id=[1, 2], part_age=[5.0, 25.0]),
            DataFrame(part_id=[1, 99], cnt_age=[10.0, 20.0]),
        )
        # Duplicate participant IDs
        @test_throws ArgumentError ContactSurvey(
            DataFrame(part_id=[1, 1], part_age=[5.0, 25.0]),
            DataFrame(part_id=[1], cnt_age=[10.0]),
        )
        # Missing participant ID
        @test_throws ArgumentError ContactSurvey(
            DataFrame(part_id=[1, missing], part_age=[5.0, 25.0]),
            DataFrame(part_id=[1], cnt_age=[10.0]),
        )
    end

    @testset "q-parameter validation" begin
        @test_throws ArgumentError BlockAssortativityParams(q=Dict(:sep => 1.5))
        @test_throws ArgumentError BlockAssortativityParams(q=Dict(:sep => -1.1))
        @test_throws ArgumentError BlockAssortativityParams(q=Dict(:sep => NaN))
        @test_throws ArgumentError BlockAssortativityParams(q=Dict(:sep => Inf))

        # The positional form validates too, whatever the dict's element type,
        # so no construction path can hand the solver a q outside [-1, 1].
        @test_throws ArgumentError BlockAssortativityParams(Dict(:sep => 5.0))
        @test_throws ArgumentError BlockAssortativityParams(Dict(:sep => 5))
        @test_throws ArgumentError BlockAssortativityParams(Dict(:sep => 3 // 2))
        @test_throws ArgumentError BlockAssortativityParams(Dict(:sep => -Inf))

        # The endpoints themselves are valid, and integer input is converted.
        @test BlockAssortativityParams(Dict(:sep => 1.0)).q[:sep] == 1.0
        @test BlockAssortativityParams(Dict(:sep => -1)).q[:sep] === -1.0
        @test isempty(BlockAssortativityParams().q)
    end

    @testset "NGM has no population rescaling; bounds bracket true R₀" begin
        part = AgePartition([0, 18, 65])
        pop = [11000.0, 33000.0, 9500.0]
        totals = [3.0e4 2.0e4 5.0e3; 2.0e4 5.0e4 1.5e4; 5.0e3 1.5e4 2.0e4]
        M = [totals[i, j] / pop[j] for i in 1:3, j in 1:3]
        cm = ContactMatrix(M, part, pop)
        # Reciprocal by construction (symmetric total contacts)
        @test M * Diagonal(pop) ≈ (M * Diagonal(pop))'
        # NGM is just scaled mean contacts — no N[i]/N[j] factor
        K = next_generation_matrix(cm; transmissibility=0.4, recovery_rate=0.5)
        @test K ≈ (0.4 / 0.5) .* M
        R0_true = R₀(cm; transmissibility=0.4, recovery_rate=0.5)
        # Detailed-balance bounds must bracket R₀ for a reciprocal matrix
        db = r0_bounds_detailed_balance(cm; info=:row, transmissibility=0.4, recovery_rate=0.5)
        @test db.lower ≤ R0_true + 1e-8
        @test db.upper ≥ R0_true - 1e-8
        # General row/col/both bounds also bracket
        for s in (:row, :col, :both)
            b = r0_bounds(cm; info=s, transmissibility=0.4, recovery_rate=0.5)
            @test b.lower ≤ R0_true + 1e-8
            @test b.upper ≥ R0_true - 1e-8
        end
    end

    @testset "q-parameter sampler bounds restricted to [-1,1]" begin
        using Random
        age = AgePartition([0, 18]; labels=["child", "adult"])
        sep = CategoricalPartition(:sep; participant_col=:part_sep,
            contact_col=:cnt_sep, levels=["low", "high"])
        prod = age × sep
        interm = SourceStratifiedContactMatrix(
            [3.0 2.0 1.0 1.5; 1.5 1.0 2.5 2.0], age, prod, [100.0, 100.0, 150.0, 150.0])
        base_pop = [200.0, 300.0]
        base_counts = [5.0 2.5; 2.5 4.5] * Diagonal(base_pop)
        base_M = ((base_counts + base_counts') / 2) * Diagonal(1.0 ./ base_pop)
        base_cm = ContactMatrix(base_M, age, base_pop)
        spec = ConstrainedGeneralizedLift(interm; source_map=PartitionMap(prod, age))
        rng = Random.MersenneTwister(1)
        @test_throws ArgumentError sample_constrained_lifts(base_cm, spec, 1; bounds=(-2.0, 2.0), rng=rng)
        @test_throws ArgumentError QParameterSpace(base_cm, spec; bounds=(-2.0, 2.0))
        @test_throws ArgumentError sample_perblock_lifts(base_cm, spec, 1; bounds=(0.0, 1.5), rng=rng)
        @test_throws ArgumentError mcmc_constrained_lifts(base_cm, spec, 1; bounds=(-1.5, 1.5), rng=rng)
        # In-range bounds still work
        @test QParameterSpace(base_cm, spec; bounds=(-1.0, 1.0)) isa QParameterSpace
    end

    @testset "Generalized lift requires reciprocal base" begin
        age = AgePartition([0, 18])
        ses = CategoricalPartition(:ses; levels=["a", "b"])
        nonrecip = ContactMatrix([2.0 0.5; 0.9 1.5], age, [100.0, 200.0])
        spec = GeneralizedLift(ses; distribution=[0.5, 0.5])
        @test_throws ArgumentError generalize(nonrecip, spec)
        # Reciprocal base succeeds
        recip = ContactMatrix([2.0 0.5; 1.0 1.5], age, [100.0, 200.0])
        @test generalize(recip, spec) isa ContactMatrix
    end

    @testset "migrate_coarsen rejects skip_invalid surveys with unset groups" begin
        survey = ContactSurvey(
            DataFrame(part_id=[1, 2], part_age=[5.0, missing]),
            DataFrame(part_id=[1], cnt_age=[10.0]),
        )
        fine = AgePartition([0, 18, 65])
        coarse = AgePartition([0, 65])
        acs = ContactSurveyACSet(survey, fine; skip_invalid=true)
        @test_throws ArgumentError migrate_coarsen(acs, PartitionMap(fine, coarse))
    end

    @testset "compose_uwd validates wiring" begin
        sharers = Dict(
            :home => ContactSharer([1.0 0.0; 0.0 1.0]),
            :work => ContactSharer([0.5 0.5; 0.5 0.5]),
        )
        good = @relation (age,) begin
            home(age)
            work(age)
        end
        @test compose_uwd(good, sharers) ≈ [1.5 0.5; 0.5 1.5]
        # A box wired to an internal junction is not additive composition
        internal = @relation (age,) begin
            home(age)
            work(other)
        end
        @test_throws ArgumentError compose_uwd(internal, sharers)
    end

    @testset "q-lift marginal preservation" begin
        age = AgePartition([0, 18]; labels=["child", "adult"])
        sep = CategoricalPartition(:sep; participant_col=:part_sep,
            contact_col=:cnt_sep, levels=["low", "high"])
        prod = age × sep
        prod_pop = [100.0, 100.0, 150.0, 150.0]
        interm_M = [3.0 2.0 1.0 1.5; 1.5 1.0 2.5 2.0]
        interm = SourceStratifiedContactMatrix(interm_M, age, prod, prod_pop)
        base_pop = [200.0, 300.0]
        base_M = [5.0 2.5; 2.5 4.5]
        base_counts = base_M * Diagonal(base_pop)
        base_counts_sym = (base_counts + base_counts') / 2
        base_M_sym = base_counts_sym * Diagonal(1.0 ./ base_pop)
        base_cm = ContactMatrix(base_M_sym, age, base_pop)
        source_to_age = PartitionMap(prod, age)
        spec = ConstrainedGeneralizedLift(interm; source_map=source_to_age)

        for q_val in [0.3, 0.5, -0.3]
            params = BlockAssortativityParams(q=Dict(:sep => q_val))
            pspec = ParameterizedConstrainedLift(spec; default_params=params)
            result = constrained_generalize(base_cm, pspec)
            C = matrix(result) * Diagonal(population(result))
            # Total contacts must be symmetric (reciprocity)
            @test C ≈ C' atol=1e-10
            # Coarsening back to base must be preserved
            @test matrix(result ↓ age) ≈ matrix(base_cm) atol=1e-10
        end
    end

    @testset "Refinement population mismatch rejected" begin
        p = AgePartition([0, 18, 65])
        M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
        pop = [1000.0, 3000.0, 500.0]
        cm = ContactMatrix(M, p, pop)
        fine = AgePartition([0, 5, 18, 30, 65])
        # 5 fine groups: [0,5), [5,18), [18,30), [30,65), 65+
        # coarse pop is [1000, 3000, 500] but this doesn't match aggregated fine pop
        bad_pop = [400.0, 500.0, 1500.0, 1200.0, 600.0]
        @test_throws ArgumentError refine(cm, fine, bad_pop)
    end

    @testset "Epidemic bounds π validation" begin
        K = [2.0 0.5; 0.5 1.5]
        # π doesn't sum to 1
        @test_throws ArgumentError solve_final_size_vector(K, [0.3, 0.3])
        # π has negative values
        @test_throws ArgumentError solve_final_size_vector(K, [0.6, -0.4])
        # π has NaN
        @test_throws ArgumentError solve_final_size_vector(K, [0.5, NaN])
        # Valid case still works
        τ = solve_final_size_vector(K, [0.6, 0.4])
        @test all(τ .>= 0)
        @test all(τ .<= 1)
    end
end

# ---------------------------------------------------------------------------
# Build/coarsen equivariance for activity-vector lifts.
#
# An activity vector `a` (mean contacts per person, per group) can be turned
# into a ContactMatrix, and a ContactMatrix can be coarsened along a fibre map
# that merges fine groups into coarse ones. These two operations should commute:
# coarsening the activity vector and then building should give the same matrix
# as building and then coarsening. If they don't, a model fitted at one spatial
# or demographic resolution can't be interpreted at another.
#
# This is distinct from the coarsen-functoriality test above ("∘ (AgeMap
# composition)"), which composes two coarsenings of the *same* matrix. Here we
# relate two different operations — build and coarsen — so the result depends
# on which lift is used, not just on coarsening itself.
#
# Three findings pinned here:
#   1. The reciprocal lift M[i,j] = a_i*N_i*a_j/D  (D = Σ_k a_k*N_k) commutes
#      unconditionally — any surjective fibre map, contiguous or not, any
#      populations. Both paths reduce to S_I*S_J/(D*N_J) with S_I = Σ_{i∈I} a_i*N_i.
#   2. The naive lift M[i,j] = a_i*a_j/D (missing the N_i factor) is not a valid
#      MeanContacts ContactMatrix: MeanContacts reciprocity requires symmetry of
#      *total* contacts, M[i,j]*N_j == M[j,i]*N_i, which fails whenever group
#      populations differ. So it is a trap under MeanContacts, not a
#      counterexample to (1). The N_i factor is required by ContACT's
#      convention that column = participant (per-capita) and row = contactee
#      (extensive).
#   3. The Britton-Ball assortative/disassortative kernels allocate contacts by
#      walking activity strata in order, so they commute with coarsening only
#      when the fibre map is order-preserving (each coarse group is a contiguous
#      run of fine strata). Merging non-adjacent strata changes the order the
#      algorithm depends on and the mismatch is O(1) in mean-contacts units, not
#      floating-point noise. Normal usage merges adjacent strata, so this is a
#      documented gap rather than a blocker.
#
# ---------------------------------------------------------------------------
@testset "activity lift: build/coarsen equivariance" begin
    function pop_weighted_mean_pushforward(a, N, assignments, n_coarse)
        num = zeros(n_coarse); pop_c = zeros(n_coarse)
        for i in eachindex(a)
            I = assignments[i]
            num[I] += N[i] * a[i]
            pop_c[I] += N[i]
        end
        num ./ pop_c, pop_c
    end

    a = [0.5, 1.5, 2.0, 1.0]
    N = [100.0, 30.0, 60.0, 200.0]
    fine = CategoricalPartition{:activity,Int}(collect(1:4))
    coarse = CategoricalPartition{:activity,Int}(collect(1:2))

    @testset "reciprocal lift commutes for any fibre map" begin
        for assignments in ([1, 1, 2, 2], [1, 2, 1, 2])   # contiguous, non-contiguous
            f = PartitionMap(fine, coarse, assignments)
            cm_fine = proportionate_mixing(a, N, fine)

            # sanity: this construction is already reciprocal (symmetrise is a no-op)
            @test matrix(↔(cm_fine)) ≈ matrix(cm_fine)

            cm_coarse_build_then_coarsen = coarsen(cm_fine, f)

            a_c, N_c = pop_weighted_mean_pushforward(a, N, assignments, 2)
            cm_coarse_coarsen_then_build = proportionate_mixing(a_c, N_c, coarse)

            @test matrix(cm_coarse_build_then_coarsen) ≈
                  matrix(cm_coarse_coarsen_then_build) atol=1e-10
            @test population(cm_coarse_build_then_coarsen) ≈ N_c
        end
    end

    @testset "naive lift without N_i factor is not valid MeanContacts" begin
        # Same D = Σ a_k*N_k as the reciprocal lift, so the N_i factor in the
        # numerator is the only thing that differs between the two.
        D = sum(a .* N)
        M_naive = (a * transpose(a)) ./ D
        cm_naive = ContactMatrix(M_naive, fine, N, MeanContacts())
        # MeanContacts reciprocity requires M[i,j]*N_j == M[j,i]*N_i; fails
        # whenever populations differ, so this is not a valid MeanContacts lift.
        recip_err = maximum(abs.(M_naive .* transpose(N) .- transpose(M_naive .* transpose(N))))
        @test recip_err > 0.5   # measured ≈0.67 on this case — not a rounding artifact
    end

    @testset "assortative/disassortative kernels need an order-preserving fibre map" begin
        row = [5.0, 15.0, 20.0, 10.0]
        col = [8.0, 12.0, 18.0, 12.0]
        contiguous = [1, 1, 2, 2]
        noncontig = [1, 2, 1, 2]

        fiber_sum(M, assignments, n_coarse) = begin
            out = zeros(n_coarse, n_coarse)
            for i in eachindex(assignments), j in eachindex(assignments)
                out[assignments[i], assignments[j]] += M[i, j]
            end
            out
        end
        coarsen_marginal(v, assignments, n_coarse) = begin
            out = zeros(n_coarse)
            for i in eachindex(assignments)
                out[assignments[i]] += v[i]
            end
            out
        end

        for kernel in (:assortative, :disassortative)
            plan_fine = activity_mixing_plan(row, col, kernel)

            # Contiguous fibre map: coarsen-then-build == build-then-coarsen.
            row_c = coarsen_marginal(row, contiguous, 2)
            col_c = coarsen_marginal(col, contiguous, 2)
            @test activity_mixing_plan(row_c, col_c, kernel) ≈
                  fiber_sum(plan_fine, contiguous, 2)

            # Non-contiguous fibre map: genuine counterexample, not noise —
            # order-based transport doesn't commute with reordering the strata.
            row_nc = coarsen_marginal(row, noncontig, 2)
            col_nc = coarsen_marginal(col, noncontig, 2)
            mismatch = maximum(abs.(
                activity_mixing_plan(row_nc, col_nc, kernel) .-
                fiber_sum(plan_fine, noncontig, 2)))
            @test mismatch > 1e-6
        end
    end
end

# ---------------------------------------------------------------------------
# Reconnect UK survey: cross-validation against R reference values
# ---------------------------------------------------------------------------
# These tests use pre-computed R reference matrices (simple day-weighted mean)
# committed as CSV fixtures. They validate that ContACT.jl's compute_matrix
# functor produces identical results when given the same data and weights.
#
# Full-data tests require downloading Zenodo data (~13 MB).
# Run with: CONTACT_RUN_RECONNECT=1 julia --project=. -e 'using Pkg; Pkg.test()'

if get(ENV, "CONTACT_RUN_RECONNECT", "") == "1"

using CSV
using Downloads: download

@testset "Reconnect UK cross-validation" begin
    # Load reference matrices from fixtures
    fixture_dir = joinpath(@__DIR__, "fixtures", "reconnect")

    function load_ref_matrix(filename)
        df = CSV.read(joinpath(fixture_dir, filename), DataFrame)
        Matrix{Float64}(df[:, 2:end])
    end

    ref_age = load_ref_matrix("reconnect_age_matrix_R.csv")
    ref_eth = load_ref_matrix("reconnect_ethnicity_matrix_R.csv")
    ref_coarse = load_ref_matrix("reconnect_age_coarse_matrix_R.csv")
    ref_home = load_ref_matrix("reconnect_age_home_matrix_R.csv")
    ref_work = load_ref_matrix("reconnect_age_work_matrix_R.csv")
    ref_school = load_ref_matrix("reconnect_age_school_matrix_R.csv")
    ref_other = load_ref_matrix("reconnect_age_other_matrix_R.csv")

    # Download Reconnect data from Zenodo
    zenodo_base = "https://zenodo.org/records/17339866/files"
    data_dir = mktempdir()

    filenames = [
        "reconnect_participant_common.csv",
        "reconnect_participant_extra.csv",
        "reconnect_contact_common.csv",
        "reconnect_contact_extra.csv",
        "reconnect_sday.csv",
    ]
    for f in filenames
        dest = joinpath(data_dir, f)
        isfile(dest) || download("$zenodo_base/$f", dest)
    end

    # Load and join data
    part_common = CSV.read(joinpath(data_dir, "reconnect_participant_common.csv"), DataFrame; missingstring="NA")
    part_extra = CSV.read(joinpath(data_dir, "reconnect_participant_extra.csv"), DataFrame; missingstring="NA")
    sday = CSV.read(joinpath(data_dir, "reconnect_sday.csv"), DataFrame; missingstring="NA")
    cnt_common = CSV.read(joinpath(data_dir, "reconnect_contact_common.csv"), DataFrame; missingstring="NA")
    cnt_extra = CSV.read(joinpath(data_dir, "reconnect_contact_extra.csv"), DataFrame; missingstring="NA")

    participants = innerjoin(part_common, part_extra, on = :part_id)
    participants = innerjoin(participants, sday, on = :part_id)

    # Drop participants with missing dayofweek
    dropmissing!(participants, :dayofweek)

    # Day type & weights (dayofweek: 0=Sunday, 6=Saturday)
    participants.day_type = ifelse.(
        participants.dayofweek .∈ Ref([0, 6]),
        "weekend", "weekday"
    )
    day_counts = combine(groupby(participants, :day_type), nrow => :n)
    total_n = nrow(participants)
    day_counts.sample_prop = day_counts.n ./ total_n
    day_counts.target_prop = ifelse.(day_counts.day_type .== "weekday", 5/7, 2/7)
    day_counts.day_weight = day_counts.target_prop ./ day_counts.sample_prop
    participants = innerjoin(participants, select(day_counts, :day_type, :day_weight), on = :day_type)

    # Contacts
    contacts = innerjoin(cnt_common, cnt_extra, on = [:cont_id, :part_id])
    rename!(participants, :part_age_exact => :part_age)
    participants.part_age = Float64.(participants.part_age)
    rename!(contacts, :cnt_age_exact => :cnt_age)
    dropmissing!(contacts, :cnt_age)
    contacts.cnt_age = Float64.(collect(contacts.cnt_age))
    valid_ids = Set(participants.part_id)
    filter!(row -> row.part_id in valid_ids, contacts)

    survey = ContactSurvey(participants, contacts)

    # -----------------------------------------------------------------------
    # Age partition (16 groups)
    # -----------------------------------------------------------------------
    age_breaks = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
    age_part = AgePartition(age_breaks)

    @testset "Age matrix matches R reference" begin
        cm = compute_matrix(survey, age_part; weights = :day_weight)
        # R convention: rows are reporter age → transpose of ContACT (cols = reporter)
        # R computes M[reporter, contact] while ContACT has M[contact, reporter]
        # The R CSV has rownames = reporter → ref_age[i,j] = mean contacts of
        # reporter group i with contact group j
        # ContACT: matrix(cm)[i,j] = mean contacts of reporter group j with contact group i
        # So ContACT matrix = transpose(R matrix)
        M_julia = matrix(cm)
        @test size(M_julia) == (16, 16)
        @test isapprox(M_julia, ref_age', atol = 1e-4)
    end

    @testset "Setting-specific matrices match R reference" begin
        settings_ref = Dict(
            "Home" => ref_home,
            "Work" => ref_work,
            "School" => ref_school,
            "Other" => ref_other,
        )
        for (setting, ref) in settings_ref
            cnt_s = filter(row -> row.cnt_location == setting, contacts)
            survey_s = ContactSurvey(participants, cnt_s)
            cm_s = compute_matrix(survey_s, age_part; weights = :day_weight)
            @test isapprox(matrix(cm_s), ref', atol = 1e-4)
        end
    end

    @testset "Setting composition holds exactly" begin
        cms_settings = ContactMatrix[]
        for setting in ["Home", "Work", "School", "Other"]
            cnt_s = filter(row -> row.cnt_location == setting, contacts)
            survey_s = ContactSurvey(participants, cnt_s)
            push!(cms_settings, compute_matrix(survey_s, age_part; weights = :day_weight))
        end
        cm_total = compute_matrix(survey, age_part; weights = :day_weight)
        cm_composed = cms_settings[1] ⊕ cms_settings[2] ⊕ cms_settings[3] ⊕ cms_settings[4]
        @test isapprox(matrix(cm_composed), matrix(cm_total), atol = 1e-10)
    end

    @testset "Ethnicity matrix matches R reference" begin
        eth_levels = ["Asian", "Black", "Mixed", "Other", "White"]
        eth_part = CategoricalPartition(
            :ethnicity, eth_levels;
            participant_col = :part_ethnicity,
            contact_col = :cnt_ethnicity,
        )
        participants_eth = filter(
            row -> !ismissing(row.part_ethnicity) && row.part_ethnicity in eth_levels,
            participants,
        )
        valid_eth_ids = Set(participants_eth.part_id)
        contacts_eth = filter(
            row -> row.part_id in valid_eth_ids &&
                   !ismissing(row.cnt_ethnicity) &&
                   row.cnt_ethnicity in eth_levels,
            contacts,
        )
        survey_eth = ContactSurvey(participants_eth, contacts_eth)
        cm_eth = compute_matrix(survey_eth, eth_part; weights = :day_weight)
        @test isapprox(matrix(cm_eth), ref_eth', atol = 1e-4)
    end

    @testset "Coarse 3×3 matrix matches R reference (direct)" begin
        coarse_part = AgePartition([0, 18, 65])
        cm_coarse = compute_matrix(survey, coarse_part; weights = :day_weight)
        @test isapprox(matrix(cm_coarse), ref_coarse', atol = 1e-4)
    end

    @testset "Symmetrisation is reciprocal and idempotent" begin
        cm = compute_matrix(survey, age_part; weights = :day_weight)
        cm_sym = ↔(cm)
        M_s = matrix(cm_sym)
        N = population(cm_sym)
        n = n_groups(cm_sym)
        # Reciprocity
        for i in 1:n, j in 1:n
            @test isapprox(M_s[i, j] * N[j], M_s[j, i] * N[i], atol = 1e-10)
        end
        # Idempotence
        cm_sym2 = ↔(cm_sym)
        @test isapprox(matrix(cm_sym2), M_s, atol = 1e-14)
    end

    @testset "NGM and R₀ are well-defined" begin
        cm = compute_matrix(survey, age_part; weights = :day_weight)
        cm_sym = ↔(cm)
        τ_val = 0.05
        K = next_generation_matrix(cm_sym; transmissibility = τ_val)
        R0_val = R₀(cm_sym; transmissibility = τ_val)
        @test R0_val > 0
        @test isapprox(R0_val, maximum(abs.(eigvals(K))), rtol = 1e-10)
        # With only ~3.8 known-age contacts and τ=0.05, R₀ is modest
        @test 0.01 < R0_val < 5.0
    end

    @testset "Epidemic bounds bracket true R₀" begin
        cm = compute_matrix(survey, age_part; weights = :day_weight)
        cm_sym = ↔(cm)
        τ_val = 0.05
        K = next_generation_matrix(cm_sym; transmissibility = τ_val)
        R0_val = R₀(cm_sym; transmissibility = τ_val)
        bounds = r0_bounds(K)
        @test bounds.lower ≤ R0_val + 1e-10
        @test R0_val ≤ bounds.upper + 1e-10
    end
end

end # CONTACT_RUN_RECONNECT

end # top-level testset
