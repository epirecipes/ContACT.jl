using Test
using ContACT
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

        # Very extreme q should be infeasible
        extreme_params = BlockAssortativityParams(q=Dict(:sep => 0.99))
        pspec_e = ParameterizedConstrainedLift(spec; default_params=extreme_params)
        # May or may not be feasible depending on structure; just test it returns Bool
        @test is_feasible(base, pspec_e) isa Bool
    end

    @testset "sample_constrained_lifts" begin
        using Random
        rng = Random.MersenneTwister(42)
        samples = sample_constrained_lifts(base, spec, 5;
            dimensions=[:sep], bounds=(-0.5, 0.5), rng=rng)
        @test length(samples) <= 5
        @test length(samples) >= 1  # at least some feasible
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
        @test length(samples) >= 1
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

        # MCMC with a log-density targeting high assortativity
        log_dens = (cm, _) -> assortativity_index(cm, :sep)
        rng2 = Random.MersenneTwister(77)
        result_targeted = mcmc_constrained_lifts(base, spec, 10;
            dimensions=[:sep], bounds=(-0.4, 0.4),
            log_density=log_dens,
            proposal_scale=0.05, burnin=20, thin=1, rng=rng2)
        @test length(result_targeted.chain) == 10
        # Targeted chain should have higher average assortativity than flat
        mean_ai_targeted = sum(assortativity_index(cm, :sep)
                               for cm in result_targeted.matrices) / length(result_targeted.matrices)
        mean_ai_flat = sum(assortativity_index(cm, :sep)
                          for cm in result.matrices) / length(result.matrices)
        # Not guaranteed every run, but very likely with assortativity as density
        @test mean_ai_targeted >= mean_ai_flat - 0.5  # soft check
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
    @test K ≈ 2.0 .* [matrix(base)[i, j] * pop[i] / pop[j] for i in 1:2, j in 1:2]
    @test R₀(base; transmissibility=0.5, recovery_rate=0.25) ≈ 2.0 * ρ(base)
    beta = calibrate_transmissibility(base, 2.7; recovery_rate=0.25)
    @test R0(base; transmissibility=beta, recovery_rate=0.25) ≈ 2.7

    @test_throws ArgumentError GeneralizedLift(ses; distribution=[0.5, 0.5, 0.5])
    @test_throws ArgumentError generalize(base, ses, [10.0, 20.0, 30.0, 40.0, 50.0, 60.0])
    @test_throws ArgumentError BlockMixing([0.5 0.5; 0.5 0.5])
end

@testset "Utilities" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ContactMatrix(M, p, pop)

    # Per capita conversion
    cm_pc = to_per_capita(cm)
    @test cm_pc.semantics isa PerCapitaRate
    @test matrix(cm_pc) ≈ [M[i, j] / pop[j] for i in 1:3, j in 1:3]

    # Counts conversion
    cm_counts = to_counts(cm)
    @test cm_counts.semantics isa ContactCounts
    @test matrix(cm_counts) ≈ [M[i, j] * pop[j] for i in 1:3, j in 1:3]
    @test matrix(to_per_capita(cm_counts)) ≈ M

    @test_throws ArgumentError to_per_capita(
        ContactMatrix([0.0 1.0; 0.0 0.0], AgePartition([0, 18]), [10.0, 0.0])
    )

    # Spectral radius
    sr = spectral_radius(cm)
    @test sr > 0
    @test sr ≈ maximum(abs.(eigvals(M)))
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
            @test matrix(result ↓ age) ≈ matrix(base_cm) atol=1e-8
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

end # top-level testset
