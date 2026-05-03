using Test
using ContACT
using DataFrames
using LinearAlgebra
using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Programs: @relation

# Resolve ambiguity: use ContACT's operators explicitly
import ContACT: ⊕, ⊗, ↓, ↑, ▷, ρ

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

    # Custom labels
    p2 = AgePartition([0, 5, 18]; labels=["child", "youth", "adult"])
    @test age_labels(p2) == ["child", "youth", "adult"]

    # Validation
    @test_throws ArgumentError AgePartition(Int[])
    @test_throws ArgumentError AgePartition([5, 5, 10])
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

    # Block structure check: diagonal blocks should have coupling[i,i] * M
    M_strat = matrix(cm_strat)
    @test M_strat[1:3, 1:3] ≈ 0.7 .* M
    @test M_strat[4:6, 1:3] ≈ 0.3 .* M
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

@testset "Utilities" begin
    p = AgePartition([0, 18, 65])
    pop = [1000.0, 3000.0, 500.0]
    M = [2.0 1.0 0.5; 1.0 3.0 1.0; 0.5 1.0 1.5]
    cm = ContactMatrix(M, p, pop)

    # Per capita conversion
    cm_pc = to_per_capita(cm)
    @test cm_pc.semantics isa PerCapitaRate

    # Counts conversion
    cm_counts = to_counts(cm)
    @test cm_counts.semantics isa ContactCounts

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

end # top-level testset
