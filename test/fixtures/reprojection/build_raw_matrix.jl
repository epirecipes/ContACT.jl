using ContACT
using CSV
using DataFrames

data_dir = joinpath(@__DIR__, "..", "..", "..", "data")

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

labels = group_labels(age_partition)
df = DataFrame(matrix(cm_age), Symbol.(labels))
insertcols!(df, 1, Symbol("") => labels)
CSV.write(joinpath(@__DIR__, "polymod_uk_age_raw_matrix.csv"), df)
