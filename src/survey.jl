"""
Survey operations: filtering and subsetting morphisms.

In the categorical framework, these are morphisms in the restricted
survey category that preserve the schema (participants + contacts + link).
"""

"""
    filter_survey(survey::ContactSurvey; kwargs...)

Filter a survey by participant attributes. Returns a new `ContactSurvey`
with only matching participants and their contacts.

# Keyword Arguments
Any column name in `participants` can be used as a filter. Values can be:
- A single value: exact match
- A vector: membership test
- A function: predicate applied to column values

# Examples
```julia
# Filter by country
uk_survey = filter_survey(survey, country="UK")

# Filter by multiple countries
euro_survey = filter_survey(survey, country=["UK", "DE", "FR"])

# Filter by age range
adult_survey = filter_survey(survey, part_age=age -> age >= 18)
```
"""
function filter_survey(survey::ContactSurvey; kwargs...)
    mask = trues(nrow(survey.participants))
    for (col, val) in kwargs
        col in propertynames(survey.participants) ||
            throw(ArgumentError("column :$col not found in participants"))
        if val isa Function
            mask .&= val.(survey.participants[!, col])
        elseif val isa AbstractVector
            mask .&= survey.participants[!, col] .∈ Ref(val)
        else
            mask .&= survey.participants[!, col] .== val
        end
    end

    filtered_parts = survey.participants[mask, :]
    part_ids = Set(filtered_parts.part_id)
    filtered_contacts = filter(row -> row.part_id in part_ids, survey.contacts)

    ContactSurvey(filtered_parts, filtered_contacts; metadata=copy(survey.metadata))
end

"""
    subset_survey(survey::ContactSurvey, part_ids)

Subset a survey to specific participant IDs. Returns a new `ContactSurvey`.
"""
function subset_survey(survey::ContactSurvey, part_ids::AbstractVector)
    id_set = Set(part_ids)
    parts = filter(row -> row.part_id in id_set, survey.participants)
    contacts = filter(row -> row.part_id in id_set, survey.contacts)
    ContactSurvey(parts, contacts; metadata=copy(survey.metadata))
end
