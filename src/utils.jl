function second_moment_of_area_closed_ended_pile(D::Real)
    π * D^4  / 64
end

function second_moment_of_area_open_ended_pile(D_outer::Real, D_inner::Real)
    @assert D_outer > D_inner
    π * (D_outer^4 - D_inner^4) / 64
end
