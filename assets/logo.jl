# The julia-style colors are selected from https://github.com/JuliaGraphics/julia-logo-graphics

using Luxor

R = 160.0
r = 120.0
θ = -60.0

pt1 = Point(R * cosd(90.0), R * sind(90.0))
pt2 = Point(R * cosd(210.0), R * sind(210.0))
pt3 = Point(R * cosd(330.0), R * sind(330.0))

pt21 = Point(pt2.x - r * sind(θ), pt2.y - r * cosd(θ))
pt22 = Point(pt2.x + r * sind(θ), pt2.y + r * cosd(θ))

pt31 = Point(pt3.x - r * sind(-θ), pt3.y - r * cosd(-θ))
pt32 = Point(pt3.x + r * sind(-θ), pt3.y + r * cosd(-θ))

@svg begin
    rotate(π)

    sethue(0.22, 0.596, 0.149)
    circle(pt1, r, :stroke)
    pie(pt1.x, pt1.y, r, π/2, π, :fill)
    pie(pt1.x, pt1.y, r, 3/2 * π, 2π, :fill)

    sethue(0.584, 0.345, 0.698)
    circle(pt2, r, :stroke)
    arc2sagitta(pt21, pt22, r/2, :fillpreserve)
    arc2sagitta(pt22, pt21, r/2, :fillpreserve)

    sethue(0.796, 0.235, 0.2)
    circle(pt3, r, :fill)
    sethue("white")
    arc2sagitta(pt31, pt32, r/2, :fillpreserve)
    arc2sagitta(pt32, pt31, r/2, :fillpreserve)

    settext(
        "<span font='33' background='royalblue' foreground='white'> JuEQ.jl </span>",
        pt3,
        halign="center",
        valign="center",
        markup=true,
        angle = -30,
    )
end 800 800 "logo"
