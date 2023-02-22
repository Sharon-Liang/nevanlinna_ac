using Test
using NevanlinnaAC
using LinearAlgebra
using FiniteDifferences

#make test fuctions
y(x) = sech(x)
dy(x) = -tanh(x) * sech(x)
d²y(x) = -sech(x)^3 + tanh(x)^2 * sech(x)
d³y(x) = 5*sech(x)^3 * tanh(x) - tanh(x)^3 * sech(x)

#mesh
xmax = 10
L = 2*xmax
nmesh = 100000

xmesh = [i for i in range(-xmax, xmax -L/nmesh, length=nmesh)]

ymesh = @. y(xmesh)

#derivatives
dy_exact = @. dy(xmesh)
dy_fft = fft_derivative(ymesh, L, 1)
dy_ngrad = map(x -> central_fdm(5,1)(y, x), xmesh)

@show(norm(@. dy_fft - dy_exact))
@show(norm(@. dy_ngrad - dy_exact))

d²y_exact = @. d²y(xmesh)
d²y_fft = fft_derivative(ymesh, L, 2)
d²y_ngrad = map(x -> central_fdm(5,2)(y, x), xmesh)
@show(norm(@. d²y_fft - d²y_exact))
@show(norm(@. d²y_ngrad - d²y_exact))

d³y_exact = @. d³y(xmesh)
d³y_fft = fft_derivative(ymesh, L, 3)
d³y_ngrad = map(x -> central_fdm(5, 3)(y, x), xmesh)
@show(norm(@. d³y_fft[1+50: 501-50] - d³y_exact[1+50: 501-50]))
@show(norm(@. d³y_ngrad - d³y_exact))

#plot
using Plots
plot(xmesh, ymesh)


plot(xmesh, dy_exact)
plot!(xmesh, dy_fft, line=(:dash))

plot(xmesh, d²y_exact)
plot!(xmesh, d²y_fft, line=(:dash))


plot(xmesh, d³y_exact)
plot!(xmesh, d³y_fft, line=(:dash))
