function runge_kutta(x0,xf,N,w,w0)
    h = (xf-x0)/N
    xr = zeros(N+1)
    xr[1] = x0
    wr = zeros(N+1)
    wr[1] = w0
    for i in 1:N
        k1 = h*w(xr[i],wr[i])
        k2 = h*w(xr[i]+h/2,wr[i] +k1/2)
        k3 = h*w(xr[i]+h/2,wr[i]+k2/2)
        k4 = h*w(xr[i]+h,wr[i]+k3/2)
        wr[i+1] = wr[i] + 1/6 * ( k1 + 2*k2 + 2*k3 + k4)
        xr[i+1] = xr[i] + h
    end
    return xr, wr

end