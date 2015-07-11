# performance benchmark

# singular value decomposition (SVD)
# @time svd(zeros(2000,2000));
# lstopo --no-legend --no-io cambridge.png

m = 500
M = zeros(m,m)

for c = 1:CPU_CORES
    blas_set_num_threads(c)

    s = 0.
    for i = 1:32
        t = @elapsed svd(M)

        if i > 2
            s += (t - s) / (i - 2)
        end
    end
    println(s)
end


