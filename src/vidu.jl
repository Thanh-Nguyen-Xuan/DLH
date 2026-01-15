using DifferentialEquations
using Plots
using LinearAlgebra
using Printf

# ------------------------------------------------------------
# 1. MÔ HÌNH CƠ BẢN: CÔNG TRÌNH SDOF VỚI DAMPER ĐIỆN TỪ
# ------------------------------------------------------------

struct ElectromagneticDamper
    # Thông số cơ học
    m::Float64      # Khối lượng công trình (kg)
    c::Float64      # Hệ số cản cơ học (Ns/m)
    k::Float64      # Độ cứng (N/m)
    
    # Thông số điện từ
    Ke::Float64     # Hệ số cảm ứng điện (Vs/m)
    Kt::Float64     # Hệ số lực điện từ (N/A)
    R::Float64      # Điện trở mạch (Ω)
    L::Float64      # Điện cảm (H)
    R_load::Float64 # Điện trở tải (Ω)
    
    # Vị trí lắp đặt
    position::Float64 # Chiều cao lắp đặt (0-1)
end

# ------------------------------------------------------------
# 2. HỆ PHƯƠNG TRÌNH VI PHÂN
# ------------------------------------------------------------

function coupled_equations!(du, u, p, t)
    # u = [x, v, i] - chuyển vị, vận tốc, dòng điện
    # du = [dx/dt, dv/dt, di/dt]
    
    x, v, i = u
    damper, F_external, wind_velocity = p
    
    # Tải trọng ngoài (gió/động đất)
    F_ext = F_external(t, x, v, wind_velocity(t))
    
    # Phương trình cơ học: m*dv/dt = -c*v - k*x - Kt*i + F_ext
    du[1] = v  # dx/dt = v
    du[2] = (F_ext - damper.c*v - damper.k*x - damper.Kt*i) / damper.m
    
    # Phương trình điện: L*di/dt = Ke*v - (R+R_load)*i
    du[3] = (damper.Ke*v - (damper.R + damper.R_load)*i) / damper.L
    
    return nothing
end

# ------------------------------------------------------------
# 3. CÁC MÔ HÌNH TẢI TRỌNG
# ------------------------------------------------------------

# 3.1 Tải trọng gió (mô hình Davenport)
function wind_load(t, x, v, wind_vel)
    # wind_vel: vận tốc gió tại thời điểm t
    ρ = 1.225        # Khối lượng riêng không khí (kg/m³)
    Cd = 1.2         # Hệ số cản
    A = 50.0         # Diện tích đón gió (m²)
    
    # Vận tốc gió tương đối (công trình có thể chuyển động)
    relative_vel = wind_vel - v
    
    # Lực gió (công thức đơn giản)
    F_wind = 0.5 * ρ * Cd * A * relative_vel * abs(relative_vel)
    
    return F_wind
end

# 3.2 Tải trọng động đất (sóng sin hoặc record thực)
function earthquake_load(t, x, v, ground_accel)
    # ground_accel: gia tốc nền tại thời điểm t
    return -damper.m * ground_accel(t)  # Lực quán tính
end

# 3.3 Vận tốc gió ngẫu nhiên (mô phỏng nhiễu)
function random_wind(t)
    # Thành phần trung bình + nhiễu
    V_mean = 20.0  # Vận tốc gió trung bình (m/s)
    V_turb = 5.0 * sin(0.5*t) + 2.0 * randn()  # Nhiễu
    return V_mean + V_turb
end

# 3.4 Gia tốc nền động đất (El Centro)
function el_centro_accel(t)
    # Dữ liệu El Centro 1940 (rút gọn)
    if t < 0.5
        return 0.0
    elseif t < 10.0
        return 0.3 * sin(2π*2*t) * exp(-0.1*t)  # Mô phỏng đơn giản
    else
        return 0.0
    end
end

# ------------------------------------------------------------
# 4. MÔ PHỎNG HỆ THỐNG
# ------------------------------------------------------------

function simulate_emd_system(;
        m=10000.0,      # Khối lượng công trình 10 tấn
        k=1e6,          # Độ cứng
        c=10000.0,      # Cản cơ học
        Ke=50.0,        # Hệ số cảm ứng
        Kt=50.0,        # Hệ số lực
        R=5.0,          # Điện trở cuộn dây
        L=0.1,          # Điện cảm
        R_load=10.0,    # Tải tiêu thụ
        simulation_time=30.0,
        load_type=:wind  # :wind hoặc :earthquake
    )
    
    # Tạo damper
    damper = ElectromagneticDamper(m, c, k, Ke, Kt, R, L, R_load, 0.7)
    
    # Chọn hàm tải trọng
    if load_type == :wind
        F_external = (t, x, v, w) -> wind_load(t, x, v, w)
        external_input = random_wind
    else  # earthquake
        F_external = (t, x, v, a) -> earthquake_load(t, x, v, a)
        external_input = el_centro_accel
    end
    
    # Điều kiện ban đầu
    u0 = [0.0, 0.0, 0.0]  # x=0, v=0, i=0
    
    # Tham số
    p = (damper, F_external, external_input)
    
    # Khoảng thời gian
    tspan = (0.0, simulation_time)
    
    # Giải hệ ODE
    prob = ODEProblem(coupled_equations!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=0.01)
    
    return sol, damper
end

# ------------------------------------------------------------
# 5. TÍNH TOÁN NĂNG LƯỢNG
# ------------------------------------------------------------

function calculate_energy(sol, damper)
    t = sol.t
    x = [sol[i][1] for i in 1:length(sol)]
    v = [sol[i][2] for i in 1:length(sol)]
    i = [sol[i][3] for i in 1:length(sol)]
    
    # Năng lượng cơ học
    KE = 0.5 .* damper.m .* v.^2          # Động năng
    PE = 0.5 .* damper.k .* x.^2          # Thế năng
    mechanical_energy = KE + PE
    
    # Năng lượng tiêu tán cơ học
    mechanical_dissipated = cumul_integral(t, damper.c .* v.^2)
    
    # Năng lượng điện tiêu tán
    electrical_power = (damper.R + damper.R_load) .* i.^2
    electrical_energy = cumul_integral(t, electrical_power)
    
    # Năng lượng thu được (trên tải)
    harvested_power = damper.R_load .* i.^2
    harvested_energy = cumul_integral(t, harvested_power)
    
    return (t=t, x=x, v=v, i=i,
            KE=KE, PE=PE, mech_energy=mechanical_energy,
            mech_diss=mechanical_dissipated,
            elec_energy=electrical_energy,
            harv_energy=harvested_energy,
            harv_power=harvested_power)
end

# Hàm tích phân số
function cumul_integral(t, y)
    result = zeros(length(y))
    for j in 2:length(y)
        dt = t[j] - t[j-1]
        result[j] = result[j-1] + 0.5 * (y[j] + y[j-1]) * dt
    end
    return result
end

# ------------------------------------------------------------
# 6. VẼ ĐỒ THỊ KẾT QUẢ
# ------------------------------------------------------------

function plot_results(results, damper, load_type)
    t = results.t
    
    p1 = plot(t, results.x, 
              title="Chuyển vị công trình",
              xlabel="Thời gian (s)", ylabel="Chuyển vị (m)",
              legend=false, linewidth=2)
    
    p2 = plot(t, results.v,
              title="Vận tốc công trình",
              xlabel="Thời gian (s)", ylabel="Vận tốc (m/s)",
              legend=false, linewidth=2)
    
    p3 = plot(t, results.i,
              title="Dòng điện cảm ứng",
              xlabel="Thời gian (s)", ylabel="Dòng điện (A)",
              legend=false, linewidth=2, color=:red)
    
    p4 = plot(t, results.harv_power,
              title="Công suất thu được",
              xlabel="Thời gian (s)", ylabel="Công suất (W)",
              legend=false, linewidth=2, color=:green)
    
    # Năng lượng tích lũy
    p5 = plot(t, results.mech_energy, label="Cơ năng",
              title="Năng lượng hệ thống",
              xlabel="Thời gian (s)", ylabel="Năng lượng (J)",
              linewidth=2)
    plot!(t, results.elec_energy, label="Điện năng tiêu tán", linewidth=2)
    plot!(t, results.harv_energy, label="Năng lượng thu được", linewidth=2)
    
    # Hiệu quả giảm dao động
    if load_type == :earthquake
        # Tính chỉ số giảm dao động
        x_max = maximum(abs.(results.x))
        plot!(title=@sprintf("Giảm dao động: x_max = %.4f m", x_max))
    end
    
    plot(p1, p2, p3, p4, p5, layout=(5,1), size=(800, 1200))
end

# ------------------------------------------------------------
# 7. VÍ DỤ CHẠY MÔ PHỎNG
# ------------------------------------------------------------

function run_example()
    println("="^60)
    println("MÔ PHỎNG HỆ CƠ-ĐIỆN TỪ: CÔNG TRÌNH VỚI DAMPER ĐIỆN TỪ")
    println("="^60)
    
    # 7.1 Trường hợp 1: Tải trọng gió
    println("\n1. MÔ PHỎNG VỚI TẢI TRỌNG GIÓ")
    sol_wind, damper_wind = simulate_emd_system(
        load_type=:wind,
        simulation_time=20.0
    )
    
    results_wind = calculate_energy(sol_wind, damper_wind)
    
    println("\nKẾT QUẢ - Tải gió:")
    println("  Chuyển vị cực đại: $(maximum(abs.(results_wind.x))) m")
    println("  Dòng điện cực đại: $(maximum(abs.(results_wind.i))) A")
    println("  Công suất cực đại: $(maximum(results_wind.harv_power)) W")
    println("  Năng lượng thu được: $(results_wind.harv_energy[end]) J")
    
    # 7.2 Trường hợp 2: Tải trọng động đất
    println("\n2. MÔ PHỎNG VỚI TẢI TRỌNG ĐỘNG ĐẤT")
    sol_eq, damper_eq = simulate_emd_system(
        load_type=:earthquake,
        simulation_time=30.0
    )
    
    results_eq = calculate_energy(sol_eq, damper_eq)
    
    println("\nKẾT QUẢ - Động đất:")
    println("  Chuyển vị cực đại: $(maximum(abs.(results_eq.x))) m")
    println("  Dòng điện cực đại: $(maximum(abs.(results_eq.i))) A")
    println("  Công suất cực đại: $(maximum(results_eq.harv_power)) W")
    println("  Năng lượng thu được: $(results_eq.harv_energy[end]) J")
    
    # 7.3 So sánh hiệu quả giảm dao động
    println("\n" * "="^60)
    println("SO SÁNH HIỆU QUẢ GIẢM DAO ĐỘNG")
    println("="^60)
    
    # Mô phỏng hệ không có damper điện từ (chỉ cản cơ học)
    damper_noEM = ElectromagneticDamper(
        10000.0, 10000.0, 1e6, 
        0.0, 0.0, 5.0, 0.1, 10.0, 0.7  # Ke=Kt=0: không có hiệu ứng điện từ
    )
    
    # Vẽ so sánh
    println("\nĐang vẽ đồ thị so sánh...")
    
    # Vẽ kết quả cho từng trường hợp
    plot_results(results_wind, damper_wind, :wind)
    savefig("wind_damper_results.png")
    
    plot_results(results_eq, damper_eq, :earthquake)
    savefig("earthquake_damper_results.png")
    
    println("\nĐã lưu đồ thị:")
    println("  - wind_damper_results.png")
    println("  - earthquake_damper_results.png")
    
    return (results_wind, results_eq)
end

# ------------------------------------------------------------
# 8. PHÂN TÍCH THAM SỐ TỐI ƯU
# ------------------------------------------------------------

function parameter_sensitivity_analysis()
    println("\n" * "="^60)
    println("PHÂN TÍCH ĐỘ NHẠY THAM SỐ")
    println("="^60)
    
    # Khảo sát ảnh hưởng của điện trở tải
    R_loads = [1.0, 5.0, 10.0, 20.0, 50.0]
    displacements = Float64[]
    harvested_energies = Float64[]
    
    for R_load in R_loads
        sol, damper = simulate_emd_system(
            R_load=R_load,
            load_type=:wind,
            simulation_time=10.0
        )
        
        results = calculate_energy(sol, damper)
        push!(displacements, maximum(abs.(results.x)))
        push!(harvested_energies, results.harv_energy[end])
        
        println(@sprintf("R_load = %5.1f Ω: x_max = %.4f m, Energy = %.2f J", 
                        R_load, displacements[end], harvested_energies[end]))
    end
    
    # Vẽ đồ thị sensitivity
    p1 = plot(R_loads, displacements,
              title="Ảnh hưởng của điện trở tải đến chuyển vị",
              xlabel="Điện trở tải (Ω)", ylabel="Chuyển vị cực đại (m)",
              marker=:circle, linewidth=2, label="")
    
    p2 = plot(R_loads, harvested_energies,
              title="Năng lượng thu được theo điện trở tải",
              xlabel="Điện trở tải (Ω)", ylabel="Năng lượng (J)",
              marker=:circle, linewidth=2, label="", color=:green)
    
    plot(p1, p2, layout=(2,1), size=(800, 600))
    savefig("parameter_sensitivity.png")
    
    println("\nĐã lưu: parameter_sensitivity.png")
end

# ------------------------------------------------------------
# 9. MÔ HÌNH NÂNG CAO: CÔNG TRÌNH MDOF
# ------------------------------------------------------------

module MDOF_EMD
using LinearAlgebra
using SparseArrays
using DifferentialEquations

export MDOFStructure, MultiDamperSystem, simulate_mdof

"""
Công trình nhiều bậc tự do với nhiều damper điện từ
"""
struct MDOFStructure
    M::Matrix{Float64}     # Ma trận khối lượng
    C::Matrix{Float64}     # Ma trận cản
    K::Matrix{Float64}     # Ma trận độ cứng
    ndof::Int              # Số bậc tự do
end

struct EMDevice
    floor::Int             # Tầng lắp đặt
    Ke::Float64           # Hệ số cảm ứng
    Kt::Float64           # Hệ số lực
    R::Float64            # Điện trở
    L::Float64            # Điện cảm
    R_load::Float64       # Điện trở tải
end

struct MultiDamperSystem
    structure::MDOFStructure
    devices::Vector{EMDevice}
    n_devices::Int
end

function mdof_equations!(du, u, p, t)
    system, F_external = p
    n = system.structure.ndof
    m = system.n_devices
    
    # u = [x1...xn, v1...vn, i1...im]
    x = u[1:n]
    v = u[n+1:2n]
    i = u[2n+1:2n+m]
    
    # Lực điện từ từ các damper
    F_em = zeros(n)
    for (idx, device) in enumerate(system.devices)
        floor = device.floor
        F_em[floor] += device.Kt * i[idx]
    end
    
    # Phương trình cơ học: M*dv/dt = F_ext - C*v - K*x - F_em
    du[1:n] = v
    du[n+1:2n] = system.structure.M \ (F_external(t) - 
                                       system.structure.C*v - 
                                       system.structure.K*x - 
                                       F_em)
    
    # Phương trình điện cho từng damper
    for (idx, device) in enumerate(system.devices)
        floor = device.floor
        du[2n+idx] = (device.Ke*v[floor] - 
                     (device.R + device.R_load)*i[idx]) / device.L
    end
end

function simulate_mdof()
    # Ví dụ: Công trình 5 tầng
    n = 5
    m_val = 10000.0  # Khối lượng mỗi tầng
    k_val = 1e6      # Độ cứng mỗi tầng
    
    M = m_val * Matrix(I, n, n)
    K = zeros(n, n)
    for i in 1:n
        K[i,i] = 2*k_val
        if i > 1
            K[i,i-1] = -k_val
            K[i-1,i] = -k_val
        end
    end
    
    # Ma trận cản Rayleigh: C = αM + βK
    α = 0.1
    β = 0.001
    C = α * M + β * K
    
    structure = MDOFStructure(M, C, K, n)
    
    # Thiết bị điện từ ở các tầng 3 và 5
    devices = [
        EMDevice(3, 50.0, 50.0, 5.0, 0.1, 10.0),
        EMDevice(5, 50.0, 50.0, 5.0, 0.1, 10.0)
    ]
    
    system = MultiDamperSystem(structure, devices, length(devices))
    
    # Điều kiện ban đầu
    u0 = zeros(2*n + length(devices))
    
    # Tải trọng gió (lực tập trung tầng mái)
    F_ext = t -> [0.0, 0.0, 0.0, 0.0, 10000*sin(2π*0.5*t)]
    
    # Giải hệ phương trình
    prob = ODEProblem(mdof_equations!, u0, (0.0, 20.0), (system, F_ext))
    sol = solve(prob, Tsit5())
    
    return sol, system
end

end  # module MDOF_EMD

# ------------------------------------------------------------
# 10. CHẠY CHƯƠNG TRÌNH CHÍNH
# ------------------------------------------------------------

function main()
    println("CHƯƠNG TRÌNH MÔ PHỎNG DAMPER ĐIỆN TỪ CHO CÔNG TRÌNH")
    println("="^60)
    
    # Chạy ví dụ cơ bản
    println("\n[PHẦN 1] MÔ PHỎNG HỆ SDOF")
    results = run_example()
    
    # Phân tích độ nhạy
    println("\n[PHẦN 2] PHÂN TÍCH ĐỘ NHẠY THAM SỐ")
    parameter_sensitivity_analysis()
    
    # Mô hình MDOF
    println("\n[PHẦN 3] MÔ HÌNH CÔNG TRÌNH NHIỀU TẦNG")
    using .MDOF_EMD
    sol_mdof, system_mdof = MDOF_EMD.simulate_mdof()
    
    # Trích xuất kết quả MDOF
    n = system_mdof.structure.ndof
    m = system_mdof.n_devices
    t_mdof = sol_mdof.t
    roof_disp = [sol_mdof[i][n] for i in 1:length(sol_mdof)]
    
    println("\nKết quả MDOF (5 tầng):")
    println("  Chuyển vị cực đại tầng mái: $(maximum(abs.(roof_disp))) m")
    
    # Vẽ kết quả MDOF
    plot(t_mdof, roof_disp,
         title="Chuyển vị tầng mái - Mô hình MDOF",
         xlabel="Thời gian (s)", ylabel="Chuyển vị (m)",
         linewidth=2)
    savefig("mdof_roof_displacement.png")
    
    println("\n" * "="^60)
    println("KẾT THÚC MÔ PHỎNG")
    println("="^60)
    println("\nCác file đã tạo:")
    println("  1. wind_damper_results.png")
    println("  2. earthquake_damper_results.png")
    println("  3. parameter_sensitivity.png")
    println("  4. mdof_roof_displacement.png")
    
    return results
end

# Chạy chương trình
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end