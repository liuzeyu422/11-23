{
    // 1) 用当前速度一致地重建体积流量通量
    phi = fvc::flux(U);

    // 2) 让 phi 在边界上与 p_rgh 的边界条件一致（5.x 提供的工具函数）
    adjustPhi(phi, U, p_rgh);

    // 3) 连续性误差输出
    #include "continuityErrs.H"
}

