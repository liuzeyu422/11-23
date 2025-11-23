// ===== UEqn.H : falling-film pyrolysis with Wallis interfacial drag (stabilised) =====

// 1) 基础动量方程
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
 ==
    fvOptions(rho, U)
);

UEqn.relax();
fvOptions.constrain(UEqn);

// ------------------------- Wallis interfacial drag (with safety limiters) -------------------------

// 运行时可调参数（存在即覆盖）
IOdictionary tp
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// α 带宽（只在界面附近生效）
scalar aMin = tp.lookupOrDefault<scalar>("dragBandAlphaMin", 0.05);
scalar aMax = tp.lookupOrDefault<scalar>("dragBandAlphaMax", 0.95);

// 全局缩放
scalar Cdrag = tp.lookupOrDefault<scalar>("CdragWallis", 0.01);

// 板宽/等效周长 W
dimensionedScalar Wdim("plateWidth", dimLength, 2.0*constant::mathematical::pi*0.3);
if (tp.found("plateWidth")) Wdim = dimensionedScalar(tp.lookup("plateWidth"));

// 表面张力（两种气相按 α2/α3 加权）
dimensionedScalar sigma12("sigma12", dimForce/dimLength, 0.035);
dimensionedScalar sigma13("sigma13", dimForce/dimLength, 0.035);
if (tp.found("sigma12")) sigma12 = dimensionedScalar(tp.lookup("sigma12"));
if (tp.found("sigma13")) sigma13 = dimensionedScalar(tp.lookup("sigma13"));

// 数值护栏（可从字典覆盖）
dimensionedScalar tinyNu ("tinyNu",  dimViscosity, tp.lookupOrDefault<scalar>("dragTinyNu", 1e-9));     // m2/s
dimensionedScalar tinyL  ("tinyL",   dimLength,    tp.lookupOrDefault<scalar>("dragTinyDelta", 1e-6));  // m
dimensionedScalar Ucap   ("Ucap",    dimVelocity,  tp.lookupOrDefault<scalar>("dragUcap", 5.0));        // m/s
scalar            ReCap  = tp.lookupOrDefault<scalar>("dragReCap", 5e3);
scalar            Cgrad  = tp.lookupOrDefault<scalar>("dragGradAlphaCap", 5.0); // |∇α| ≤ Cgrad/δ

// 2) 面量与带宽
const surfaceScalarField alpha1f = fvc::interpolate(alpha1);
const surfaceScalarField band    = pos(alpha1f - aMin)*pos(aMax - alpha1f);

//const surfaceVectorField Uf   = fvc::interpolate(U);
const surfaceScalarField aGas = max(scalar(1) - alpha1f, scalar(SMALL));

// 建一个面场 dU，默认 0
surfaceVectorField dU
(
    IOobject("dUslip", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedVector("zero", dimVelocity, vector::zero)
);

// 内部面：用 owner/neighbour 两侧取样
const labelUList& own = mesh.owner();
const labelUList& nei = mesh.neighbour();

vectorField&       dUif = dU.primitiveFieldRef();     // internal faces
const scalarField& a1c  = alpha1.internalField();     // cell-centered alpha1
const vectorField& Uc   = U.internalField();          // cell-centered U

forAll(own, fI)
{
    // 显式屏蔽：非界面带宽内不组装滑移
    if (band.primitiveField()[fI] == 0)
    {
        dUif[fI] = vector::zero;
        continue;
    }
    
    const label c0 = own[fI];
    const label c1 = nei[fI];

    // 把 alpha1 较大的一侧当作“液侧”，较小的一侧当作“气侧”
    const label cL = (a1c[c0] >= a1c[c1]) ? c0 : c1;  // liquid-side cell
    const label cG = (cL == c0) ? c1 : c0;            // gas-side cell

    // 相对速度（气相 - 液膜表面）
    dUif[fI] = Uc[cG] - Uc[cL];
}

// 边界面：保持 0（曳力主要在内部界面起作用；band 亦会屏蔽非界面处）
forAll(dU.boundaryFieldRef(), patchi)
{
    dU.boundaryFieldRef()[patchi] = vector::zero;
}

// 3) 物性与几何尺度（面）
const surfaceScalarField rhoGf =
    fvc::interpolate(alpha2*rho2 + alpha3*rho3)/(aGas + scalar(SMALL));
const surfaceScalarField rhoLf =
    fvc::interpolate(alpha1*rho1)/(alpha1f + scalar(SMALL));

// 液相运动黏度（界面）
const surfaceScalarField nuLf = max(mixture->nu1f(), tinyNu);

// 局部膜厚 δ ← 到 heatedWall 的几何距离 y（createFields.H 已给出 y）
// —— 可选：膜厚夹持上下限（若未在字典给出，则取经验默认值）——
dimensionedScalar dMin("dragDeltaMin", dimLength, tp.lookupOrDefault<scalar>("dragDeltaMin", 1e-4));  // ~ 100 μm
dimensionedScalar dMax("dragDeltaMax", dimLength, tp.lookupOrDefault<scalar>("dragDeltaMax", 2e-3));  
// 局部膜厚 δ：以几何壁距为基准，并夹持到 [max(tinyL,dMin), dMax]
const surfaceScalarField deltaF =
    min( max(fvc::interpolate(y), max(tinyL, dMin)), dMax );

//const surfaceScalarField deltaF = max(fvc::interpolate(y), tinyL);

// 重力、σ 加权
const dimensionedScalar gMag = mag(g);
const surfaceScalarField sigmaLGf =
    fvc::interpolate(alpha2)*sigma12 + fvc::interpolate(alpha3)*sigma13;

// 毛细长度 L_l = sqrt( σ / ((ρ_l-ρ_g)|g|) )
const surfaceScalarField Ll =
    sqrt
    (
        max
        (
            sigmaLGf
           /( max(rhoLf - rhoGf,
                  dimensionedScalar("epsR", rhoLf.dimensions(), SMALL)
                 ) * gMag ),
            dimensionedScalar("epsL", dimLength*dimLength, SMALL)
        )
    );

// --- Wallis 组合量（加限幅以防极端值） ---
// Re_e = |dU| δ / ν_l
const surfaceScalarField dUMag    = mag(dU);
const dimensionedScalar tinyU("tinyU", dimVelocity, SMALL);
const surfaceScalarField dUMagLim = min(dUMag, Ucap);                            // |ΔU| ≤ Ucap
const surfaceVectorField dUdir    = dU/(dUMag + tinyU);                          // 单位向量
const surfaceVectorField dUsqLim  = dUdir * sqr(dUMagLim);                       // 方向 × |ΔU|^2

const surfaceScalarField ReeRaw   = dUMag*deltaF/max(nuLf, tinyNu);
const surfaceScalarField Ree      = min(max(ReeRaw, scalar(1)), scalar(ReCap));  // 1 ≤ Re_e ≤ ReCap

const surfaceScalarField deltaStar = max(deltaF/Ll, scalar(1e-6));
const surfaceScalarField Wstar     = max(Wdim/Ll,   scalar(1e-6));
const surfaceScalarField Nexp      = 3.95/(1.8 + 3.0/Wstar);

// f = 0.079 Re_e^{-0.25} (1 + 115 δ*^N)
const surfaceScalarField fcoeff =
    0.079*pow(Ree, -0.25)*(scalar(1) + 115.0*pow(deltaStar, Nexp));

// |∇α1|：再加一个基于 δ 的上限，|∇α1| ≤ Cgrad/δ
const surfaceScalarField gradAlphaRaw = mag(fvc::snGrad(alpha1));
const surfaceScalarField gradAlphaCap = min(gradAlphaRaw, scalar(Cgrad)/deltaF);

// （N/m^3），只在 α 带宽内生效，并做全局缩放 Cdrag
const surfaceVectorField FdragSf =
    Cdrag * band * fcoeff * rhoGf * gradAlphaCap * dUsqLim;



// 体力密度（N/m3）：把面量平均到单元中心（不做额外/乘法）
const volVectorField Fdrag
(
    IOobject("Fdrag", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    fvc::average(FdragSf)   // 仅做面到体的平均，不乘面积、不除体积
);

// 动量预测器
if (pimple.momentumPredictor())
{
    solve( UEqn == Fdrag );
    fvOptions.correct(U);
}

// ===== end UEqn.H =====

