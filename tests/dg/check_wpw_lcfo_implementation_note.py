#!/usr/bin/env python3
"""Textual contract for the divided-SCF/WF+PW LCFO implementation note."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
NOTE = ROOT / "docs/notes/wpw_lcfo_divided_scf_implementation.tex"

assert NOTE.is_file(), f"missing implementation note: {NOTE}"
text = NOTE.read_text(encoding="utf-8")

required = {
    "LuaLaTeX Japanese font support": r"\usepackage{luatexja-fontspec}",
    "standard Japanese body font": r"\setmainjfont[BoldFont={HaranoAjiMincho-Bold}]{HaranoAjiMincho-Regular}",
    "standard Japanese heading font": r"\setsansjfont[BoldFont={HaranoAjiGothic-Bold}]{HaranoAjiGothic-Regular}",
    "sans serif section headings": r"\titleformat{\section}{\Large\bfseries\sffamily}",
    "route limitation upfront heading": r"\section*{最初に確認する計算ルートと適用限界}",
    "not fast iterative route": "高速Hamiltonian作用と反復的基底状態計算を用いる別ルートではない",
    "retained basis count upfront": r"M=\sum_f M_f^{\rm keep}",
    "row distributed dense storage": "各rankは自分が所有する行について全$M$列を保持する",
    "logical dense matrices": r"論理的にはdenseな$M\times M$行列",
    "quadratic global storage": r"O(M^2)",
    "cubic dense diagonalization": r"O(M^3)",
    "distribution scaling caveat": "分散化しても全体の漸近計算量は変わらない",
    "large scale limitation": "大規模系の最終基底状態計算には使用できない",
    "separate route scope exclusion": r"O(N\log N)",
    "upfront data flow basis ownership": "保持基底のglobal IDを一つのrankが一意に所有する",
    "upfront data flow LCFO rows": "自rankが所有するLCFO行",
    "upfront data flow ScaLAPACK": "二次元block-cyclic配置へ再分配する",
    "buffer restricted to immediate neighbors": "第二近接fragmentの格子点を参照することはない",
    "collective communication is not distant dependence": "全rankがcollective通信に参加することと、遠方fragmentの値に物理的に依存することは異なる",
    "DG WF PW SCF heading": r"\subsection{WF+PW基底へDG Hamiltonianを射影して密度を更新する}",
    "fragment eigensystem goal stated first": r"$H_{0,f}+H_f^{\rm DG}$のfragment固有値とfragment固有状態を求める",
    "fragment eigensystem goal equation": r"\left(H_{0,f}+H_f^{\rm DG}\right)c_{fn}=S_fc_{fn}\epsilon_{fn}",
    "complete DG weak form in SCF": "volume項、interface flux項、随伴flux項、penalty項を全て含む",
    "fragment DG state scope": "2.4で求めるのはfragment DG固有状態であり、全体系の固有状態ではない",
    "LCFO global coupling scope": "2.5で初めて、保持したfragment DG固有状態をLCFO基底として全体系を結合する",
    "fragment solve density purpose": "得られた占有fragment固有状態から次のSCF密度を作る",
    "legacy route invariant": "DC+LCFO/Wannier90",
    "single LCFO invariant": "exactly-once",
    "no density re-SCF invariant": "LCFO後の密度再SCFを行わない",
    "no density gate invariant": "LCFO後の密度収束判定を行わない",
    "Hamiltonian": r"\hat H",
    "metric": r"S_{",
    "window": r"w_f",
    "density": r"\rho",
    "generalized eigenproblem": r"HC=SC",
    "physical nearsightedness": "電子密度の近視性",
    "variational principle": "Rayleigh--Ritz変分原理",
    "window kinetic term": r"\nabla w_f",
    "boundary incompleteness": "fragment境界の基底不完全性",
    "inter-fragment hybridization": "fragment間hybridization",
    "no-rescf physical approximation": "凍結密度近似",
    "double counting explanation": "二重計数",
    "DG derivation heading": r"\section{DG領域分割の導出}",
    "explicit WF PW wavefunction expansion": r"\sum_{\alpha=1}^{M_f^{\rm WF}}\Phi_{f\alpha}(\boldsymbol r)a_{f\alpha n}^{(k)}",
    "implementation shorthand after derivation": r"X_f=[\,\Phi_f\;P_f\,]",
    "fragment combined basis definition": r"B_{f\mu}(\boldsymbol r)=\begin{cases}",
    "WF PW block Hamiltonian": r"H_f^{\mathrm{WF,WF}}",
    "general DG heading": r"\section{一般的なdiscontinuous Galerkin法}",
    "DG jump": r"\llbracket \psi\rrbracket",
    "DG average": r"\{\!\{\nabla \psi\}\!\}",
    "SIPG penalty": r"\frac{\eta_e}{h_e}",
    "numerical flux": "数値流束",
    "1D DG example": r"\Omega=(0,L)",
    "element partition": r"0=x_0<x_1<\cdots<x_{N_e}=L",
    "left right traces": r"\psi_h(x_j^-)",
    "H1 explanation": r"H^1(\Omega)",
    "DG motivation": "自由度を共有しない",
    "physicist wavefunction notation": r"\hat H\psi=E\psi",
    "trial wavefunction notation": r"\varphi^*",
    "kinetic interpretation": "運動エネルギー期待値",
    "probability amplitude explanation": "確率振幅",
    "piecewise DG ansatz": r"\left.\psi_h\right|_{K_j}",
    "element residual": r"\mathcal R_j",
    "strong substitution limitation": "強形式へ大域的に代入できない",
    "test function role": "試験波動関数は新しい物理状態ではない",
    "Galerkin test choice": r"\varphi=\phi_{j\beta}",
    "residual projection": r"\langle\varphi|\mathcal R\rangle=0",
    "DG Hamiltonian transition": r"\subsection{区分基底からDG Hamiltonianへ}",
    "independent coefficient reminder": "異なる要素に属する係数の間へ等値条件を課していない",
    "piecewise test function wording": "全ての区分的試験関数",
    "continuous x explanation": "$x$は区間内を連続に動く座標",
    "discrete j explanation": "$j$は要素と境界を数える離散的な番号",
    "same interface point": "別々の格子点ではない",
    "left limit notation": r"\lim_{x\to x_j^-}",
    "right limit notation": r"\lim_{x\to x_j^+}",
    "element residual projection": r"\int_{K_j}\phi_{j\beta}^*(x)\mathcal R_j(x)\,\dd x=0",
    "element boundary term": r"\left[\phi_{j\beta}^*\nabla\phi_{j\alpha}\right]",
    "boundary contribution definition": "境界$x_j$が弱形式へ与える総寄与を$B_j$と定義する",
    "left element endpoint": r"T_{j\to x_j}=-\frac12\varphi^{-*}(\nabla\psi_h)^-",
    "right element endpoint": r"T_{j+1\to x_j}=+\frac12\varphi^{+*}(\nabla\psi_h)^+",
    "boundary sum identity": r"B_j=T_{j\to x_j}+T_{j+1\to x_j}",
    "equation family explanation": "式(8)は一つの要素だけに対する一回限りの式ではない",
    "left equation instance": r"\ell=j",
    "right equation instance": r"\ell=j+1",
    "Kj integration interval": r"\int_{K_j}=\int_{x_{j-1}}^{x_j}",
    "Kj1 integration interval": r"K_{j+1}=(x_j,x_{j+1})",
    "endpoint sign explanation": "上端は負符号、下端は正符号",
    "neighbor exclusion": "$K_{j-1}$は境界$x_j$に接していない",
    "previous boundary accounting": "境界$x_{j-1}$を共有する二要素と後で組にする",
    "interface regrouping": r"\sum_{j=1}^{N_e-1}B_j",
    "global element sum": r"\sum_{\ell=1}^{N_e}\left[-\frac12\left[\varphi^*\nabla\psi_h\right]",
    "no double addition": "同じ項を二重に加える操作ではない",
    "external boundary definition": r"外側二端から来る総寄与を$B_{\rm ext}$と定義する",
    "external boundary formula": r"B_{\rm ext}=+\frac12\varphi^*(x_0^+)\nabla\psi_h(x_0^+)",
    "right external endpoint": r"-\frac12\varphi^*(x_{N_e}^-)\nabla\psi_h(x_{N_e}^-)",
    "local DG heading": r"\paragraph{$j$番目の領域のDG方程式を再定義する。}",
    "numerical derivative flux": r"\widehat{\nabla\psi}_{\,j}=\mathcal F",
    "current global DG equation": r"a_h^{\rm flux}(\varphi,\psi_h)=E_hm_h(\varphi,\psi_h)",
    "current global flux coupling": r"\left(\varphi^{-*}(x_j)-\varphi^{+*}(x_j)\right)",
    "provisional DG equation caveat": "SIPGを選ぶ前の暫定的なDG方程式",
    "global local-basis expansion": r"\psi_h(x)=\sum_{k=1}^{N_e}\sum_{\alpha=1}^{M_k}c_{k\alpha}\phi_{k\alpha}(x)",
    "same Galerkin basis": r"\varphi(x)=\phi_{j\beta}(x)",
    "same basis set explanation": "同じ基底集合から一つを選ぶ",
    "basis-substituted flux equation": r"\sum_{k=1}^{N_e}\sum_{\alpha=1}^{M_k}c_{k\alpha}",
    "provisional DG matrix element": r"H^{\rm flux}_{j\beta,k\alpha}=a_h^{\rm flux}",
    "provisional overlap element": r"S_{j\beta,k\alpha}=m_h",
    "1D gradient average": r"\{\nabla\psi_h\}_j=\frac12\left[(\nabla\psi_h)^-+(\nabla\psi_h)^+\right]",
    "1D wavefunction jump": r"[\psi_h]_j=\psi_h^- -\psi_h^+",
    "SIPG numerical gradient flux": r"\widehat{\nabla\psi}_{\,j}=\{\nabla\psi_h\}_j-\frac{\eta_j}{h_j}[\psi_h]_j",
    "central flux limitation": "平均だけではjumpを小さくする作用がない",
    "adjoint consistency caveat": "数値勾配fluxを定めるだけではHermitian性はまだ保証されない",
    "flux substitution boundary term": r"-\frac12[\varphi^*]_j\widehat{\nabla\psi}_{\,j}",
    "consistency penalty expansion": r"-\frac12[\varphi^*]_j\{\nabla\psi_h\}_j",
    "adjoint consistency term": r"-\frac12\{\nabla\varphi^*\}_j[\psi_h]_j",
    "1D SIPG derivation": r"a_{\rm kin}^{\rm SIPG}(\varphi,\psi_h)",
    "basis-expanded SIPG kinetic action": r"a_{\rm kin}^{\rm SIPG}(\phi_{p\beta},\psi_h)",
    "SIPG kinetic matrix": r"K^{\rm SIPG}_{p\beta,k\alpha}",
    "SIPG matrix consistency term": r"[\phi_{p\beta}^*]_j\{\nabla\phi_{k\alpha}\}_j",
    "SIPG matrix adjoint term": r"\{\nabla\phi_{p\beta}^*\}_j[\phi_{k\alpha}]_j",
    "SIPG matrix penalty term": r"[\phi_{p\beta}^*]_j[\phi_{k\alpha}]_j",
    "1D DG Hamiltonian matrix": r"H^{\rm DG}_{p\beta,k\alpha}=K^{\rm SIPG}_{p\beta,k\alpha}+V_{p\beta,k\alpha}",
    "1D DG generalized eigenproblem": r"H^{\rm DG}_{p\beta,k\alpha}c_{k\alpha}",
    "pointwise DG residual": r"\mathcal R_h(x)=\hat H\psi_h(x)-E_h\psi_h(x)",
    "projected DG residual": r"\langle\phi_{p\beta}|\mathcal R_h\rangle_{\rm DG}=0",
    "residual not pointwise zero": "点ごとに0であることを要求していない",
    "basis versus coefficient roles": "基底は近似空間を決め、係数はその空間内で射影残差を0にする",
    "orthogonality distinction": "基底同士の直交性と、残差の基底空間への直交性は別の条件",
    "orthonormal overlap consequence": r"S_{IJ}=\langle\phi_I|\phi_J\rangle=\delta_{IJ}",
    "complete basis residual limit": r"\mathcal R_h\longrightarrow0",
    "time dependent DG heading": r"\paragraph{時間依存DG方程式の導出。}",
    "time dependent Schrodinger equation": r"\mathrm i\,\partial_t\psi(x,t)=\hat H(t)\psi(x,t)",
    "time dependent basis coefficients": r"c_{k\alpha}(t)\phi_{k\alpha}(x)",
    "time derivative expansion": r"\partial_t\psi_h(x,t)=\sum_{k,\alpha}\dot c_{k\alpha}(t)\phi_{k\alpha}(x)",
    "time dependent projected DG equation": r"\mathrm i\,m_h(\varphi,\partial_t\psi_h)",
    "time dependent DG matrix equation": r"\mathrm i\,S\dot C(t)=H^{\rm DG}(t)C(t)",
    "fixed basis assumption": "局所基底は時間に依存しない",
    "moving basis connection term": r"D_{IJ}(t)=\langle\phi_I(t)|\partial_t\phi_J(t)\rangle",
    "DG norm conservation": r"\partial_t\!\left[C^\dagger(t)SC(t)\right]",
    "full-grid time-dependent KS equation": r"\mathrm i\,\partial_t\boldsymbol\psi_n(t)",
    "continuous spatial inner product": r"\int_\Omega u^*(\boldsymbol r)v(\boldsymbol r)\,\dd^3r",
    "spatial integral inner product": r"\langle u|v\rangle=\int_\Omega u^*(\boldsymbol r)v(\boldsymbol r)\,\dd^3r",
    "Hamiltonian bra-ket matrix element": r"\langle b_i|\hat H[\rho]|b_j\rangle",
    "buffer bra-ket definition": r"\langle u|v\rangle_{\Omega_f^{\rm buf}}",
    "full-grid time-dependent orthonormality": r"\int_\Omega\psi_m^*(\boldsymbol r,t)\psi_n(\boldsymbol r,t)\,\dd^3r=\delta_{mn}",
    "full-grid Hermiticity": r"\mathcal H^\dagger=\mathcal H",
    "LCFO to WF preparation heading": r"\subsection{既存DC+LCFO状態からWFを準備する}",
    "LCFO basis expansion": r"\psi_n^{(f)}(\boldsymbol r)=\sum_{\mu=1}^{M_f}",
    "LCFO coefficient role": r"C_{\mu n}^{(f)}$は展開係数",
    "Wannier transformation": r"&=\sum_{n=1}^{N_f^{\rm band}}\psi_n^{(f)}",
    "WF coefficient product": r"B_{\mu a}^{(f)}=\sum_n C_{\mu n}^{(f)}U_{na}^{(f)}",
    "C alone insufficient": r"$C^{(f)}$だけからWFを作るのではない",
    "symmetry transported WF": r"\hat g\,w_a=\sum_b w_bD_{ba}(g)",
    "Wannier gauge freedom": "位相、列の並び順、縮退部分空間内のユニタリ混合",
    "covariance residual": r"\Delta_g=\left\|\hat gW-WD(g)\right\|",
    "tolerance based symmetry guarantee": "許容値以下であることを検証する",
    "WF preparation order summary": r"\Phi^{(f)}C^{(f)}U^{(f)}",
    "occupied LCFO rank": r"N_{\rm occ}=\left\lceil N_e/2\right\rceil",
    "not all LCFO eigenstates": "LCFOで計算可能な全固有状態ではない",
    "complete sp complement": "complete $s+p$ projector channel",
    "Wannier target rank": r"N_{\rm target}=N_{\rm occ}+N_{s+p}",
    "projector construction heading": r"\paragraph{complete $s+p$ projectorをSALMON側で作る。}",
    "projector radial angular definition": r"p_{I,lm}(\boldsymbol r)=R_{I,l}(r_I)Y_{lm}^{\rm real}(\widehat{\boldsymbol d}_I)",
    "projector displacement definition": r"\boldsymbol d_I=\boldsymbol r-\boldsymbol R_I",
    "periodic nearest image": "最も近い周期像",
    "pseudo atomic orbital table": r"\code{pp\%upptbl\_ao}",
    "not nonlocal projector": "非局所擬ポテンシャル演算子のprojectorそのものではない",
    "four channels per atom": r"\{p_{I,s},p_{I,p_x},p_{I,p_y},p_{I,p_z}\}",
    "SALMON computed A matrix": r"A_{na}=\langle\chi_n|p_a\rangle",
    "projector overlap integral": r"\int_\Omega\chi_n^*(\boldsymbol r)p_a(\boldsymbol r)\,\dd^3r",
    "polar unitary gauge": r"A_{\rm seed}=U_AV_A^\dagger",
    "not Wannier internal projection": "Wannier90に同じprojectorを再生成させるのではない",
    "precomputed A matrix": r"\code{precomputed\_a\_matrix}",
    "symmetry adapted seed space": "閉じるように適応",
    "fragment LCFO contribution": r"\psi_n^{(f)}(\boldsymbol r)",
    "global LCFO assembly": r"\widetilde\psi_n^{\rm LCFO}(\boldsymbol r)=\sum_{f=1}^{N_f}w_f(\boldsymbol r)\psi_n^{(f)}(\boldsymbol r)",
    "partition weight definition": r"\sum_f w_f(\boldsymbol r)=1",
    "raw global seed definition": r"q_a(\boldsymbol r)=\begin{cases}",
    "orthonormal global band definition": r"\chi_n(\boldsymbol r)=\sum_{a=1}^{N_{\rm target}}q_a^{\rm sym}(\boldsymbol r)Z_{an}",
    "global band orthonormality": r"Z^\dagger S^{(q)}Z=I",
    "chi implementation mapping": r"\code{global\_closed\_core}",
    "chi is full cell": r"$\chi_n$は全セル上の関数",
    "LCFO bands already orthonormal": "LCFOで得た全系bandはすでに直交規格化されている",
    "ideal reconstruction identity": r"\widetilde\psi_n^{\rm LCFO}=\psi_n^{\rm LCFO}",
    "no separate chi necessity": r"占有LCFO blockだけをWannier90へ渡すなら、別の物理状態として$\chi_n$を導入する必要はない",
    "chi equals LCFO choice": r"\chi_n=\psi_n^{\rm LCFO}",
    "distributed reconstruction purpose": "分散全セル表現へmaterializeし直す",
    "orthogonalization is numerical safeguard": "新しい物理部分空間を作る操作ではない",
    "occupied span preserved": r"\operatorname{span}\{\chi_n^{\rm occ}\}=\operatorname{span}\{\psi_n^{\rm LCFO}\}",
    "combined chi distinction": r"projectorを加えた全入力空間の直交基底も同じ記号$\chi_n$で表している",
    "DMN filename": r"\code{overlapping\_wannier\_mlwf.dmn}",
    "DMN before Wannier90": "Wannier90を実行する前に",
    "band symmetry representation": r"D^{\rm band}(g)",
    "Wannier symmetry representation": r"D^{\rm WF}(g)",
    "symmetry intertwining condition": r"D^{\rm band}(g)U=UD^{\rm WF}(g)",
    "post Wannier covariance validation": "Wannier90が返した$U$について同じ共変性をもう一度検証",
    "Wannier interface heading": r"\paragraph{Wannier90へ渡す情報の構造。}",
    "no real space array handoff": r"実空間配列$\chi_n(\boldsymbol r)$そのものをWannier90へ渡すのではない",
    "M matrix definition": r"M_{mn}^{(b)}=\int_\Omega\chi_m^*(\boldsymbol r)",
    "M matrix phase": r"e^{-2\pi\mathrm i\boldsymbol n_b\cdot\boldsymbol\xi(\boldsymbol r)}",
    "M matrix shape": r"N_b\times N_b\times N_{\rm ntot}",
    "A matrix shape": r"N_b\times N_w",
    "Gamma one k point": r"N_k=1",
    "zero Wannier energies": r"\varepsilon_n^{\rm W90}=0",
    "rank local partial integrals": "各rankは自分が所有する実空間点だけから部分積分を計算",
    "rank zero Wannier call": "rank 0へreduceしてWannier90 libraryを呼ぶ",
    "Wannier return matrices": r"U^{\rm opt}\in\mathbb C^{N_b\times N_w}",
    "Wannier final transform": r"U^{\rm W90}=U^{\rm opt}U",
    "Wannier transform broadcast": "全MPI rankへbroadcast",
    "rank zero input list": "Wannier90呼出し直前にrank 0が保持する入力",
    "DMN file input distinction": r"\code{.dmn}は配列引数ではなく",
    "Wannier lwindow return": r"\code{lwindow}$(N_b,1)$",
    "Wannier total spread return": r"\Omega_{\rm total}",
    "wrapper internal U distinction": r"$U^{\rm opt}$と$U$はwrapper内部の一時配列",
    "wrapper exposed outputs": "wrapper外へ返して全rankへbroadcastするのは",
    "no Ng handoff explicit": r"$N_g$に比例する波動関数配列をWannier90へ渡さない",
    "global weak equation heading": "式(8)を全要素について足した全系の式",
    "notation heading": r"\section*{本書全体の記号規約}",
    "complex conjugate definition": "複素共役",
    "domain definition": r"\Omega\subset\mathbb R^d",
    "gradient definition": "空間勾配",
    "index ranges": "添字範囲",
    "physical point ownership": "physical-point ID",
    "fragment communicator": r"\texttt{icomm\_frag}",
    "total communicator": r"\texttt{icomm\_tot}",
    "implemented status": "実装済み",
    "unit-tested status": "単体試験済み",
    "pending physical validation": "物理検証保留",
    "MPI condition": "MPI 8 ranks",
    "OpenMP condition": r"OMP\_NUM\_THREADS=1",
    "no timeout condition": "時間打ち切りなし",
    "main implementation path": r"src/common/dg\_hybrid\_production\_pw\_basis.f90",
    "route implementation path": r"src/gs/main\_dft.f90",
    "LCFO implementation path": r"src/gs/dc/dg\_hybrid\_lcfo.f90",
    "one-shot implementation path": r"src/gs/dc/dg\_hybrid\_generalized\_eigensystem.f90",
    "checkpoint implementation path": r"src/rt/dg/rt\_dg\_hybrid\_checkpoint.f90",
    "test path": r"tests/dg/run\_dg\_hybrid\_si64\_divided\_lcfo.py",
}

missing = [label for label, token in required.items() if token not in text]
assert not missing, "implementation note contract missing: " + ", ".join(missing)
assert len(text.splitlines()) >= 60, "implementation note is still only a skeleton"
dg = text.index(r"\section{DG領域分割の導出}")
general_dg = text.index(r"\section{一般的なdiscontinuous Galerkin法}")
notation = text.index(r"\section*{本書全体の記号規約}")
route_limit = text.index(r"\section*{最初に確認する計算ルートと適用限界}")
scope = text.index(r"\section{適用範囲と絶対不変条件}")
ks = text.index(r"\section{連続Kohn--Sham問題から離散弱形式へ}")
assert route_limit < notation < general_dg < dg < scope and dg < ks, (
    "general DG, SALMON DG, and WF+PW derivations must precede implementation narrative"
)
general_body = text[general_dg:dg]
for premature in ("WF", "PW", "SALMON"):
    assert premature not in general_body, f"{premature} appears before the general DG derivation is complete"
for old_derivative in (r"\frac{\dd", r"\prime", "psi_h'", "varphi'", r"\nabla^2"):
    assert old_derivative not in general_body, f"mixed derivative notation remains: {old_derivative}"
assert "q_j" not in general_body and r"\widehat q" not in general_body, (
    "unnecessary q shorthand remains in the general DG derivation"
)
for obsolete_multidimensional_repetition in (
    r"\subsection{要素別部分積分と未決定な境界値}",
    r"\label{eq:element-parts}",
    r"\boldsymbol q=\nabla \psi",
    r"\widehat{\boldsymbol q}",
):
    assert obsolete_multidimensional_repetition not in general_body, (
        "obsolete repeated multidimensional derivation remains: "
        + obsolete_multidimensional_repetition
    )
flux_abstract = text.index(r"\widehat{\nabla\psi}_{\,j}=\mathcal F")
current_global_dg = text.index(r"a_h^{\rm flux}(\varphi,\psi_h)=E_hm_h(\varphi,\psi_h)")
basis_expansion = text.index(r"\psi_h(x)=\sum_{k=1}^{N_e}\sum_{\alpha=1}^{M_k}c_{k\alpha}\phi_{k\alpha}(x)")
same_basis_test = text.index(r"\varphi(x)=\phi_{j\beta}(x)", basis_expansion)
provisional_matrix = text.index(r"H^{\rm flux}_{j\beta,k\alpha}=a_h^{\rm flux}")
gradient_average = text.index(r"\{\nabla\psi_h\}_j=\frac12")
sipg_flux = text.index(r"\widehat{\nabla\psi}_{\,j}=\{\nabla\psi_h\}_j")
adjoint_caveat = text.index("数値勾配fluxを定めるだけではHermitian性はまだ保証されない")
assert flux_abstract < current_global_dg < basis_expansion < same_basis_test < provisional_matrix < gradient_average < sipg_flux < adjoint_caveat, (
    "the DG equation must be projected on the same basis set before SIPG is selected"
)
flux_substitution = text.index(r"-\frac12[\varphi^*]_j\widehat{\nabla\psi}_{\,j}")
adjoint_term = text.index(r"-\frac12\{\nabla\varphi^*\}_j[\psi_h]_j")
sipg_1d = text.index(r"a_{\rm kin}^{\rm SIPG}(\varphi,\psi_h)")
assert adjoint_caveat < flux_substitution < adjoint_term < sipg_1d, (
    "SIPG must be derived directly from numerical-flux substitution"
)
basis_sipg_action = text.index(r"a_{\rm kin}^{\rm SIPG}(\phi_{p\beta},\psi_h)")
sipg_matrix = text.index(r"K^{\rm SIPG}_{p\beta,k\alpha}", basis_sipg_action)
assert sipg_1d < basis_sipg_action < sipg_matrix, (
    "the local-basis expansion must be substituted immediately after the 1D SIPG form"
)
dg_eigenproblem = text.index(r"\label{eq:1d-dg-generalized-eigenproblem}")
pointwise_residual = text.index(r"\mathcal R_h(x)=\hat H\psi_h(x)-E_h\psi_h(x)")
projected_residual = text.index(r"\langle\phi_{p\beta}|\mathcal R_h\rangle_{\rm DG}=0")
assert dg_eigenproblem < pointwise_residual < projected_residual, (
    "pointwise and projected residuals must be explained after the DG eigenproblem"
)
td_dg = text.index(r"\paragraph{時間依存DG方程式の導出。}")
td_matrix = text.index(r"\mathrm i\,S\dot C(t)=H^{\rm DG}(t)C(t)")
assert projected_residual < td_dg < td_matrix, (
    "the time-dependent DG equation must follow the stationary residual discussion"
)
global_grid = text.index(r"\label{eq:global-grid}")
global_grid_td = text.index(r"\label{eq:global-grid-time-dependent}")
wf_preparation = text.index(r"\subsection{既存DC+LCFO状態からWFを準備する}")
wf_basis_expansion = text.index(r"\psi_n^{(f)}(\boldsymbol r)=\sum_{\mu=1}^{M_f}", wf_preparation)
wf_seed_selection = text.index(r"N_{\rm target}=N_{\rm occ}+N_{s+p}", wf_basis_expansion)
wf_projector = text.index(r"\paragraph{complete $s+p$ projectorをSALMON側で作る。}", wf_seed_selection)
wf_projector_definition = text.index(r"p_{I,lm}(\boldsymbol r)=R_{I,l}(r_I)Y_{lm}^{\rm real}", wf_projector)
wf_a_matrix = text.index(r"A_{na}=\langle\chi_n|p_a\rangle", wf_projector_definition)
wf_dmn = text.index(r"\code{overlapping\_wannier\_mlwf.dmn}", wf_a_matrix)
wf_transform = text.index(r"&=\sum_{n=1}^{N_f^{\rm band}}\psi_n^{(f)}", wf_dmn)
wf_symmetry = text.index(r"\hat g\,w_a=\sum_b w_bD_{ba}(g)", wf_dmn)
wf_covariance = text.index(r"\Delta_g=\left\|\hat gW-WD(g)\right\|", wf_symmetry)
local_trial_projection = text.index(r"\subsection{WF+PW基底へDG Hamiltonianを射影して密度を更新する}")
assert global_grid < global_grid_td < wf_preparation < wf_basis_expansion < wf_seed_selection < wf_projector < wf_projector_definition < wf_a_matrix < wf_dmn < wf_transform < wf_symmetry < wf_covariance < local_trial_projection, (
    "SALMON-side projector and A-matrix construction must precede DMN-constrained Wannier transformation and post-validation"
)
assert r"\subsection{WF+PWブロック行列}" not in text, (
    "the WF+PW block derivation must be integrated into the fragment SCF subsection"
)
for unimplemented_rt_claim in (
    "実時間計算の初期占有係数として保存する",
    "実時間計算では、このfragment局所固有値問題を時間stepごとに解き直さない",
    "保存された初期係数から大域WF+PW係数の時間発展方程式を積分する",
):
    assert unimplemented_rt_claim not in text, (
        "unimplemented final-LCFO-to-RT claim remains: " + unimplemented_rt_claim
    )
wfpw_explicit = text.index(r"\sum_{\alpha=1}^{M_f^{\rm WF}}\Phi_{f\alpha}(\boldsymbol r)a_{f\alpha n}^{(k)}", local_trial_projection)
wfpw_block = text.index(r"\label{eq:wfpw-block}", wfpw_explicit)
wfpw_shorthand = text.index(r"X_f=[\,\Phi_f\;P_f\,]", wfpw_block)
assert local_trial_projection < wfpw_explicit < wfpw_block < wfpw_shorthand, (
    "WF and PW must be derived explicitly before introducing the X_f implementation shorthand"
)
next_fragment_section = text.index(r"\subsection{fragment解から全体LCFOへの変換}", local_trial_projection)
fragment_dg_body = text[local_trial_projection:next_fragment_section]
assert r"a^{\rm DG}" not in fragment_dg_body, (
    "the concrete WF+PW calculation must use Hamiltonian matrix elements directly"
)
for obsolete_global_fragment_claim in (
    "fragmentごとの独立な固有値問題には分解できない",
    r"H_{f\mu,g\nu}^{\rm DG}",
    "全系の状態$n$を",
):
    assert obsolete_global_fragment_claim not in fragment_dg_body, (
        "section 2.4 still describes the global DG eigensystem: "
        + obsolete_global_fragment_claim
    )
assert r"\paragraph{実際の計算手順。}" not in text, (
    "obsolete generic halo/core five-step procedure remains before the WF+PW projection"
)
assert r"\label{eq:core-buffer}" not in text, (
    "the core/buffer explanation must be procedural rather than a set equation"
)
for obsolete_core_buffer_abstraction in (
    r"R_f^\dagger",
    r"\sum_fD_f=I",
    r"\sum_f W_f=I",
    r"\label{eq:core-pou}",
    r"\label{eq:window-pou}",
    r"\label{eq:buffer-exactness}",
):
    assert obsolete_core_buffer_abstraction not in text, (
        "obsolete core/buffer abstraction remains: " + obsolete_core_buffer_abstraction
    )
for obsolete_quadrature_matrix in (
    r"Q_f",
    r"Q=\operatorname{diag}",
    r"q_g",
    "求積行列",
    r"\Delta V",
):
    assert obsolete_quadrature_matrix not in text, (
        "quadrature-matrix notation remains despite the uniform Cartesian mesh: "
        + obsolete_quadrature_matrix
    )
for obsolete_spatial_grid_sum in (
    r"\sum_{g=1}^{N_g}u_g^*v_g",
    r"\sum_gb_i",
    r"\sum_{g\text{ on buffer }f}",
    r"\sum_g\rho_g",
    r"\sum_{g=1}^{N_g}\Xi",
):
    assert obsolete_spatial_grid_sum not in text, (
        "spatial grid sum remains where a spatial integral is required: "
        + obsolete_spatial_grid_sum
    )
assert "後の式\\eqref{eq:sipg}は、この三項をまとめて" not in general_body, (
    "SIPG derivation is still deferred instead of connected to the numerical flux"
)
for obsolete_local_flux_step in (
    r"\label{eq:local-dg-action}",
    r"\label{eq:local-dg-equation}",
    r"規則$\mathcal F$の具体形をまだ選んでいない",
):
    assert obsolete_local_flux_step not in general_body, (
        "obsolete intermediate local-flux equation remains: " + obsolete_local_flux_step
    )
for obsolete_broken_equation in (
    r"\label{eq:two-element-broken}",
    r"\label{eq:broken-space}",
    r"\bigoplus_{K\in\mathcal T_h}V_h(K)",
    r"V_h^{\mathrm{br}}",
):
    assert obsolete_broken_equation not in general_body, (
        "obsolete broken-space equation remains: " + obsolete_broken_equation
    )
ansatz = text.index(r"\left.\psi_h\right|_{K_j}")
weak = text.index(r"\label{eq:continuous-weak-first}")
assert ansatz < weak, "the piecewise wavefunction ansatz must precede the weak-form derivation"
mesh_label = text.index(r"\label{eq:1d-mesh}")
mesh_reference = text.index(r"\eqref{eq:1d-mesh}")
assert mesh_label < mesh_reference, "the 1D mesh equation must be defined before it is referenced"
element_weak = text.index(r"\label{eq:element-weak-with-boundary}")
global_weak = text.index(r"\label{eq:global-weak-before-regrouping}")
left_endpoint = text.index(r"\label{eq:left-element-interface-piece}")
assert element_weak < global_weak < left_endpoint, (
    "the global element sum must follow the one-element equation before interface expansion"
)
b_definition = text.index("境界$x_j$が弱形式へ与える総寄与を$B_j$と定義する")
b_sum = text.index(r"\sum_{j=1}^{N_e-1}B_j")
assert b_definition < b_sum, "B_j must be defined before the compact interface sum uses it"

print("WF+PW LCFO implementation note contract: PASS")
