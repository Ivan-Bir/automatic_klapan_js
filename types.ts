export interface MathModelStateCoeffs {
    pa_vh1?: number,
    pb_vh1?: number,
    pa_12?: number,
    pb_12?: number,
    pa_23?: number,
    pb_23?: number,
    pa_2a?: number,
    pb_2a?: number,
    pa_14?: number,
    pb_14?: number,
    pa_4v?: number,
    pb_4v?: number,
    pa_45?: number,
    pb_45?: number,
    pa_5a?: number,
    pb_5a?: number,
    Tvh1?: number,
    T12?: number,
    T23?: number,
    Ta2?: number,
    T45?: number,
    T14?: number,
    Tv4?: number,
    Ta5?: number,
    fkl?: number,
    g1vh?: number,
    gvh1?: number,
    g12?: number,
    g21?: number,
    g23?: number,
    g32?: number,
    g2a?: number,
    ga2?: number,
    g14?: number,
    g41?: number,
    g4v?: number,
    gv4?: number,
    g45?: number,
    g54?: number,
    g5a?: number,
    ga5?: number,
    v1?: number,
    v3?: number,
    v3p?: number,
    v1p?: number,
    s1?: number,
    s3?: number,
    f5?: number,
    f7?: number,
    f2?:number,
    f4?:number,

    f1:number,
    f3:number,
    f6:number,
    f8:number,
    fkp: number,
    m2:number,
    m4:number,
};

export interface MathModelCharacteristics {
    P1: number,
    T1: number,
    P2: number,
    T2: number,
    P3: number,
    T3: number,
    P4: number,
    T4: number,
    P5: number,
    T5: number,
    Pa: number,
    Ta: number,
    Xp: number,
    Vp: number,
    _Xp: number,
    _Vp: number,
};

export interface MathModelObject {
    functionalParams: MathModelCharacteristics,
    options: MathModelStateCoeffs,
};

export interface RK_factors {
    P1: number,
    T1: number,
    P2: number,
    T2: number,
    P3: number,
    T3: number,
    P4: number,
    T4: number,
    P5: number,
    T5: number,
    Pa: number,
    Ta: number,
    Xp: number,
    Vp: number,
    _Xp:number,
    _Vp: number,
}

export interface Solver {
    (params: MathModelCharacteristics, mutableOptions: MathModelStateCoeffs) : RK_factors
}

export interface State {
    isClosed: 0|1,
    isOpened: 0|1,
    isRabProcess: 0|1,
    fromOpenToClose: 0|1,
    fromCloseToOpen: 0|1,
    isDone: 0|1,
}

